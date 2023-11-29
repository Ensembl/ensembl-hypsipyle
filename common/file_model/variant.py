from typing import Any, Mapping, List, Union
import re
import os
import json
import operator
from functools import reduce
from common.file_model.variant_allele import VariantAllele

def reduce_allele_length(allele_list: List):
    allele_length = -1
    for allele in allele_list:
        if len(allele.value) > allele_length:
            allele_length = len(allele.value)
    return allele_length



class Variant ():
    def __init__(self, record: Any, header: Any) -> None:
        self.name = record.ID[0]
        self.record = record 
        self.header = header
        self.chromosome = record.CHROM         ###TODO: convert the contig name in the file to match the chromosome id given in the payload 
        self.position = record.POS
        self.alts = record.ALT
        self.ref = record.REF
        self.info = record.INFO
        self.type = "Variant"
        self.vep_version = re.search("v\d+", self.header.get_lines("VEP")[0].value).group()
    
    def get_alternative_names(self) -> List:
        return []
    
    def get_primary_source(self) -> Mapping:
        """
        Fetches source from variant INFO columns
        Fallsback to fetching from the header
        """
        
        try:
            if "SOURCE" in self.info:
                source = self.info["SOURCE"]
            else:
                source = self.header.get_lines("source")[0].value

            if re.search("^dbSNP", source):
                source_id = "dbSNP"
                source_name = "dbSNP"
                source_description = "NCBI db of human variants"
                source_url = "https://www.ncbi.nlm.nih.gov/snp/"
                source_url_id = source_url
                source_release = 156
                variant_id = self.name
            
            elif re.search("^EVA", source):
                source_id = "EVA"
                source_name = "EVA"
                source_description = "European Variation Archive"
                source_url = "https://www.ebi.ac.uk/eva"
                source_url_id = "https://www.ebi.ac.uk/eva/?variant&accessionID="
                source_release = "release_5"
                variant_id = self.name

            elif re.search("^Ensembl", source):
                source_id = "Ensembl"
                source_name = "Ensembl"
                source_description = "Ensembl"
                source_url = "https://beta.ensembl.org"
                source_url_id = "https://beta.ensembl.org/"
                source_release = "110" # to be fetched from the file
                variant_id = f"{self.chromosome}:{self.position}:{self.name}"
            
            

        except Exception as e:
            return None 

        return {
            "accession_id": self.name,
            "name": self.name,
            "description": f"{source_description}",
            "assignment_method": {
                                "type": "DIRECT",
                                "description": "A reference made by an external resource of annotation to an Ensembl feature that Ensembl imports without modification"
                            },
            "url": f"{source_url_id}{variant_id}",
            "source": {

                        "id" : f"{source_id}",
                        "name": f"{source_name}",
                        "description": f"{source_description}",
                        "url":  f"{source_url}",
                        "release": f"{source_release}"
                        }
        }


        

    def set_allele_type(self, alt_one_bp: bool, ref_one_bp: bool, ref_alt_equal_bp: bool):         
        match [alt_one_bp, ref_one_bp, ref_alt_equal_bp]:
            case [True, True, True]: 
                allele_type = "SNV" 
                SO_term = "SO:0001483"

            case [True, False, False]: 
                allele_type = "deletion" 
                SO_term = "SO:0000159"

            case [False, True, False]: 
                allele_type = "insertion"
                SO_term = "SO:0000667"

            case [False, False, False]: 
                allele_type = "indel"
                SO_term = "SO:1000032"

            case [False, False, True]: 
                allele_type = "substitution"
                SO_term = "SO:1000002"   
        return allele_type, SO_term
    
    def get_allele_type(self, allele: Union[str, List]) -> Mapping :
        if isinstance(allele, str):
            if allele == self.ref:
                allele_type, SO_term = "biological_region","SO:0001411"
            else:
                allele_type, SO_term = self.set_allele_type(len(allele)<2, len(self.ref)<2, len(allele) == len(self.ref))
        elif isinstance(allele, list):
            alt_length = reduce_allele_length(allele)
            allele_type, SO_term = self.set_allele_type(alt_length < 2 , len(self.ref)<2, alt_length == len(self.ref))

        return {
            "accession_id": allele_type,
            "value": allele_type,
            "url": f"http://sequenceontology.org/browser/current_release/term/{SO_term}",
            "source": {
                    "id": "",
                    "name": "Sequence Ontology",
                    "url": "www.sequenceontology.org",
                    "description": "The Sequence Ontology..."
                    }

        }  
    
    def get_slice(self, allele: Union[str, List] ) -> Mapping :

        start = self.position
        length = len(self.ref)
        end = start + length -1
        if allele != self.ref:
            allele_type = self.get_allele_type(allele)
            if allele_type["accession_id"] == "insertion":
                length = 0
                end = start + 1
        
        return {
            "location": {
                "start": start,
                "end": end,
                "length": length
            },
            "region": {
                "name": self.chromosome,
                "code": "chromosome",
                "topology": "linear",
                "so_term": "SO:0001217"
            },
            "strand": {
                "code": "forward",
                "value": 1
            }
        }
    
    
    
    def get_alleles(self) -> List:
        variant_allele_list = []

        
        for index,alt in enumerate(self.alts):
            if index+1 <= len(self.alts):
                variant_allele = VariantAllele(index+1,  alt.value, self)
                variant_allele_list.append(variant_allele)
        reference_allele = VariantAllele(0, self.ref, self)
        variant_allele_list.append(reference_allele)
        return variant_allele_list
    
    def get_most_severe_consequence(self) -> Mapping:
        consequence_index = self.get_info_key_index("Consequence")
        consequence_map = {}
        directory = os.path.dirname(__file__)
        with open(os.path.join(directory,'variation_consequence_rank.json')) as rank_file:
            consequence_rank = json.load(rank_file)
        for csq_record in self.info["CSQ"]:
            csq_record_list = csq_record.split("|")
            for cons in csq_record_list[consequence_index].split("&"):
                rank = consequence_rank[cons]
                consequence_map[rank] = cons
        return{
                    "result": consequence_map[min(consequence_map.keys())]  ,
                    "analysis_method": {
                        "tool": "Ensembl VEP",
                        "qualifier": "most severe consequence"
                    }
        } 

    def get_gerp_score(self) -> Mapping:
        csq_record = self.info["CSQ"]
        csq_record_list = csq_record[0].split("|")
        if self.get_info_key_index("Conservation") is not None:
            gerp_index = self.get_info_key_index("Conservation") 
            gerp_prediction_result = {
                    "result": csq_record_list[gerp_index]  ,
                    "analysis_method": {
                        "tool": "GERP",
                        "qualifier": "GERP"
                    }
                } if csq_record_list[gerp_index] else {}
            return gerp_prediction_result
    
    def get_ancestral_allele(self) -> Mapping:
        csq_record = self.info["CSQ"]
        csq_record_list = csq_record[0].split("|")
        if self.get_info_key_index("AA") is not None:
            aa_index = self.get_info_key_index("AA")
            aa_prediction_result = {
                    "result": csq_record_list[aa_index] ,
                    "analysis_method": {
                        "tool": "AncestralAllele",
                        "qualifier": "",
                        "version": "110" #self.vep_version
                    }
                } if csq_record_list[aa_index] and csq_record_list[aa_index]!="."  else {}
            return aa_prediction_result
    
    def get_info_key_index(self, key: str, info_id: str ="CSQ") -> int:
            info_field = self.header.get_info_field_info(info_id).description
            csq_list = info_field.split("Format: ")[1].split("|")
            for index, value in enumerate(csq_list):
                if value == key:
                    return index   
    
    def set_frequency_flags(self, allele_list: List):
        """
        Calculates minor allele frequency by iterating through each allele 
        Assumption: Considers that population is only gnomAD (genomes) for now
        Sets the maf as hpmaf as gnomAD is the only population at the moment
        """
        maf_frequency = -1
        maf_index = -1
        highest_frequency = -1
        highest_frequency_index = -1
        highest_maf_frequency = -1
        highest_maf_frequency_index = -1
        maf_map = {}
        for allele_index, allele in enumerate(allele_list):
            if(len(allele["population_frequencies"]) > 0):
                pop = allele["population_frequencies"][0]
                pop_allele_frequency = float(pop["allele_frequency"])
                if ( pop_allele_frequency > maf_frequency and pop_allele_frequency < highest_frequency ):
                    maf_frequency = pop_allele_frequency
                    maf_index = allele_index
                elif ( pop_allele_frequency > highest_frequency ):
                    maf_frequency = highest_frequency
                    maf_index = highest_frequency_index
                    highest_frequency = pop_allele_frequency
                    highest_frequency_index = allele_index

        if maf_frequency>=0:
            allele_list[maf_index]["population_frequencies"][0]["is_minor_allele"]  = True
            allele_list[maf_index]["population_frequencies"][0]["is_hpmaf"]  = True  

