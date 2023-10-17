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
    
    def get_alternative_names(self) -> List:
        return []
    
    def get_primary_source(self) -> Mapping:
        try:
            source = self.header.get_lines("source")[0].value
            if re.search("^dbSNP", source):
                source_id = "dbSNP"
                source_name = "dbSNP"
                source_description = "NCBI db of human variants"
                source_url = "https://www.ncbi.nlm.nih.gov/snp/"
                source_release =154

            elif re.search("^ClinVar", source):
                source_id = "ClinVar"
                source_name = "ClinVar"
                source_description = "ClinVar db of human variants"
                source_url = "https://www.ncbi.nlm.nih.gov/clinvar/variation/"
                source_release = ""

        except:
            source_id = "test"
            source_name = "test"
            source_description = "test db of human variants"
            source_url = "https://www.ncbi.nlm.nih.gov/test/variation/"
            source_release = ""

        return {
            "accession_id": self.name,
            "name": self.name,
            "description": "",
            "assignment_method": {
                                "type": "DIRECT",
                                "description": "A reference made by an external resource of annotation to an Ensembl feature that Ensembl imports without modification"
                            },
            "url": f"{source_url}{self.name}",
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
        # self.set_frequency_flags(variant_allele_list)
        return variant_allele_list

    def create_variant_allele(self, alt: str) -> Mapping:
        """
        The function currently does not include parent-child resolvers.
        Ideally, each field should be a function called upon demand and does lazy loading
        This needs to be rewritten similar to Variant
        """
        name = f"{self.chromosome}:{self.position}:{self.ref}:{alt}"
        # min_alt = self.minimise_allele(alt)
        return {
            "name": name,
            "allele_sequence": alt,
            "reference_sequence": self.ref,
            "type": "VariantAllele"
            # "alternative_names": self.get_alternative_names(),
            # "allele_type": self.get_allele_type(alt),
            # "slice": self.get_slice(alt),
            # "phenotype_assertions": info_map[min_alt]["phenotype_assertions"] if min_alt in info_map else [],
            # "predicted_molecular_consequences": info_map[min_alt]["predicted_molecular_consequences"] if min_alt in info_map else [],
            # "prediction_results": info_map[min_alt]["prediction_results"] if min_alt in info_map else [],
            # "population_frequencies": self.get_population_allele_frequencies(frequency_map, allele_index)
        }

    
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

    
    def get_info_key_index(self, key: str, info_id: str ="CSQ") -> int:
            info_field = self.header.get_info_field_info(info_id).description
            csq_list = info_field.split("Format: ")[1].split("|")
            for index, value in enumerate(csq_list):
                if value == key:
                    return index   
    

