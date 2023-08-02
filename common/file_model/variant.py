from typing import Any, Mapping, List, Union
import re
import operator
from functools import reduce

def reduce_allele(acc, alt):
    return (len(acc.value) < 2) & (len(alt.value) < 2)    

def reduce_allele_length(acc, alt):
    return len(acc.value) == len(alt.value)

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
    
    def get_primary_source(self) -> Mapping:
        source = self.header.get_lines("source")[0].value
        if re.search("^dbSNP", source):
            source_id = "dbSNP"
            source_name = "dnSNP"
            source_description = "NCBI db of human variants"
            source_url = "https://www.ncbi.nlm.nih.gov/snp/"
            source_release =154
            
        elif re.search("^ClinVar", source):
            source_id = "ClinVar"
            source_name = "ClinVar"
            source_description = "ClinVar db of human variants"
            source_url = "https://www.ncbi.nlm.nih.gov/clinvar/variation/"
            source_release = ""

        else:
            source_id = "test"
            source_name = "test"
            source_description = "test db of human variants"
            source_url = "https://www.ncbi.nlm.nih.gov/test/variation/"
            source_release = ""
        ### TODO: This is not elegant, need to serialise?
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
        if isinstance(allele, str) :
            allele_type, SO_term = self.set_allele_type(len(allele)<2, len(self.ref)<2, len(allele) == len(self.ref))
        elif isinstance(allele, list):
            allele_type, SO_term = self.set_allele_type(reduce(reduce_allele,allele), len(self.ref)<2, reduce(reduce_allele_length,allele) == len(self.ref))
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
    
    
    
    def get_slice(self) -> Mapping :
        allele_type = self.get_allele_type(self.alts)
        start = self.position
        length = len(self.ref)
        end = start + length -1
        if allele_type == "insertion":
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
        for alt in self.alts:
            variant_allele = self.create_variant_allele(alt.value)
            variant_allele_list.append(variant_allele)
        return variant_allele_list

    def create_variant_allele(self, alt: str) -> Mapping:
        name = f"{self.chromosome}:{self.position}:{self.ref}:{alt}"
        return {
            "name": name,
            "allele_sequence": alt,
            "reference_sequence": self.ref,
            "type": "VariantAllele",
            "allele_type": self.get_allele_type(alt),
            "slice": self.get_slice()
        }
        

