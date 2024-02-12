from typing import Any, Mapping, List, Union
import re
import os
import json
import operator
from functools import reduce

class VariantAllele():
    def __init__(self, allele_index: str, alt: str, variant:dict) -> None:
        
        self.name = f"{variant.chromosome}:{variant.position}:{variant.ref}:{alt}"
        self.variant = variant
        self.allele_index = allele_index
        self.alt = alt
        self.type = "VariantAllele"
        self.allele_sequence = alt
        self.reference_sequence = variant.ref
        self.population_map = []
        self.info_map = self.traverse_csq_info()

    def get_allele_type(self):
        #TODO: change this to VariantAllele level
        return self.variant.get_allele_type(self.alt)

    def get_alternative_names(self):
        return self.variant.get_alternative_names()

    def get_slice(self):
        #TODO: review this to change to VariantAllele level
        return self.variant.get_slice(self.alt)
    
    def get_phenotype_assertions(self):
        min_alt = self.minimise_allele(self.alt)
        return self.info_map[min_alt]["phenotype_assertions"] if min_alt in self.info_map else []

    def get_predicted_molecular_consequences(self):
        ## compute info_map only for these methods, first check if non-empty and compute only then
        min_alt = self.minimise_allele(self.alt)
        return self.info_map[min_alt]["predicted_molecular_consequences"] if min_alt in self.info_map else []
    
    def get_prediction_results(self):
        min_alt = self.minimise_allele(self.alt)
        return self.info_map[min_alt]["prediction_results"] if min_alt in self.info_map else []
    
    def get_population_allele_frequencies(self):
        min_alt = self.minimise_allele(self.alt)
        population_map = self.variant.set_frequency_flags()
        return population_map[min_alt].values() if min_alt in population_map else []


    def traverse_csq_info(self) -> Mapping:
        """
        This function is to traverse the CSQ record and extract columns
        corresponding to Consequence, SIFT, PolyPhen, CADD
        """
        column_list = ["Allele", "PHENOTYPES", "Feature_type", "Feature", "Consequence", "SIFT", "PolyPhen", "SPDI", "CADD_PHRED","Conservation"]
        prediction_index_map = {}
        for col in column_list:
                if self.variant.get_info_key_index(col) is not None:
                    prediction_index_map[col.lower()] = self.variant.get_info_key_index(col) 

        info_map = {}
        for csq_record in self.variant.info["CSQ"]:
            csq_record_list = csq_record.split("|")
            allele = csq_record_list[prediction_index_map["allele"]]

            if allele not in info_map.keys():
                info_map[allele] = {"phenotype_assertions": [], "predicted_molecular_consequences": [], "prediction_results": []} 

            # parse and form phenotypes
            if not info_map[allele]["phenotype_assertions"]:
                phenotypes = csq_record_list[prediction_index_map["phenotypes"]].split("&") if "phenotypes" in prediction_index_map.keys() else []   
                for phenotype in phenotypes:
                    phenotype_assertions = self.create_allele_phenotype_assertion(phenotype) if phenotype else []
                    if (phenotype_assertions):
                        info_map[allele]["phenotype_assertions"].append(phenotype_assertions)
            
            # parse and form molecular consequences
            predicted_molecular_consequences = self.create_allele_predicted_molecular_consequence(csq_record_list, prediction_index_map)
            if (predicted_molecular_consequences):
                info_map[allele]["predicted_molecular_consequences"].append(predicted_molecular_consequences)
            
            # parse and form prediction results
            current_prediction_results = info_map[allele]["prediction_results"] 
            info_map[allele]["prediction_results"] += self.create_allele_prediction_results(current_prediction_results, csq_record_list, prediction_index_map)
        return info_map
     
    def create_allele_prediction_results(self, current_prediction_results: Mapping, csq_record: List, prediction_index_map: dict) -> list:
        prediction_results = []
        if "cadd_phred" in prediction_index_map.keys():
            if not self.prediction_result_already_exists(current_prediction_results, "CADD"):
                cadd_prediction_result = {
                        "score": csq_record[prediction_index_map["cadd_phred"]] ,
                        "analysis_method": {
                            "tool": "CADD",
                            "qualifier": "CADD"
                        }

                } if csq_record[prediction_index_map["cadd_phred"]] else None
                if cadd_prediction_result:
                    prediction_results.append(cadd_prediction_result)

        
        return prediction_results
    
    def prediction_result_already_exists(self, current_prediction_results: Mapping, tool: str) -> bool:
        for prediction_result in current_prediction_results:
            if prediction_result["analysis_method"]["tool"] == tool:
                return True

        return False
    
    def create_allele_predicted_molecular_consequence(self, csq_record: List, prediction_index_map: dict) -> Mapping:
        consequences_list = []
        if "consequence" in prediction_index_map.keys():
            for cons in csq_record[prediction_index_map["consequence"]].split("&"):
                consequences_list.append(
                    {
                        "accession_id": cons
                    }
                )

        prediction_results = []
        if "sift" in prediction_index_map.keys():
            sift_score = csq_record[prediction_index_map["sift"]]
            if sift_score:
                (result, score) = self.format_sift_polyphen_output(sift_score)
                if result is not None and score is not None:
                    sift_prediction_result = {
                            "result": result,
                            "score": score,
                            "analysis_method": {
                                "tool": "SIFT",
                                "qualifier": "SIFT"
                            }
                        }
                    prediction_results.append(sift_prediction_result)

        if "polyphen" in prediction_index_map.keys():
            polyphen_score = csq_record[prediction_index_map["polyphen"]]
            if polyphen_score:
                (result, score) = self.format_sift_polyphen_output(polyphen_score)
                if result is not None and score is not None:
                    polyphen_prediction_result = {
                            "result": result,
                            "score": score,
                            "analysis_method": {
                                "tool": "PolyPhen",
                                "qualifier": "PolyPhen"
                            }
                        }
                    prediction_results.append(polyphen_prediction_result)
        if consequences_list:
            return {
                "allele_name": csq_record[prediction_index_map["allele"]],
                "feature_stable_id": csq_record[prediction_index_map["feature"]],
                "feature_type": {
                    "accession_id": csq_record[prediction_index_map["feature_type"]]               
                } ,
                "consequences": consequences_list,
                "prediction_results": prediction_results
            }

    
    def format_sift_polyphen_output(self, output: str) -> tuple:
        try:
            (result, score) = re.split(r"[()]", output)[:2]
        except:
            return (None, None)

        if result not in [
            'probably damaging',
            'possibly damaging',
            'benign',
            'unknown',
            'tolerated',
            'deleterious',
            'tolerated - low confidence',
            'deleterious - low confidence',
        ]:
            result = None

        try:
            score = float(score)
        except:
            # need to log something here
            score = None

        return (result, score)
    
    def create_allele_phenotype_assertion(self, phenotype: str) -> Mapping:
        splits = phenotype.split("+")
        if len(splits) != 3 or splits[2].startswith("ENS"):
            return None

        phenotype_name,source,feature = splits
        evidence_list = []
        if re.search("^ENS.*G\d+", feature):
            feature_type = "Gene"
        elif re.search("^ENS.*T\d+", feature):
            feature_type = "Transcript"
        elif re.search("^rs", feature):
            feature_type = "Variant"  
        else:
            feature_type = None
        
        if phenotype:
            return {
                "feature": feature,
                "feature_type": {
                    "accession_id": feature_type                
                } ,
                "phenotype": {
                    "name": phenotype_name
                },
                "evidence": evidence_list
            }

    
    def minimise_allele(self, alt: str) -> str:
        """
        VCF file has the representation without anchoring bases
        for prediction scores in INFO column. This function is useful
        in matching the SPDI format in VCF with the allele in memory
        """
        minimised_allele_string = alt
        if len(alt) > len(self.variant.ref):
            minimised_allele_string = alt[1:] 
        elif len(alt) < len(self.variant.ref):
            minimised_allele_string = "-"
        return minimised_allele_string
    