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
        self.info_map = self.traverse_csq_info()

    def get_allele_type(self):
        return self.variant.get_allele_type(self.alt)

    def get_alternative_names(self):
        return self.variant.get_alternative_names()

    def get_slice(self):
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
    
    
            
    def traverse_csq_info(self) -> Mapping:
        """
        This function is to traverse the CSQ record and extract columns
        corresponding to Consequence, SIFT, PolyPhen, CADD
        """
        allele_index = self.variant.get_info_key_index("Allele")
        phenotypes_index = self.variant.get_info_key_index("VAR_SYNONYMS")
        feature_type_index = self.variant.get_info_key_index("Feature_type")
        feature_index = self.variant.get_info_key_index("Feature")
        sift_index = self.variant.get_info_key_index("SIFT")
        polyphen_index = self.variant.get_info_key_index("PolyPhen")
        consequence_index = self.variant.get_info_key_index("Consequence")
        spdi_index = self.variant.get_info_key_index("SPDI")
        cadd_index = self.variant.get_info_key_index("CADD_PHRED")
        gerp_index = self.variant.get_info_key_index("Conservation")
        info_map = {}
        for csq_record in self.variant.info["CSQ"]:
            csq_record_list = csq_record.split("|")
            phenotype = None
            if re.search("--OMIM",csq_record_list[phenotypes_index]):
                phenotype = csq_record_list[phenotypes_index].split("--")[1]
            if csq_record_list[allele_index] in info_map.keys():
                if phenotype:
                    info_map[csq_record_list[allele_index]]["phenotype_assertions"].append(self.create_allele_phenotype_assertion(csq_record_list[feature_index], csq_record_list[feature_type_index], phenotype ))
                info_map[csq_record_list[allele_index]]["predicted_molecular_consequences"].append(self.create_allele_predicted_molecular_consequence(csq_record_list[spdi_index], 
                                                                                                                                               csq_record_list[feature_index], 
                                                                                                                                               csq_record_list[feature_type_index], 
                                                                                                                                               csq_record_list[consequence_index],
                                                                                                                                               csq_record_list[sift_index], 
                                                                                                                                               csq_record_list[polyphen_index]
                                                                                                                                   ))
                current_prediction_results = info_map[csq_record_list[allele_index]]["prediction_results"]
                info_map[csq_record_list[allele_index]]["prediction_results"] += self.create_allele_prediction_results(current_prediction_results, 
                                                                                                                                               csq_record_list[spdi_index], 
                                                                                                                                               csq_record_list[cadd_index],
                                                                                                                                               csq_record_list[gerp_index]
                                                                                                                                               )
                
            else:
                info_map[csq_record_list[allele_index]] = {"phenotype_assertions": [], "predicted_molecular_consequences": [], "prediction_results": []} 
                if phenotype:
                    info_map[csq_record_list[allele_index]]["phenotype_assertions"].append(self.create_allele_phenotype_assertion(csq_record_list[feature_index], csq_record_list[feature_type_index], phenotype))
                info_map[csq_record_list[allele_index]]["predicted_molecular_consequences"].append(self.create_allele_predicted_molecular_consequence(csq_record_list[spdi_index], 
                                                                                                                                               csq_record_list[feature_index], 
                                                                                                                                               csq_record_list[feature_type_index], 
                                                                                                                                               csq_record_list[consequence_index],
                                                                                                                                               csq_record_list[sift_index], 
                                                                                                                                               csq_record_list[polyphen_index]
                                                                                                                                               ))
                info_map[csq_record_list[allele_index]]["prediction_results"] += self.create_allele_prediction_results([], csq_record_list[spdi_index],
                                                                                                                                               csq_record_list[cadd_index],
                                                                                                                                               csq_record_list[gerp_index]
                                                                                                                                               )

        return info_map
    
    def create_allele_prediction_results(self, current_prediction_results: Mapping, allele: str, cadd_score: str, gerp_score: str) -> list:
        prediction_results = []
        if cadd_score:
            if not self.prediction_result_already_exists(current_prediction_results, "CADD"):
                cadd_prediction_result = {
                        "result": cadd_score ,
                        "analysis_method": {
                            "tool": "CADD",
                            "qualifier": "CADD"
                        }

                }
                prediction_results.append(cadd_prediction_result)
        if gerp_score:
            if not self.prediction_result_already_exists(current_prediction_results, "GERP"):
                gerp_prediction_result = {
                        "result": gerp_score ,
                        "analysis_method": {
                            "tool": "GERP",
                            "qualifier": "GERP"
                        }
                    }
                prediction_results.append(gerp_prediction_result)
        
        return prediction_results
    
    def prediction_result_already_exists(self, current_prediction_results: Mapping, tool: str) -> bool:
        for prediction_result in current_prediction_results:
            if prediction_result["analysis_method"]["tool"] == tool:
                return True

        return False
    
    def create_allele_predicted_molecular_consequence(self, allele: str, feature: str, feature_type: str, consequences: str, sift_score: str, polyphen_score: str) -> Mapping:
        consequences_list = []
        for cons in consequences.split("&"):
            consequences_list.append(
                {
                    "accession_id": cons
                }
            )

        prediction_results = []
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

        
        return {
            "allele_name": allele,
            "feature_stable_id": feature,
            "feature_type": {
                "accession_id": feature_type                
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
    
    def create_allele_phenotype_assertion(self, feature: str, feature_type: str , phenotype: str) -> Mapping:
        return {
            "feature": feature,
            "feature_type": {
                "accession_id": feature_type                
            } ,
            "phenotype": {
                "name": phenotype
            },
            "evidence": []

        }
    
    
    def minimise_allele(self, alt: str):
        """
        VCF file has the representation without anchoring bases
        for prediction scores in INFO column. This function is useful
        in matching the SPDI format in VCF with the allele in memory
        """
        minimised_allele = alt
        if len(alt) > len(self.variant.ref):
            minimised_allele = alt[1:] 
        elif len(alt) < len(self.variant.ref):
            minimised_allele = "-"
        return minimised_allele
    
    def format_frequency(self, raw_frequency_list: List) -> Mapping:
        freq_map = {}
        for freq in raw_frequency_list:
            key = freq.split(":")[0]
            freq_list = freq.split(":")[1].split(",")
            freq_map[key] = freq_list
        return freq_map
    
    def get_population_allele_frequencies(self) -> List:
        frequency_map = {}
        if "FREQ" in self.variant.info:
            frequency_map = self.format_frequency(",".join(map(str,self.variant.info["FREQ"])).split("|"))
        population_allele_frequencies = []
        for key, pop_list in frequency_map.items():
            ## Adding only GnomAD population
            if key == "GnomAD":
                if pop_list[self.allele_index] not in ["None", "."]:
                    population_allele_frequencies.append({
                        "population": key,
                        "allele_frequency": pop_list[self.allele_index],
                        "is_minor_allele": False,
                        "is_hpmaf": False
                    })
        return population_allele_frequencies

  
