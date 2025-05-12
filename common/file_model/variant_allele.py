"""
.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at
       http://www.apache.org/licenses/LICENSE-2.0
   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

from typing import Any, Mapping, List, Union, Tuple
import re
import os
import json
import operator
from functools import reduce
from common.file_model.utils import minimise_allele

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
        min_alt = minimise_allele(self.alt, self.reference_sequence)
        return self.info_map[min_alt]["phenotype_assertions"] if min_alt in self.info_map else []

    def get_predicted_molecular_consequences(self):
        ## compute info_map only for these methods, first check if non-empty and compute only then
        min_alt = minimise_allele(self.alt, self.reference_sequence)
        return self.info_map[min_alt]["predicted_molecular_consequences"] if min_alt in self.info_map else []
    
    def get_prediction_results(self):
        min_alt = minimise_allele(self.alt, self.reference_sequence)
        return self.info_map[min_alt]["prediction_results"] if min_alt in self.info_map else []
    
    def get_population_allele_frequencies(self):
        min_alt = minimise_allele(self.alt, self.reference_sequence)
        population_map = self.variant.set_frequency_flags()
        return population_map[min_alt].values() if min_alt in population_map else []

    def get_web_display_data(self) -> Mapping:
        return self.variant.get_statistics_info()[self.allele_sequence]
    
    def traverse_csq_info(self) -> Mapping:
        """
        This function is to traverse the CSQ record and extract columns
        corresponding to Consequence, SIFT, PolyPhen, CADD
        """
        column_list = ["Allele", "PHENOTYPES", "Feature_type", "Feature", "Consequence",
                        "SIFT", "PolyPhen", "SpliceAI_pred_DP_AG", "SpliceAI_pred_DP_AL", 
                        "SpliceAI_pred_DP_DG", "SpliceAI_pred_DP_DL", "SpliceAI_pred_DS_AG", 
                        "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL", 
                        "SpliceAI_pred_SYMBOL", "SPDI", "CADD_PHRED","Conservation", "Gene", 
                        "SYMBOL","BIOTYPE","cDNA_position", "CDS_position", "Protein_position", 
                        "Amino_acids", "Codons"]
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

            # parse and form phenotypes - adding phenotype from any of the csq record would be enough for adding only variant-linked phenotypes
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
                            "qualifier": 
                             {
                                "result_type": "CADD Phred score"
                             },
                             "reference_data": []
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
    
    def parse_position(self, position: str)->Tuple:
        position_list  = position.split("-")
        position_start = position_list[0]
        position_end = position_start if len(position_list) < 2 else position_list[1]
        allele_type = self.get_allele_type()["accession_id"]
        if (position_start == "?" or position_end == "?"):
            position_start = position_start if position_start !="?" else None
            position_end = position_end if position_end !="?" else None
            position_length = None
        elif (allele_type == "insertion"):
            position_length = 0  # consistent logic for insertion
        else:
            position_length = int(position_end) - int(position_start) + 1
        return (position_start,position_end,position_length)
    
    def create_allele_predicted_molecular_consequence(self, csq_record: List, prediction_index_map: dict) -> Mapping:
        feature_type = csq_record[prediction_index_map["feature_type"]]
        consequences_list = []
        if "consequence" in prediction_index_map.keys():
            for cons in csq_record[prediction_index_map["consequence"]].split("&"):
                if cons in ["downstream_gene_variant", "upstream_gene_variant", "intergenic_variant", "regulatory_region_variant", "TF_binding_site_variant"]:
                    consequences_list = []
                    break
                consequences_list.append(
                    {
                        "value": cons
                    }
                )

        prediction_results = []
        if "sift" in prediction_index_map.keys():
            sift_score = csq_record[prediction_index_map["sift"]]
            if sift_score:
                (label, score) = self.format_sift_polyphen_output(sift_score)
                if label is not None and score is not None:
                    classification = {
                        "label": label,
                        "definition": ""
                    }
                    sift_prediction_result = {
                            "classification": classification,
                            "score": score,
                            "analysis_method": {
                                "tool": "SIFT",
                                "version": "6.2.1", ## hardcoding until next data generation
                                "qualifier": {
                                    "result_type": "SIFT score",
                                    "modes": []
                                },
                                "reference_data": []
                            }
                        }
                    prediction_results.append(sift_prediction_result)
        
        if "spliceai_pred_dp_dg" in prediction_index_map:
            spliceai = {}
            spliceai["spliceai_pred_dp_dg"] = {
                "score": csq_record[prediction_index_map["spliceai_pred_dp_dg"]],
                "qualifier": self.variant.header.get_lines("SpliceAI_pred_DP_DG")[0].value.split(". ")[1]
            }
            spliceai["spliceai_pred_dp_dl"] = {
                "score": csq_record[prediction_index_map["spliceai_pred_dp_dl"]],
                "qualifier": self.variant.header.get_lines("SpliceAI_pred_DP_DL")[0].value.split(". ")[1]
            }
            spliceai["spliceai_pred_ds_dg"] = {
            "score": csq_record[prediction_index_map["spliceai_pred_ds_dg"]],
            "qualifier": self.variant.header.get_lines("SpliceAI_pred_DS_DG")[0].value.split(". ")[1]
            }
            spliceai["spliceai_pred_ds_dl"] = {
            "score": csq_record[prediction_index_map["spliceai_pred_ds_dl"]],
            "qualifier": self.variant.header.get_lines("SpliceAI_pred_DS_DL")[0].value.split(". ")[1]
            }
            spliceai["spliceai_pred_dp_ag"] = {
            "score": csq_record[prediction_index_map["spliceai_pred_dp_ag"]],
            "qualifier": self.variant.header.get_lines("SpliceAI_pred_DP_AG")[0].value.split(". ")[1]
            }
            spliceai["spliceai_pred_dp_al"] = {
            "score": csq_record[prediction_index_map["spliceai_pred_dp_al"]],
            "qualifier": self.variant.header.get_lines("SpliceAI_pred_DP_AL")[0].value.split(". ")[1]
            }
            spliceai["spliceai_pred_ds_ag"] = {
            "score": csq_record[prediction_index_map["spliceai_pred_ds_ag"]],
            "qualifier": self.variant.header.get_lines("SpliceAI_pred_DS_AG")[0].value.split(". ")[1]
            }
            spliceai["spliceai_pred_ds_al"] = {
            "score": csq_record[prediction_index_map["spliceai_pred_ds_al"]],
            "qualifier": self.variant.header.get_lines("SpliceAI_pred_DS_AL")[0].value.split(". ")[1]
            }

            for spliceai_score in spliceai.values():
                if spliceai_score["score"]:
                    spliceai_prediction_result  = {
                            "score": spliceai_score["score"],
                            "analysis_method": {
                                "tool": "SpliceAI",
                                "version": "1.3.1",  ## hardcoding until next data generation
                                "qualifier": {
                                    "result_type": spliceai_score["qualifier"],
                                    "modes": ["masked scores"] ## hardcoding until next data generation

                                },
                                "reference_data": []
                            }
                        }
                    prediction_results.append(spliceai_prediction_result)

        if "polyphen" in prediction_index_map.keys():
            polyphen_score = csq_record[prediction_index_map["polyphen"]]
            if polyphen_score:
                (label, score) = self.format_sift_polyphen_output(polyphen_score)
                if label is not None and score is not None:
                    classification = {
                        "label": label,
                        "definition": ""
                    }
                    polyphen_prediction_result = {
                            "classification": classification,
                            "score": score,
                            "analysis_method": {
                                "tool": "PolyPhen",
                                "version": "2.2.3", ## hardcoding until next data generation
                                "qualifier": {
                                    "result_type": "PolyPhen score",
                                    "modes": ["HumVar"] ## hardcoding until next data generation
                                },
                                "reference_data": []
                            }
                        }
                    prediction_results.append(polyphen_prediction_result)
        
        cdna_location = cds_location = protein_location = None

        ###parse cdna location
        cdna_position = csq_record[prediction_index_map["cdna_position"]]
        codons = csq_record[prediction_index_map["codons"]]
        ref_cdna_sequence = alt_cdna_sequence = None
        if cdna_position:
            cdna_start, cdna_end, cdna_length = self.parse_position(cdna_position)
            if (cdna_start != None and cdna_end != None):
                if codons:
                    ref_cdna_sequence = re.sub('([a-z])','',codons.split("/")[0]) 
                    alt_cdna_sequence = re.sub('([a-z])','',codons.split("/")[1])
                    
                else:
                    # TODO: Handle when strand is reverse
                    ref_cdna_sequence = minimise_allele(self.reference_sequence, self.reference_sequence)
                    alt_cdna_sequence = csq_record[prediction_index_map["allele"]]
                ref_cdna_sequence = ref_cdna_sequence if ref_cdna_sequence else "-"
                alt_cdna_sequence = alt_cdna_sequence if alt_cdna_sequence else "-"
            cdna_location = {
                "start": cdna_start,
                "end": cdna_end, 
                "length": cdna_length, 
                "ref_sequence": ref_cdna_sequence,
                "alt_sequence": alt_cdna_sequence
            }

        ###parse cds location
        cds_position = csq_record[prediction_index_map["cds_position"]]
        codons = csq_record[prediction_index_map["codons"]]
        ref_cds_sequence = alt_cds_sequence = None
        if cds_position:
            cds_start, cds_end, cds_length = self.parse_position(cds_position)
            if (cds_start != None and cds_end != None):
                ref_cds_sequence = codons.split("/")[0]
                alt_cds_sequence = codons.split("/")[1]
            cds_location = {
                "start": cds_start,
                "end": cds_end, 
                "length": cds_length, 
                "ref_sequence": ref_cds_sequence,
                "alt_sequence": alt_cds_sequence
            }

        ###parse protein location
        protein_position = csq_record[prediction_index_map["protein_position"]]
        amino_acids = csq_record[prediction_index_map["amino_acids"]]
        ref_protein_sequence = alt_protein_sequence = None
        if protein_position:
            protein_start, protein_end, protein_length = self.parse_position(protein_position)
            if (protein_start != None and protein_end != None):
                amino_acids_array = amino_acids.split("/")
                ref_protein_sequence = amino_acids_array[0]
                alt_protein_sequence = amino_acids_array[1] if len(amino_acids_array)>1 else amino_acids_array[0]
            protein_location = {
                "start": protein_start,
                "end": protein_end, 
                "length": protein_length, 
                "ref_sequence": ref_protein_sequence,
                "alt_sequence": alt_protein_sequence
            }

        if consequences_list and feature_type == "Transcript":
            return {
                "allele_name": csq_record[prediction_index_map["allele"]],
                "stable_id": csq_record[prediction_index_map["feature"]],
                "feature_type": {  
                    "value":    feature_type  
                } ,
                "consequences": consequences_list,
                "gene_stable_id" : csq_record[prediction_index_map["gene"]], 
                "gene_symbol": csq_record[prediction_index_map["symbol"]], 
                "transcript_biotype": csq_record[prediction_index_map["biotype"]], 
                "prediction_results": prediction_results,
                "cdna_location": cdna_location,
                "cds_location": cds_location,
                "protein_location": protein_location
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

        feature_type=None
        splits = phenotype.split("+")
        if splits[2].startswith("ENS"):
            return None
        if len(splits)==3:
            phenotype_name,source_id,feature_id = splits
        elif len(splits)==5:
            phenotype_name,source_id,feature_id,feature_type,clinvar_clin_sig = splits
        else:
            return None

        if not feature_type:
            if re.search("^rs", feature_id):
                feature_type = "Variation"  
            else:
                feature_type = None

        evidence_list = []
            
        if phenotype:
            return {
                "feature": feature_id,

                "feature_type": {
                    "value": feature_type   

                } ,
                "phenotype": {
                    "name": phenotype_name,
                    "source": {
                        "id": source_id,
                        "name": source_id.replace("_"," ")
                    }
                },
                "evidence": evidence_list
            }

    
