"""
Structural variant allele model.
"""

from typing import Mapping, List


class StructuralVariantAllele:
    def __init__(self, allele_index: int, alt: str, variant: dict) -> None:
        ref_len = len(variant.ref) if variant.ref is not None else 0

        # Prefer SVLEN entry from INFO when available, else fallback to actual alt sequence length.
        alt_len = 0
        if hasattr(variant, "info") and isinstance(variant.info, dict):
            raw_svlen = variant.info.get("SVLEN")
            svlen_value = None
            if isinstance(raw_svlen, (list, tuple)):
                if allele_index > 0 and allele_index - 1 < len(raw_svlen):
                    svlen_value = raw_svlen[allele_index - 1]
            else:
                svlen_value = raw_svlen

            if svlen_value is not None:
                try:
                    alt_len = abs(int(svlen_value))
                except (TypeError, ValueError):
                    pass

        allele_name = None
        if hasattr(variant, "info") and isinstance(variant.info, dict):
            raw_allele_name = variant.info.get("ALLELE_NAME")
            if isinstance(raw_allele_name, (list, tuple)):
                # allele_index of 1 means first ALT; 0 refers to ref
                if allele_index > 0 and allele_index - 1 < len(raw_allele_name):
                    allele_name = raw_allele_name[allele_index - 1]
                elif allele_index == 0 and variant.ref is not None:
                    # no ref-specific name in list, keep fallback
                    allele_name = None
            elif isinstance(raw_allele_name, str):
                allele_name = raw_allele_name

        if allele_name:
            self.name = allele_name
        else:
            self.name = f"{variant.chromosome}:{variant.position}:{ref_len}:{alt_len}"

        self.variant = variant
        self.allele_index = allele_index
        self.alt = alt
        self.type = "StructuralVariantAllele"

    def get_allele_type(self) -> Mapping:
        return self.variant.get_allele_type()

    def get_alternative_names(self) -> list:
        return self.variant.get_alternative_names()

    def get_slice(self) -> Mapping:
        return self.variant.get_slice(self.alt)

    def get_phenotype_assertions(self) -> list:
        return []

    def get_predicted_molecular_consequences(self) -> list:
        return []

    def get_prediction_results(self) -> list:
        return []

    def get_population_allele_frequencies(self) -> list:
        return []

    def get_web_display_data(self) -> Mapping:
        return {}

    def create_allele_prediction_results(
        self,
        current_prediction_results: Mapping,
        csq_record: List,
        prediction_index_map: dict,
    ) -> list:
        """Creates prediction results for the allele based on a CSQ record.

        Args:
            current_prediction_results (Mapping): Existing prediction results.
            csq_record (List): The CSQ record split into fields.
            prediction_index_map (dict): A mapping from annotation keys to their indices.

        Returns:
            list: A list of new prediction results.
        """
        prediction_results = []
        if "cadd_phred" in prediction_index_map.keys():
            if not self.prediction_result_already_exists(
                current_prediction_results, "CADD"
            ):
                cadd_prediction_result = (
                    {
                        "score": csq_record[prediction_index_map["cadd_phred"]],
                        "analysis_method": {
                            "tool": "CADD",
                            "qualifier": {"result_type": "CADD Phred score"},
                            "reference_data": [],
                        },
                    }
                    if csq_record[prediction_index_map["cadd_phred"]]
                    else None
                )
                if cadd_prediction_result:
                    prediction_results.append(cadd_prediction_result)

        return prediction_results

    def prediction_result_already_exists(
        self, current_prediction_results: Mapping, tool: str
    ) -> bool:
        """Checks if a prediction result for a specific tool already exists.

        Args:
            current_prediction_results (Mapping): The current prediction results.
            tool (str): The analysis tool name.

        Returns:
            bool: True if a result exists, False otherwise.
        """
        for prediction_result in current_prediction_results:
            if prediction_result["analysis_method"]["tool"] == tool:
                return True

        return False
