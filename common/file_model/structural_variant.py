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

from typing import Any, List

from vcfpy import SymbolicAllele

from common.file_model.base_variant import BaseVariant
from common.file_model.structural_variant_allele import StructuralVariantAllele
from common.file_model.utils import SVTYPE_TO_TERM


class StructuralVariant(BaseVariant):
    """StructuralVariant model inherits shared behaviour from BaseVariant."""

    @staticmethod
    def _normalize_spdi_name(name: str | None) -> str | None:
        if not name:
            return name

        parts = name.split(":")
        if len(parts) < 4:
            return name

        if parts[0] and parts[1] and parts[2] == "":
            parts[2] = "0"
            return ":".join(parts)

        if parts[-1] == "":
            parts[-1] = "0"
            return ":".join(parts)

        if len(parts) >= 5 and parts[-2] == "":
            parts[-2] = "0"
            return ":".join(parts)

        return name

    def __init__(
        self,
        record: Any,
        header: Any,
        genome_uuid: str,
        synonyms_by_allele_id: dict | None = None,
    ) -> None:
        """Initialise SV-specific attributes and delegate shared setup.

        Sequence allele lengths use minimised REF/ALT sequences. Symbolic alleles
        use VCF INFO `SVLEN` or `END`, falling back to zero.
        """
        super().__init__(record, header, genome_uuid)
        self.type = "StructuralVariant"
        self.synonyms_by_allele_id = synonyms_by_allele_id or {}
        self.length = self._compute_length()

    def _compute_length(self) -> int:
        if not self.alts:
            return 0

        has_symbolic_alt = any(isinstance(alt, SymbolicAllele) for alt in self.alts)
        if not has_symbolic_alt:
            lengths = []
            for alt in self.alts:
                alt_value = alt.value if hasattr(alt, "value") else str(alt)
                ref_value = self.ref or ""
                # Remove shared bases, including the VCF padding base.
                start = 0
                while (
                    start < min(len(ref_value), len(alt_value))
                    and ref_value[start] == alt_value[start]
                ):
                    start += 1
                ref_end, alt_end = len(ref_value), len(alt_value)
                while (
                    ref_end > start
                    and alt_end > start
                    and ref_value[ref_end - 1] == alt_value[alt_end - 1]
                ):
                    ref_end -= 1
                    alt_end -= 1
                ref_length = ref_end - start
                alt_length = alt_end - start
                # Insertions and indels use the effective sequence added;
                # deletions use the reference sequence removed. Equal-length
                # substitutions (including SNPs) have the same length either way.
                lengths.append(alt_length if alt_length else ref_length)
            return max(lengths)

        svlen = None
        if "SVLEN" in self.info:
            svlen = self.info.get("SVLEN") if isinstance(self.info, dict) else None
        elif "END" in self.info:
            svlen = self.info.get("END") - self.position + 1 if isinstance(self.info, dict) else None

        try:
            if isinstance(svlen, (list, tuple)) and svlen:
                return abs(max(map(int, svlen), key=abs))
            if svlen is not None:
                return abs(int(svlen))
        except Exception:
            pass

        return 0

    def set_synonyms_by_allele_id(self, synonyms_by_allele_id: dict | None) -> None:
        self.synonyms_by_allele_id = synonyms_by_allele_id or {}

    def get_alternative_names(self, allele_id: str | None = None) -> list:
        if allele_id:
            return [
                self._build_synonym(synonym["synonym"], synonym.get("source"))
                for synonym in self.synonyms_by_allele_id.get(allele_id, [])
            ]

        synonyms = []
        seen = set()
        for allele_synonyms in self.synonyms_by_allele_id.values():
            for synonym in allele_synonyms:
                synonym_id = synonym["synonym"]
                source = synonym.get("source")
                key = (
                    synonym_id,
                    source,
                )
                if key not in seen:
                    seen.add(key)
                    synonyms.append(self._build_synonym(synonym_id, source))

        return synonyms

    @staticmethod
    def _build_synonym(synonym: str, source: str) -> dict:
        source = source or "dbVar"
        return {
            "accession_id": synonym,
            "name": synonym,
            "description": f"Structural variant synonym from {source}",
            "assignment_method": {
                "type": "DIRECT",
                "description": (
                    "A reference made by an external resource of annotation to an "
                    "Ensembl feature that Ensembl imports without modification"
                ),
            },
            "url": None,
            "source": {
                "id": source,
                "name": source,
                "description": "",
                "url": None,
                "release": None,
            },
        }

    @staticmethod
    def _build_allele_type_payload(allele_type: str, so_term: str) -> dict:
        return {
            "accession_id": allele_type,
            "value": allele_type,
            "url": f"http://sequenceontology.org/browser/current_release/term/{so_term}",
            "source": {
                "id": "",
                "name": "Sequence Ontology",
                "url": "www.sequenceontology.org",
                "description": "The Sequence Ontology...",
            },
        }

    def get_slice(self) -> dict:
        return self._get_slice_for_type(self.get_allele_type()["value"])

    def _get_slice_for_type(self, allele_type: str) -> dict:
        start = self.position
        length = self.length
        if allele_type == "insertion":
            length = 0
            end = start
        else:
            end = start + length - 1
        return {
                    "location": {"start": start, "end": end, "length": length},
                    "region": {
                        "name": self.chromosome,
                        "code": "chromosome",
                        "topology": "linear",
                        "so_term": "SO:0001217",
                    },
                    "strand": {"code": "forward", "value": 1},
                }

    def get_allele_type(self) -> dict:
        """Classify the whole variant across all alternate alleles."""
        if any(isinstance(alt, SymbolicAllele) for alt in self.alts or []):
            alts = {
                alt.value if isinstance(alt, SymbolicAllele) else str(alt)
                for alt in self.alts
            }
            if {"DUP", "DEL"} <= alts:
                return self._build_allele_type_payload(*SVTYPE_TO_TERM["CNV"])
            for svtype in ("DUP", "DEL", "INV", "CNV", "INS", "BND"):
                if svtype in alts:
                    return self._build_allele_type_payload(*SVTYPE_TO_TERM[svtype])
            return self._build_allele_type_payload("structural_variant", "SO:0001537")

        if not self.alts:
            return self._build_allele_type_payload("structural_variant", "SO:0001537")
        return super().get_allele_type()

    def get_length(self) -> int:
        return self.length

    def get_alleles(self) -> List:
        variant_allele_list = []
        alts = self.alts or []
        for index, alt in enumerate(alts):
            alt_value = alt.value if hasattr(alt, "value") else str(alt)
            variant_allele_list.append(StructuralVariantAllele(index + 1, alt_value, self))
        if self.ref:
            variant_allele_list.append(StructuralVariantAllele(0, self.ref, self))
        return variant_allele_list

    def get_prediction_results(self) -> list:
        return []

    def get_ensembl_website_display_data(self) -> dict:
        return {}

    def get_web_display_data(self) -> dict:
        return {}
