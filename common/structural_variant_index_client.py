"""
DuckDB-backed structural variant index client.
"""

from common.base_file_client import BaseFileClient


class StructuralVariantIndexClient(BaseFileClient):
    """Client for structural variant DuckDB indexes."""

    def __init__(self, sv_dir: str) -> None:
        self.sv_dir = sv_dir
        self.db_path = self._get_index_db_path()

    def exists(self) -> bool:
        return self.db_path is not None

    def get_variant_locations(self, structural_variant_id: str) -> list:
        rows = self.fetch_duckdb_rows(
            self.db_path,
            """
            SELECT chr, start, source_file
            FROM variants
            WHERE id = ?
            """,
            [structural_variant_id],
        )

        return [
            {
                "chr": str(chromosome),
                "start": int(start),
                "source_file": source_file,
            }
            for chromosome, start, source_file in rows
        ]

    def get_synonyms_by_allele_id(self, structural_variant_id: str) -> dict:
        rows = self.fetch_duckdb_rows(
            self.db_path,
            """
            SELECT DISTINCT allele_id, alias, source
            FROM variant_aliases
            WHERE id = ? AND alias IS NOT NULL
            ORDER BY alias
            """,
            [structural_variant_id],
        )

        synonyms_by_allele_id = {}
        for allele_id, synonym, source in rows:
            synonyms_by_allele_id.setdefault(allele_id, []).append(
                {"synonym": synonym, "source": source}
            )

        return synonyms_by_allele_id

    def load_structural_variant_synonyms(self, structural_variant) -> None:
        structural_variant.set_synonyms_by_allele_id(
            self.get_synonyms_by_allele_id(structural_variant.name)
        )

    def _get_index_db_path(self) -> str | None:
        return self.get_index_db_path(self.sv_dir)
