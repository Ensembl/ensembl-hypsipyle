"""
Shared file and index helpers for variant clients.
"""

import io
import os
import shutil
import subprocess
from common.file_model.variant import Variant
from common.file_model.structural_variant import StructuralVariant

import vcfpy


class BaseFileClient:
    """Common helpers for loading VCF records and locating related files."""

    INDEX_FILENAMES = ("variants.duckdb", "variants.db")

    def split_variant_id(self, variant_id: str):
        """Splits a variant identifier in the form contig:position:identifier."""
        return variant_id.split(":")

    def get_record_from_file(
        self,
        datafile: str,
        contig: str,
        pos: int,
        id: str,
        genome_uuid: str,
        record_class,
        fallback_to_iteration: bool = False,
        error_prefix: str = None,
    ) -> None | Variant| StructuralVariant:
        """Read a matching record from a VCF file using vcfpy."""
        if not os.path.exists(datafile):
            return None

        reader = None
        try:
            reader = vcfpy.Reader.from_path(datafile)
            try:
                for rec in reader.fetch(contig, pos - 1, pos):
                    if self.record_matches(rec, contig, pos, id):
                        return record_class(rec, reader.header, genome_uuid)
            except (NotImplementedError, TypeError):
                if not fallback_to_iteration:
                    return None
                self.close_reader(reader)
                reader = vcfpy.Reader.from_path(datafile)
                for rec in reader:
                    if self.record_matches(rec, contig, pos, id):
                        return record_class(rec, reader.header, genome_uuid)
        except Exception as e:
            if error_prefix:
                print(f"{error_prefix}: {str(e)}")
            return None
        finally:
            self.close_reader(reader)

        return None

    def get_record_with_bcftools(
        self,
        datafile: str,
        contig: str,
        pos: int,
        id: str,
        genome_uuid: str,
        record_class,
        source_name: str = None,
    ) -> None | Variant | StructuralVariant :
        """Read a matching record from a VCF file using bcftools view.
        bcftools gives us raw VCF text, but the rest of this codebase 
        expects parsed VCF record objects
        """
        region = f"{contig}:{pos}-{pos}"
        try:
            result = subprocess.run(
                ["bcftools", "view", "-r", region, datafile],
                capture_output=True,
                text=True,
                timeout=10,
                check=False,
            )
        except subprocess.TimeoutExpired:
            print(f"bcftools timeout for {source_name or datafile}")
            return None

        if result.returncode != 0:
            return None

        reader = vcfpy.Reader(io.StringIO(result.stdout))
        for rec in reader:
            if self.record_matches(rec, contig, pos, id):
                return record_class(rec, reader.header, genome_uuid)

        return None

    def get_full_record(
        self,
        datafile: str,
        contig: str,
        pos: int,
        id: str,
        genome_uuid: str,
        record_class,
    ) -> None | Variant | StructuralVariant:
        """Read a full VCF record using bcftools when available, otherwise vcfpy."""
        try:
            return self.get_record_from_file(
                datafile,
                contig,
                pos,
                id,
                genome_uuid,
                record_class,
                fallback_to_iteration=True,
            )
        except Exception as e:
            print(f"Error fetching full record from {datafile}: {str(e)}")
            return None

    def get_variant_source_files(
        self, genome_dir: str, source_name: str = None
    ) -> list[str]:
        """Find small-variant source VCFs under supported genome layouts."""
        bases = [genome_dir, os.path.join(genome_dir, "variation")]
        candidates = []

        if source_name:
            source_names = [source_name]
            lower_source = source_name.lower()
            if lower_source != source_name:
                source_names.append(lower_source)

            for base in bases:
                for name in source_names:
                    candidates.append(os.path.join(base, name, "variation.vcf.gz"))
                    candidates.append(os.path.join(base, f"variation_{name}.vcf.gz"))

                if os.path.isdir(base):
                    for entry in os.listdir(base):
                        if entry.lower() == lower_source:
                            candidates.append(
                                os.path.join(base, entry, "variation.vcf.gz")
                            )
                        if entry.lower() == f"variation_{lower_source}.vcf.gz":
                            candidates.append(os.path.join(base, entry))
        else:
            for base in bases:
                if not os.path.isdir(base):
                    continue

                for entry in sorted(os.listdir(base)):
                    path = os.path.join(base, entry)
                    if entry.startswith("variation_") and entry.endswith(".vcf.gz"):
                        candidates.append(path)
                    elif os.path.isdir(path) and entry != "structural-variation":
                        candidates.append(os.path.join(path, "variation.vcf.gz"))

        return self.existing_paths(candidates)

    def get_vcf_files(
        self, directory: str, prefix: str = "variation_", suffix: str = ".vcf.gz"
    ) -> list[str]:
        if not os.path.isdir(directory):
            return []

        return sorted(
            filename
            for filename in os.listdir(directory)
            if filename.startswith(prefix) and filename.endswith(suffix)
        )

    def get_index_db_path(self, directory: str) -> str | None:
        for filename in self.INDEX_FILENAMES:
            db_path = os.path.join(directory, filename)
            if os.path.exists(db_path):
                return db_path

        return None

    def fetch_duckdb_rows(self, db_path: str, query: str, parameters: list) -> list:
        if not db_path:
            return []

        try:
            import duckdb

            conn = duckdb.connect(db_path, read_only=True)
            try:
                return conn.execute(query, parameters).fetchall()
            finally:
                conn.close()
        except Exception:
            return []

    def existing_paths(self, candidates: list[str]) -> list[str]:
        existing = []
        seen = set()
        for path in candidates:
            if path not in seen and os.path.exists(path):
                existing.append(path)
                seen.add(path)
        return existing

    def record_matches(self, rec, contig: str, pos: int, id: str) -> bool:
        return rec.CHROM == contig and rec.POS == pos and rec.ID and rec.ID[0] == id

    def close_reader(self, reader) -> None:
        if reader is None:
            return

        close = getattr(reader, "close", None)
        if close:
            close()
