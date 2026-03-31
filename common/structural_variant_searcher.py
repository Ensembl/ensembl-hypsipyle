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

import vcfpy
import os
import subprocess
import shutil
import io
from common.file_model.structural_variant import StructuralVariant


class StructuralVariantSearcher:
    """Handles searching for structural variants across indexed VCF files."""

    def search_in_file(
        self,
        datafile: str,
        contig: str,
        pos: int,
        id: str,
        genome_uuid: str,
        source_name: str = None,
    ) -> StructuralVariant | None:
        """Search for a structural variant in a specific VCF file.

        Args:
            datafile (str): Path to the VCF file.
            contig (str): The contig/chromosome name.
            pos (int): The position.
            id (str): The variant identifier.
            genome_uuid (str): The genome UUID.
            source_name (str, optional): The source name for logging purposes.

        Returns:
            StructuralVariant or None: A StructuralVariant instance if found; otherwise, None.
        """
        if not os.path.exists(datafile):
            return None

        try:
            return self._search_with_bcftools(
                datafile, contig, pos, id, genome_uuid, source_name
            )
        except Exception as e:
            print(f"bcftools search failed for {source_name or datafile}: {str(e)}")
            return None

    def search_all_files(
        self,
        sv_dir: str,
        vcf_files: list,
        contig: str,
        pos: int,
        id: str,
        genome_uuid: str,
    ) -> StructuralVariant | None:
        """Search for a structural variant across all VCF files.

        This helper uses a DuckDB index (variants.duckdb) for speed.


        Args:
            sv_dir (str): Path to the structural-variation directory.
            vcf_files (list): List of VCF filenames to search.
            contig (str): The contig/chromosome name.
            pos (int): The position.
            id (str): The variant identifier.
            genome_uuid (str): The genome UUID.

        Returns:
            StructuralVariant or None: A StructuralVariant instance if found; otherwise, None.
        """
        try:
            db_path = os.path.join(sv_dir, "variants.duckdb")

            # if only a single file is supplied, we can short-circuit into the
            # single-file searcher rather than using the generic multi-file scan.
            # however we still honour the DuckDB index if one exists.
            if len(vcf_files) == 1:
                single_file = os.path.join(sv_dir, vcf_files[0])
                # attempt duckdb-based lookup first if index exists
                if os.path.exists(db_path):
                    variant = self._search_with_duckdb(
                        db_path, sv_dir, contig, pos, id, genome_uuid
                    )
                    if variant:
                        return variant
                    # index was present but gave no match; fall back to single-file search
                    print(
                        "DuckDB index search returned no matching record, falling back to single-file bcftools search"
                    )
                else:
                    print(
                        f"DuckDB index not found ({db_path}), using single-file search"
                    )
                # perform the one-file search and return whatever we get
                return self.search_in_file(single_file, contig, pos, id, genome_uuid)

            # attempt duckdb-based lookup first if index exists
            if os.path.exists(db_path):
                variant = self._search_with_duckdb(
                    db_path, sv_dir, contig, pos, id, genome_uuid
                )
                if variant:
                    return variant
                print("DuckDB index search returned no matching record")
            else:
                print(f"DuckDB index not found ({db_path})")
            return None
        except Exception as e:
            print(f"Error searching all structural variation files: {str(e)}")
            return None

    def _search_with_duckdb(
        self,
        db_path: str,
        sv_dir: str,
        contig: str,
        pos: int,
        id: str,
        genome_uuid: str,
    ) -> StructuralVariant | None:
        """Search for a structural variant using DuckDB index.

        Args:
            db_path (str): Path to the variants.duckdb database.
            sv_dir (str): Path to the structural-variation directory.
            contig (str): The contig/chromosome name.
            pos (int): The position.
            id (str): The variant identifier.
            genome_uuid (str): The genome UUID.

        Returns:
            StructuralVariant or None: A StructuralVariant instance if found; otherwise, None.
        """
        try:
            # Import lazily so environments without DuckDB can still use bcftools fallback.
            import duckdb

            conn = duckdb.connect(db_path, read_only=True)

            # Query the index by variant ID
            results = conn.execute(
                "SELECT chr, start, source_file FROM variants WHERE id = ?",
                [id],
            )
            results = results.fetchall()
            conn.close()

            if not results:
                print(f"Variant {id} not found in DuckDB index")
                return None

            # print(f"DuckDB lookup for {id}: found {len(results)} result(s)")

            # Check if any result matches the contig and position
            for result_chrom, result_pos, filename in results:
                # print(
                #     f"  Checking: chrom={result_chrom} (type={type(result_chrom).__name__}), pos={result_pos} (type={type(result_pos).__name__}), file={filename}"
                # )
                # print(
                #     f"  Against:  contig={contig} (type={type(contig).__name__}), pos={pos} (type={type(pos).__name__})"
                # )

                # Convert to same types for comparison
                result_chrom = str(result_chrom)
                result_pos = int(result_pos)

                if result_chrom == contig and result_pos == pos:
                    # print("  ✓ Coordinates match!")
                    # Found matching variant, load the full record
                    datafile = os.path.join(sv_dir, filename)
                    variant = self._get_full_record(
                        datafile, contig, pos, id, genome_uuid
                    )
                    if variant:
                        return variant
                else:
                    print("  ✗ Coordinates don't match")

            print(f"Variant {id} found in index but no coordinates matched")
            return None

        except Exception as e:
            print(f"Error searching DuckDB index: {str(e)}")
            import traceback

            traceback.print_exc()
            return None

    def _search_with_bcftools(
        self,
        datafile: str,
        contig: str,
        pos: int,
        id: str,
        genome_uuid: str,
        source_name: str = None,
    ) -> StructuralVariant | None:
        """Search for a structural variant using bcftools view.

        Uses bcftools view to extract the exact region and then parses the output.

        Args:
            datafile (str): Path to the VCF file.
            contig (str): The contig/chromosome name.
            pos (int): The position.
            id (str): The variant identifier.
            genome_uuid (str): The genome UUID.
            source_name (str, optional): The source name for logging purposes.

        Returns:
            StructuralVariant or None: A StructuralVariant instance if found; otherwise, None.
        """
        try:
            # Extract the exact region using bcftools view
            region = f"{contig}:{pos}-{pos}"
            result = subprocess.run(
                ["bcftools", "view", "-r", region, datafile],
                capture_output=True,
                text=True,
                timeout=10,
                check=False,
            )

            if result.returncode != 0:
                return None

            # Parse VCF output to find matching variant
            reader = vcfpy.Reader(io.StringIO(result.stdout))
            for rec in reader:
                if (
                    rec.CHROM == contig
                    and rec.POS == pos
                    and rec.ID
                    and rec.ID[0] == id
                ):
                    return StructuralVariant(rec, reader.header, genome_uuid)

            return None
        except subprocess.TimeoutExpired:
            print(f"bcftools timeout for {source_name or datafile}")
            return None
        except Exception:
            raise

    def _get_full_record(
        self, datafile: str, contig: str, pos: int, id: str, genome_uuid: str
    ) -> StructuralVariant | None:
        """Get the full VCF record for a variant using bcftools view or vcfpy fallback.

        Args:
            datafile (str): Path to the VCF file.
            contig (str): The contig/chromosome name.
            pos (int): The position.
            id (str): The variant identifier.
            genome_uuid (str): The genome UUID.

        Returns:
            StructuralVariant or None: A StructuralVariant instance if found; otherwise, None.
        """
        try:
            # Try bcftools first
            if shutil.which("bcftools"):
                return self._get_full_record_bcftools(
                    datafile, contig, pos, id, genome_uuid
                )
            else:
                # Fallback to vcfpy if bcftools is not available
                print(f"bcftools not found, falling back to vcfpy for {datafile}")
                return self._get_full_record_vcfpy(
                    datafile, contig, pos, id, genome_uuid
                )
        except Exception as e:
            print(f"Error fetching full record from {datafile}: {str(e)}")
            return None

    def _get_full_record_bcftools(
        self, datafile: str, contig: str, pos: int, id: str, genome_uuid: str
    ) -> StructuralVariant | None:
        """Get the full VCF record using bcftools view.

        Args:
            datafile (str): Path to the VCF file.
            contig (str): The contig/chromosome name.
            pos (int): The position.
            id (str): The variant identifier.
            genome_uuid (str): The genome UUID.

        Returns:
            StructuralVariant or None: A StructuralVariant instance if found; otherwise, None.
        """
        try:
            region = f"{contig}:{pos}-{pos}"
            result = subprocess.run(
                ["bcftools", "view", "-r", region, datafile],
                capture_output=True,
                text=True,
                timeout=10,
                check=False,
            )

            if result.returncode != 0:
                return None

            reader = vcfpy.Reader(io.StringIO(result.stdout))
            for rec in reader:
                if (
                    rec.CHROM == contig
                    and rec.POS == pos
                    and rec.ID
                    and rec.ID[0] == id
                ):
                    return StructuralVariant(rec, reader.header, genome_uuid)

            return None
        except Exception as e:
            print(f"Error in bcftools record retrieval: {str(e)}")
            return None

    def _get_full_record_vcfpy(
        self, datafile: str, contig: str, pos: int, id: str, genome_uuid: str
    ) -> StructuralVariant | None:
        """Get the full VCF record using vcfpy.

        Args:
            datafile (str): Path to the VCF file.
            contig (str): The contig/chromosome name.
            pos (int): The position.
            id (str): The variant identifier.
            genome_uuid (str): The genome UUID.

        Returns:
            StructuralVariant or None: A StructuralVariant instance if found; otherwise, None.
        """
        try:
            reader = vcfpy.Reader.from_path(datafile)
            # Search for the variant using fetch if supported, otherwise iterate all records
            try:
                for rec in reader.fetch(contig, pos - 1, pos):
                    if rec.ID and rec.ID[0] == id:
                        return StructuralVariant(rec, reader.header, genome_uuid)
            except (NotImplementedError, TypeError):
                # If fetch is not supported or fails, iterate through all records
                reader.close()
                reader = vcfpy.Reader.from_path(datafile)
                for rec in reader:
                    if (
                        rec.CHROM == contig
                        and rec.POS == pos
                        and rec.ID
                        and rec.ID[0] == id
                    ):
                        return StructuralVariant(rec, reader.header, genome_uuid)

            reader.close()
            return None
        except Exception as e:
            print(f"Error in vcfpy record retrieval from {datafile}: {str(e)}")
            return None
