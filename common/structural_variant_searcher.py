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
from common.file_model.structural_variant import StructuralVariant


class StructuralVariantSearcher:
    """Handles searching for structural variants across VCF files using bcftools."""

    def search_in_file(
        self, datafile: str, contig: str, pos: int, id: str, genome_uuid: str, source_name: str = None
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
            return self._search_with_bcftools(datafile, contig, pos, id, genome_uuid, source_name)
        except Exception as e:
            print(f"bcftools search failed for {source_name or datafile}: {str(e)}")
            return None

    def search_all_files(
        self, sv_dir: str, vcf_files: list, contig: str, pos: int, id: str, genome_uuid: str
    ) -> StructuralVariant | None:
        """Search for a structural variant across all VCF files at once using bcftools.

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
            # Build file paths for all VCF files
            file_paths = [os.path.join(sv_dir, f) for f in vcf_files]
            print(id, file_paths)

            # Use bcftools to search all files at once with ID filter
            result = subprocess.run(
                ["bcftools", "query", "-i", f'ID="{id}"', "-f", "%CHROM\t%POS\t%ID\t%FILENAME\n"]
                + file_paths,
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )

            if result.returncode != 0 or not result.stdout.strip():
                return None

            # Parse the output to find matching variant
            for line in result.stdout.strip().split("\n"):
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) >= 4:
                    query_chrom, query_pos, query_id, query_file = parts[0], int(parts[1]), parts[2], parts[3]
                    if query_chrom == contig and query_pos == pos and query_id == id:
                        # Found it, now load the full record with bcftools view
                        variant = self._get_full_record(query_file, contig, pos, id, genome_uuid)
                        if variant:
                            return variant

            return None
        except subprocess.TimeoutExpired:
            print("bcftools timeout searching all structural variation files")
            return None
        except Exception as e:
            print(f"Error searching all structural variation files: {str(e)}")
            return None

    def _search_with_bcftools(
        self, datafile: str, contig: str, pos: int, id: str, genome_uuid: str, source_name: str = None
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
            reader = vcfpy.Reader.from_string(result.stdout)
            for rec in reader:
                if rec.CHROM == contig and rec.POS == pos and rec.ID and rec.ID[0] == id:
                    return StructuralVariant(rec, reader.header, genome_uuid)

            return None
        except subprocess.TimeoutExpired:
            print(f"bcftools timeout for {source_name or datafile}")
            return None
        except Exception as e:
            raise

    def _get_full_record(
        self, datafile: str, contig: str, pos: int, id: str, genome_uuid: str
    ) -> StructuralVariant | None:
        """Get the full VCF record for a variant using bcftools view.

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

            reader = vcfpy.Reader.from_string(result.stdout)
            for rec in reader:
                if rec.CHROM == contig and rec.POS == pos and rec.ID and rec.ID[0] == id:
                    return StructuralVariant(rec, reader.header, genome_uuid)

            return None
        except Exception as e:
            print(f"Error fetching full record from {datafile}: {str(e)}")
            return None
