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
from pathlib import Path
from common.file_model.variant import Variant
from common.file_model.structural_variant import StructuralVariant
from common.structural_variant_searcher import StructuralVariantSearcher


class FileClient:
    """Client to load file into memory.

    This class provides methods to retrieve and process variant data from VCF files.
    """

    def __init__(self, config):
        """Initialises a FileClient instance.

        Args:
            config (dict): A configuration dictionary containing at least the "data_root" key.
        """
        self.data_root = config.get("data_root")
        self.sv_searcher = StructuralVariantSearcher()

    def get_variant_record(self, genome_uuid: str, variant_id: str):
        """Retrieves a variant entry using the specified variant identifier.

        This method loads the VCF file from the configured data root and genome UUID,
        then attempts to fetch the variant record that matches the given variant_id.
        The variant_id should be formatted as 'contig:position:identifier'.

        Args:
            genome_uuid (str): The UUID for the genome.
            variant_id (str): The variant identifier in the format 'contig:position:identifier'.

        Returns:
            Variant or None: A Variant instance if a matching record is found; otherwise, None.
        """
        datafile = os.path.join(self.data_root, genome_uuid, "variation.vcf.gz")
        if datafile:
            self.collection = vcfpy.Reader.from_path(datafile)
            self.header = self.collection.header
        else:
            print("Please check the directory path for the given genome uuid")

        try:
            [contig, pos, id] = self.split_variant_id(variant_id)
            pos = int(pos)
        except:
            # TODO: This needs to go to thoas logger
            # TODO: Exception needs to be caught appropriately
            print(
                "Please check that the variant_id is in the format: contig:position:identifier"
            )
        data = {}
        variant = None
        try:
            for rec in self.collection.fetch(contig, pos - 1, pos):
                if rec.ID[0] == id:
                    variant = Variant(rec, self.header, genome_uuid)
                    break
            return variant
        except:
            # Return None when variant cannot be fetched
            return
        
    def get_structural_variant_record(self, genome_uuid: str, variant_id: str, source_name: str = None) -> StructuralVariant|None:
        """Retrieves a structural variant entry using the specified variant identifier.

        This method loads the VCF file from the configured data root, genome UUID and source name,
        then attempts to fetch the structural variant record that matches the given structural_variant_id.
        The structural_variant_id should be formatted as 'contig:position:identifier'.
        
        If source_name is not provided, searches across all available structural variation source files at once.

        Args:
            genome_uuid (str): The UUID for the genome.
            variant_id (str): The variant identifier in the format 'contig:position:identifier'.
            source_name (str, optional): The name of the source for the structural variant.
                                        If None, searches all available sources.

        Returns:
            StructuralVariant or None: A StructuralVariant instance if a matching record is found; otherwise, None.
        """
        try:
            [contig, pos, id] = self.split_variant_id(variant_id)
            pos = int(pos)
        except Exception as e:
            print(f"Invalid variant_id format '{variant_id}': {str(e)}. Expected format: contig:position:identifier")
            return None
        
        # If source_name is provided, search that specific file
        if source_name:
            datafile = os.path.join(self.data_root, genome_uuid, f"structural-variation/variation_{source_name}.vcf.gz")
            variant = self.sv_searcher.search_in_file(datafile, contig, pos, id, genome_uuid, source_name)
            if variant:
                return variant
            print(f"Structural variant not found in source '{source_name}'")
            return None
        else:
            # Search across all available structural variation sources at once using bcftools
            sv_dir = os.path.join(self.data_root, genome_uuid, "structural-variation")
            if not os.path.isdir(sv_dir):
                print(f"Structural variation directory not found: {sv_dir}")
                return None
            
            # Get all valid VCF files
            vcf_files = sorted([
                f for f in os.listdir(sv_dir)
                if f.startswith("variation_") and f.endswith(".vcf.gz")
            ])
            
            if not vcf_files:
                print(f"No structural variation VCF files found in {sv_dir}")
                return None
            
            # Search all files at once using bcftools
            variant = self.sv_searcher.search_all_files(sv_dir, vcf_files, contig, pos, id, genome_uuid)
            print("Searching variant across all sources using bcftools...")
            if variant:
                return variant
        
        return None

    def split_variant_id(self, variant_id: str):
        """Splits the variant identifier into its constituent parts.

        The variant identifier is expected to be formatted as 'contig:position:identifier'.

        Args:
            variant_id (str): The variant identifier string.

        Returns:
            list: A list containing the contig, position, and identifier.
        """
        return variant_id.split(":")
