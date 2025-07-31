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
import json
import os
import glob
from common.file_model.variant import Variant

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
            print("Please check that the variant_id is in the format: contig:position:identifier")
        data = {}
        variant = None
        try:
            for rec in self.collection.fetch(contig, pos-1, pos):
                if rec.ID[0] == id:
                    variant = Variant(rec, self.header, genome_uuid)
                    break
            return variant
        except:
            # Return None when variant cannot be fetched
            return
        
    def split_variant_id(self, variant_id: str):
        """Splits the variant identifier into its constituent parts.

        The variant identifier is expected to be formatted as 'contig:position:identifier'.

        Args:
            variant_id (str): The variant identifier string.

        Returns:
            list: A list containing the contig, position, and identifier.
        """
        return variant_id.split(":")


