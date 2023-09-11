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
    """
    Client to load file in-memory
    """
    def __init__(self, config, datafile=None):
        mapping_file = config.get("mapping_file")
        self.data_root = config.get("data_root")
        with open(mapping_file) as f:
            self.mapping = json.load(f)
    
    def get_variant_record(self, genome_uuid: str, variant_id: str):
        """
        Get a variant entry from variant_id
        """
        genome_metadata = next(d for d in self.mapping if d['genome_uuid'] == genome_uuid)
        datafile = glob.glob(os.path.join(self.data_root,genome_metadata["data_folder"],"*.vcf.gz"))[0]
        if datafile:
            self.collection = vcfpy.Reader.from_path(datafile)
            self.header = self.collection.header
        else:
            print("Please check the directory path for the given genome uuid")

        try: 
            [contig, pos, id] = self.split_variant_id(variant_id)
            pos = int(pos)
        except:
            #TODO: This needs to go to thoas logger
            #TODO: Exception needs to be caught appropriately
            print("Please check that the variant_id is in the format: contig:position:identifier")
        data = {}
        variant = None
        try:
            for rec in self.collection.fetch(contig, pos-1, pos):
                if rec.ID[0] == id:
                    variant = Variant(rec, self.header)
                    break
            return variant
        except:
            # Return None when variant cannot be fetched
            return
        
        
    def split_variant_id(self, variant_id: str):
        """
        Splits variant_id into separate fields
        """
        return variant_id.split(":")


