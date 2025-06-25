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
from pathlib import Path
from common.file_model.variant import Variant

class FileClient:
    """
    Client to load file in-memory
    """
    def __init__(self, config):
        self.data_root = config.get("data_root")
        self.collection = {}
        
    
    def get_variant_record(self, genome_uuid: str, variant_id: str, track_name: str) -> Variant:
        """
        Get a variant entry from variant_id
        """
        track_name = track_name or "dbSNP"
        base_file = vcfpy.Reader.from_path(os.path.join(self.data_root, genome_uuid, track_name, f"variation.vcf.gz"))
        base_record = self.fetch_record(base_file, variant_id)

        ## Initialise a variant object
        variant = Variant(base_record, base_file.header, genome_uuid)

        self.load_annotation_files(genome_uuid)
        for _, datafile in self.collection.items():
            record = self.fetch_record(datafile, variant_id)
            if record:
                variant.add_vcf_record(record, datafile.header)
        return variant
            
    def load_annotation_files(self, genome_uuid):

        root_dir = Path(os.path.join(self.data_root, genome_uuid))
        ## files can be vcf or text, need a parser class
        datafiles = root_dir.glob("annotation/frequencies/*.vcf.gz")
        for datafile in datafiles:
            self.collection[os.path.basename(datafile)] = vcfpy.Reader.from_path(datafile)
        
                                                           
    def fetch_record(self, datafile: str, variant_id: str):
        try:
            [contig, pos, id] = self.split_variant_id(variant_id)
            pos = int(pos)
            for rec in datafile.fetch(contig, pos-1, pos):
                if rec.ID[0] == id:
                    break # No duplicate variants
            return rec
        except:
            # Return None when variant cannot be fetched
            return

    def split_variant_id(self, variant_id: str):
        """
        Splits variant_id into separate fields
        """
        return variant_id.split(":")


