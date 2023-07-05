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
from pysam import VariantFile
   
class FileClient:
    """
    Client to load file in-memory
    """
    def __init__(self, config, datafile=None):
        datafile = config.get("datafile")
        self.collection = VariantFile(datafile)
    
    def get_variant_record(self, variant_id: str):
        """
        Get a variant entry from variant_id
        """
        try: 
            [contig,pos,ref,id] = self.split_variant_id(variant_id)
            pos=int(pos)
        except:
            #TODO: This needs to go to thoas logger
            #TODO: Exception needs to be caught appropriately
            print("Something wrong")
        data = {}
        for rec in self.collection.fetch('1', pos-1, pos):
            data["name"] = rec.ref
        return data
        
        
    def split_variant_id(self, variant_id: str):
        """
        Splits variant_id into separate fields
        """
        return variant_id.split(":")


