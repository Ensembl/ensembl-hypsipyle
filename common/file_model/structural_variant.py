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

from typing import Any, Mapping, List, Union
import re
import os
import json


class StructuralVariant():
    variant_sources = {}  # used to cache source information, class attribute

    def __init__(self, record: Any, header: Any, genome_uuid: str) -> None:
        """Initialises a Variant instance.

        Args:
            record (Any): The variant record.
            header (Any): The header information.
            genome_uuid (str): The genome UUID.
        """
        self.genome_uuid = genome_uuid
        self.name = record.ID[0]
        self.record = record 
        self.header = header
        self.chromosome = record.CHROM         ###TODO: convert the contig name in the file to match the chromosome id given in the payload 
        self.position = record.POS
        self.info = record.INFO
        self.type = "StructuralVariant"
        self.vep_version = re.search("v\d+", self.header.get_lines("VEP")[0].value).group()
        self.population_map = {}