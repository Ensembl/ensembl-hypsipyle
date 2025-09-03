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

import json
from utils import parse_ini, get_genome_uuids

# Open and read the JSON file
with open("seed-files/population_metadata.json", "r") as file:
    data = json.load(file)

server = parse_ini("db.ini", "metadata")

population_map = {}

# Convert common name to genome uuid
for species_name, species in data.items():
    genome_uuids = get_genome_uuids(server, species_name)
    for genome_uuid in genome_uuids:
        population_map[genome_uuid] = species

# Write population-data.json
with open("test-metadata.json", "w") as file:
    json.dump(population_map, file, indent=4)
