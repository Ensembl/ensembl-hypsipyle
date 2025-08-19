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
from .utils import parse_ini, get_genome_uuids

# Open and read the JSON file
with open("seed-files/populations.json", "r") as file:
    data = json.load(file)

server = parse_ini("db.ini", "metadata")

population_map = {}

# For mouse, file list (for SNP and indels), output fields are per files
exception_list = ["MGP"]

# Convert common name to genome uuid
for species_name, species in data.items():
    population_frequencies = {}
    for population_source in species:
        if population_source["name"] != "UNSPECIFIED":
            pop_freq_name = population_source["name"]
            population_frequencies[pop_freq_name] = []
            for population_file in population_source["files"]:
                prefix = population_file["short_name"]
                for sub_population in population_file["include_fields"]:
                    for field_key, field_val in sub_population["fields"].items():
                        sub_population["fields"][field_key] = f"{prefix}_{field_val}"
                    # population can come from multiple files, skip duplicates
                    sub_population_item = next(
                        (
                            sub_pop_item
                            for sub_pop_item in population_frequencies[pop_freq_name]
                            if sub_pop_item["name"] == sub_population["name"]
                        ),
                        None,
                    )
                    if (
                        sub_population_item
                        and sub_population["name"] not in exception_list
                    ):
                        if sub_population_item != sub_population:
                            print(f"Conflicting fields for {sub_population['name']}")
                        continue
                    else:
                        population_frequencies[pop_freq_name].append(sub_population)
        else:
            for population_file in population_source["files"]:
                for sub_population in population_file["include_fields"]:
                    for field_key, field_val in sub_population["fields"].items():
                        sub_population["fields"][field_key] = f"UNSPECIFIED_{field_val}"
                        # population can come from multiple files, skip duplicates
                        if sub_population["name"] in population_frequencies:
                            if (
                                population_frequencies[sub_population["name"]]
                                != sub_population
                            ):
                                print(
                                    f"Conflicting fields for {sub_population['name']}"
                                )
                            continue
                        else:
                            population_frequencies[sub_population["name"]] = [
                                sub_population
                            ]

        genome_uuids = get_genome_uuids(server, species_name)
        for genome_uuid in genome_uuids:
            population_map[genome_uuid] = population_frequencies

# Write population-data.json
with open("test-data.json", "w") as file:
    json.dump(population_map, file, indent=4)
