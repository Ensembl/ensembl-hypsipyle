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
import configparser
import subprocess

def parse_ini(ini_file: str, section: str = "database") -> dict:
    config = configparser.ConfigParser()
    config.read(ini_file)
    
    if not section in config:
        print(f"[ERROR] Could not find '{section}' config in ini file - {ini_file}")
        exit(1)
    else:
        host = config[section]["host"]
        port = config[section]["port"]
        user = config[section]["user"]
        database = config[section]["database"]
    
    
    return {
        "host": host, 
        "port": port, 
        "user": user,
        "database": database
    }

def get_genome_uuids(server: dict, production_name: str) -> str:
    if production_name == "homo_sapiens_gca\\d{9}v\\d{1}":
        query = f"SELECT genome_uuid FROM genome WHERE PRODUCTION_NAME LIKE \"homo_sapiens_gca%\";"
    else:
        query = f"SELECT genome_uuid FROM genome WHERE PRODUCTION_NAME LIKE \"{production_name}\" ORDER BY GENEBUILD_DATE DESC LIMIT 1;"
    process = subprocess.run(["mysql",
            "--host", server["host"],
            "--port", server["port"],
            "--user", server["user"],
            "--database", server["database"],
            "-N",
            "--execute", query
        ],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    genome_uuids = process.stdout.decode().strip().split()
    return genome_uuids


# Open and read the JSON file
with open('seed-files/population_metadata.json', 'r') as file:
    data = json.load(file)

server = parse_ini("db.ini","metadata")

population_map = {}

# Convert common name to genome uuid
for species_name, species in data.items():
    
    genome_uuids = get_genome_uuids(server,species_name)
    for genome_uuid in genome_uuids:
        population_map[genome_uuid] = species

# Write population-data.json
with open('test-metadata.json', 'w') as file:
    json.dump(population_map, file, indent=4)  



