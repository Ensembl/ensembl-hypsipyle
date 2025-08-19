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
import csv
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
    
    
    return {
        "host": host, 
        "port": port, 
        "user": user
    }


def fetch_populations(server: dict, variation_db_name: str) -> str:
    query = f"SELECT name, description, collection, population_id from population;"
    process = subprocess.run(["mysql",
            "--host", server["host"],
            "--port", server["port"],
            "--user", server["user"],
            "--database",variation_db_name,
            "-N",
            "--execute", query
        ],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    populations = process.stdout.decode().strip()
    return populations

def fetch_super_population(server: dict, variation_db_name: str, pop_id: str):
    query = f"SELECT name FROM population WHERE population_id = \
        (SELECT ps.super_population_id FROM population_structure ps, population p  \
        WHERE p.population_id = ps.sub_population_id AND p.population_id = {pop_id})"
    process = subprocess.run(["mysql",
            "--host", server["host"],
            "--port", server["port"],
            "--user", server["user"],
            "--database",variation_db_name,
            "-N",
            "--execute", query
        ],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    super_population = process.stdout.decode().strip()
    return super_population

def fetch_sub_population(server: dict, variation_db_name: str, pop_id: str):
    query = f"SELECT name FROM population WHERE population_id IN \
        (SELECT ps.sub_population_id FROM population_structure ps, population p  \
        WHERE p.population_id = ps.super_population_id AND p.population_id = {pop_id})"
    process = subprocess.run(["mysql",
            "--host", server["host"],
            "--port", server["port"],
            "--user", server["user"],
            "--database",variation_db_name,
            "-N",
            "--execute", query
        ],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    sub_populations = process.stdout.decode().strip()
    return sub_populations

server = parse_ini("db.ini","variation")  
with open("assemblies.txt") as fd:
    rd = csv.reader(fd, delimiter="\t", quotechar='"')
    for row in rd:
        populations=[]
        populations_data= fetch_populations(server, row[0]).split("\n")
        for pop in populations_data:
            pop = pop.split("\t")
            p_map={}
            super_population = fetch_super_population(server, row[0], pop[3])
            sub_populations = fetch_sub_population(server, row[0], pop[3])
            p_map["name"] = pop[0]
            p_map["description"] = pop[0]
            p_map["type"] = "regional"
            p_map["is_global"] = True
            p_map["display_group_name"] = pop[0]
            p_map["super_population"] = {"name" : super_population }  if super_population else None
            if super_population:
                p_map["is_global"] = False
            p_map["sub_populations"] = []
            for sub_pop in sub_populations.split("\n"):
                if sub_pop:
                    p_map["sub_populations"].append ({"name": sub_pop})
            populations.append(p_map)


        # Write population-data.json
        with open(f"{row[0]}.json", 'w') as file:
            json.dump(populations, file, indent=4)
