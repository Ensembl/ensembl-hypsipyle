# GraphQL for Ensembl Variation

A [GraphQL trial](https://graphql.org/) for [Ensembl](https://www.ensembl.org) to reduce the need for RESTful web services.

This application is implemented with [Ariadne](https://ariadnegraphql.org/), a schema-first graphql framework for Python

GraphQL requires a schema (in /common) and implementation of resolver functions that know how to interpret specific parts of a GraphQL query. Resolvers are found in /resolver, and may also make use of "data loaders" to overcome inherent deficiencies in GraphQL implementations.


## Installation
Requires Python 3.8+.  

To install dependencies, run:

`pip install -r requirements.txt` for just the API.  Use this when deploying the service.

`pip install -r requirements-dev.txt` installs everything including dev dependencies like pytest, mypy etc.

## Running the API locally

The deployment assumes that we have a directory mounted with the convention `<data_root>/<genome_uuid>/variation.vcf.gz`
The code currently fetches a single vcf file within the genome_uuid folder. In the future, we plan to fetch data from multiple VCFs for a single genome_uuid.

Add path to the datafile in `./connections.conf`. 

The file follows the following template:
```
data_root = /usr/data/variation_vcfs
```

This command will start the server:

```uvicorn --workers 1 --host=0.0.0.0 graphql_service.server:APP```



## Containerisation for dev

Build the image using `./Dockerfile.dev`:

`docker build -t $NAME:$VERSION -f ./Dockerfile.dev .`

Run a container with the image (`--publish` below is exposing the container's ports to the host network):

`docker container run --publish 0.0.0.0:80:80/tcp --publish 0.0.0.0:8000:8000/tcp -ti -v <ensembl-hysipile-dir>:/app $NAME:$VERSION`


## Containerisation for prod
Build the image using `./Dockerfile.prod`:

`docker build -t $NAME:$VERSION  -f ./Dockerfile.prod .`

Run a container with the image (`--publish` below is exposing the container's ports to the host network):

`docker container run --publish 0.0.0.0:80:80/tcp --publish 0.0.0.0:8000:8000/tcp -ti $NAME:$VERSION`

The connection configuration is assumed to exist in the repo as the file `./connections.conf` and gets built into the Docker 
image. 
