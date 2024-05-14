FROM flybase/harvdev-docker:latest

WORKDIR /src

RUN mkdir /src/input/
RUN mkdir /src/logs/
RUN mkdir /src/output/
RUN mkdir /src/temp/

# This Dockerfile combines stuff from two repos.
# The ADD commands assume that repos are in specific folders in GoCD.

# 1. FlyBase scripts.
ADD git_alliance-linkml-flybase/requirements.txt                              /src/requirements.txt
ADD git_alliance-linkml-flybase/src/*.py                                      /src/

# 2. Alliance LinkML schema and validation script.
ADD git_agr_curation_schema/util/validate_agr_schema.py                       /src/validate_agr_schema.py
ADD git_agr_curation_schema/generated/jsonschema/allianceModel.schema.json    generated/jsonschema/allianceModel.schema.json

# Install required modules.
RUN pip3 install -r requirements.txt --no-cache-dir

ENTRYPOINT [ "/bin/bash" ]
