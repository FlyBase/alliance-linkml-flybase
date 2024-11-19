[![Build Status](https://travis-ci.com/FlyBase/alliance-linkml-flybase.svg?branch=main)](https://travis-ci.com/FlyBase/alliance-linkml-flybase)

# alliance-linkml-flybase

<!-- toc -->

- [Overview](#Overview)
- [CodeManagement](#CodeManagement)
- [AllianceResources](#AllianceResources)
  * [LinkML](#LinkML)
  * [Credentials](#Credentials)
  * [PersistentStore](#PersistentStore)
  * [MoreInfo](#MoreInfo)
- [FileGeneration](#FileGeneration)
  * [LocalEnvironment](#LocalEnvironment)
  * [GoCDPipeline](#GoCDPipeline)
- [DataSubmission](#DataSubmission)
  * [FileUpload](#FileUpload)
  * [FinalCheck](#FinalCheck)
  * [References](#References)

## Overview
This code controls the export of FlyBase data to Alliance LinkML-compliant JSON
format. FlyBase and Alliance data types are represented by their own classes in
the `fb_datatypes.py` and `agr_datatypes.py` files, respectively.
There are various handlers dedicated to the transformation of specific
datatypes: e.g., `gene_handler.py`.
Specific scripts are then used to call various handlers. A script may call on
multiple handlers for the generation of a single export file: e.g., Alliance
AGMs encompass both FB strains and genotypes. A script may also output multiple
files: e.g., `AGM_data_retrieval_curation_gene.py` generates both gene and
allele-gene-association files.
The export files are then uploaded to the Alliance persistent store.

## CodeManagement
Modifications this repo should be made in new "release" branches corresponding to Alliance LinkML schema release branches: e.g., `v1.3.1`.  
These release branches should be used to spawn git releases/tags.  
These scripts are intended to be run in Docker using GoCD pipelines.
Files should be validated locally before uploading to the Alliance.
Use the [Alliance Curation API](https://curation.alliancegenome.org/swagger-ui/#/API%20Version/get_api_version) (details below) to find the latest supported LinkML version.  
- Write your code to this LinkML version.  
- Note that different objects may be have different version support as the persistent store develops.  

## AllianceResources
### LinkML
The Alliance LinkML-based data schema is held in the [agr_curation_schemas](https://github.com/alliance-genome/agr_curation_schema) repo.
### Credentials
Get an OKTA account (through Chris Grove). You will need this to interact with the persistent store.
### PersistentStore
The [Alliance Curation Site](https://curation.alliancegenome.org/#/) is the main hub.
Once logged into the curation site, you'll have access to persistent store data:
- View/update entries.
- View file loading.
- Get your personal `Curation API Token` - required for file uploads and API use.
You can also access the persistent store using the [Alliance Curation API](https://curation.alliancegenome.org/swagger-ui/#/API%20Version/get_api_version).
- You will need to log into the API with your personal `Curation API Token`.
- You can find the `Curation API Token` under your profile (top right corner) in the [Alliance Curation Site](https://curation.alliancegenome.org/#/)
### MoreInfo
See the Alliance A-Team [Confluence Site](https://agr-jira.atlassian.net/wiki/spaces/ATEAM/pages/476053508) for more info.

## FileGeneration
Use the `Alliance_LinkML_Submission` pipeline in the `Reporting_Build` pipeline group to generate, validate and upload files.  
With each run, files generated are stored locally in a directory within `/data/alliance/`.  
- The directory name will include the LinkML version and date (`MMMDD`), which correspond to the `LINKML_VERSION` and `ITERATION` pipeline values.  
### LocalEnvironment
1. Specify pipeline-specific variables.  
  - Sometimes the FB release used for LinkML submission will not be the latest release build.  
  - This happens when a release build has just been sent to FB-Indiana, but the public GFF file is not yet available for that release.  
  - For this reason, do not rely upon the `Reporting_Build` environment values (except for `ALLIANCETOKEN`).  
  - Specify these variables at the `Alliance_LinkML_Submission` pipeline level.  
    - `SERVER` - flysql machine where the reporting db is located: e.g., `flysql25`  
    - `DATABASE` - the name of the db to use (with `audit_chado`): e.g., `production_chado` or `fb_2022_05_reporting_audit`  
    - `RELEASE` - ensure that the release matches the db used: e.g., `2022_05` (reporting) or `2022_05_EP3` (production).  
    - `ITERATION` - the date `MMMDD` on which the pipeline is run.  
    - `LINKML_VERSION` - the LinkML version for the data: e.g., `v1.3.1`.  
### GoCDPipeline
The `Alliance_LinkML_Submission` pipeline automates these steps:  
1. Builds directory for data and log output in `/data/alliance/` folder.  
2. Gets HarvDev docker container.  
3. Builds docker image with FB `alliance-linkml-flybase` and Alliance `agr_curation_schemas` repos (the latter for validation).  
4. Fetches data using `alliance-linkml-flybase` scripts from a specified branch or release tag.  
5. Validates data files using the Alliance `validate_agr_schema.py` against the `jsonschema/allianceModel.schema` file.  
6. Sends e-mail that files are ready for review (size, contents), then pauses the pipeline.  
7. Has a separate, manually initiated phase for each file upload (due to dependencies).  

## DataSubmission
### FileUpload
As mentioned above, the GoCD pipeline will upload files to the Alliance.  
Details of how to upload files are [here](https://github.com/alliance-genome/agr_curation#submitting-data).  
- Your personal `API Curation Token` is required (available from your profile the [Alliance Curation Site](https://curation.alliancegenome.org/#/)).  
Expect an `OK` message upon completion of a successful file upload.
There are dependencies: e.g., the ALLELE file should be fully loaded (takes many hours) at the Alliance before uploading the disease file.
### ConfirmUploads
Go to the [Alliance Curation Site](https://curation.alliancegenome.org/#/).  
- From the left options, choose `Other Links` (at the bottom), then `Data Loads`.  
- For each data type uploaded, there should be a `Direct (LinkML) * Loads` selection on the right.  
- For each data type, all files ever uploaded are listed.  
- Click on the most recent one to see if it's been loaded (it can take hours).  
- Review any errors in loading any of the data objects - may require code or data fix.  
- Note - there can be load order dependencies: e.g., make sure all alleles are loaded before loading disease annotations.  
  - If you suspect that load order created problems, click on the pencil icon at the right to initiate a re-load of the file.
  - If you want to re-load a file, might be best to first check with the A-Team though.
  - There are dependencies on references too (see below).
### References
Loading of most LinkML data is dependent on presence of FB references at the Alliance.  
These references are in a separate database, the Alliance Bibliography Central (ABC), handled by the Alliance Blue Team.  
The pub schema is based on the distinct [agr_schemas](https://github.com/alliance-genome/agr_schemas/tree/master/ingest/resourcesAndReferences) model.  
As such, FB code of pub submission is in the related [alliance-flybase](https://github.com/FlyBase/alliance-flybase) repo.  
- See the README for this repo for details about submission and load issues.  
Submissions to the ABC are handled by the `Alliance_Pub_Export` pipeline (runs as part of the `Epicycle` pipeline group).  
