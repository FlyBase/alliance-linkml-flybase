[![Build Status](https://travis-ci.com/FlyBase/alliance-flybase.svg?token=7Nvc5gEdzuNraK13EL3s&branch=master)](https://travis-ci.com/FlyBase/alliance-flybase)

# alliance-flybase

<!-- toc -->

- [Overview](#Overview)
- [Credentials](#Credentials)
  * [S3Access](#S3Access)
  * [SubmissionAccessToken](#SubmissionAccessToken)
- [FileGeneration](#FileGeneration)
  * [ReleaseTag](#ReleaseTag)
  * [LocalEnvironment](#LocalEnvironment)
  * [GoCDPipeline](#GoCDPipeline)
- [DataSubmission](#DataSubmission)
  * [FileUpload](#FileUpload)
  * [FinalCheck](#FinalCheck)
- [PersistentStore](#PersistentStore)

## Overview
FlyBase in-house scripts repackage data from FlyBase for export to the Alliance project.  
Modifications should be made in new "release" branches corresponding to Alliance schema release branches: e.g., 1.0.1.0.  
These release branches should be used to spawn git releases/tags.  
Once an Alliance version is publicly released, merge the release branch into master.  
Scripts are run by a GoCD pipeline within docker. See [GoCDPipeline](#GoCDPipeline) section.  
Some scripts (not all yet) have a "-l" (local) or "-c" option (config) so that script can be run on local machine for development.  
Validation depends on schema in another Alliance repo, `agr_schemas`.  
Data is uploaded manually to AWS S3 bucket using curl.  

## Credentials
Data are stored in AWS S3 bucket. Need credentials to view the bucket.  
Need another a different API access token to upload files (via curl).  
  
### S3Access
This is a one time set-up.  
1. Obtain credentials to access Amazon S3 bucket for data dumps.  
  - Ask Gil or Chris for generic FlyBase developer credentials.  
  - Can obtain personal credentials from Adam Wright.  
2. Install AWS command line interface (CLI) software.  
```
bash
workon <python_virtual_env_of_your_choice>
pip install awscli
```  
3. Configure AWS user account (supply keys when prompted).  
  - Use the `aws configure` command.  
  - Supply keys when prompted.  
  - Note: "region name" is "us-east-1"  
4. Test credentials by listing files in the datadumps S3 bucket:  
```
bash
workon <python_virtual_env_of_your_choice>
aws s3 ls s3://mod-datadumps
```  
  - Note: you should see files with names like <mod>_<release>.tar.gz : e.g., `FB_1.0.4_4.tar.gz`  

### SubmissionAccessToken
You will need a separate API access token to submit files to the S3 bucket via `curl`.  
Typically each dev gets their own token - ask Olin Blodgett for one if needed.  

## FileGeneration

### ReleaseTag
Start any data submission by making a git release tag (a snapshot of the repo used for data generation).  
- View this repo, `alliance-flybase` in GitHub.  
- Go to "release" tab, and choose "Draft a new release".  
- Choose the appropriate release branch as target source for the new release.  
- Apply tag name corresponding to Alliance release number: e.g, "1.0.0.3".  
- Apply title like "Alliance 1.0.0.8".  

### LocalEnvironment
1. Check for the latest Panther orthology file (user: anonymous, password: password).  
```
cd /data/ortholog/panther/
ftp ftp.pantherdb.org
cd sequence_classifications/current_release/PANTHER_Sequence_Classification_files/
ls *fruit_fly*
get PTRH<release>_fruit_fly
```  
  - Save file here if it's new to FlyBase: `/data/ortholog/panther/`
  - If a new file was downloaded, update the filename in the `/src/AGR_data_retrieval_genes.py` script.  
  - *Eventually, code should be implemented to automate this step!*  
2. Double check the HTP metadata tags file available here:
  - `http://download.alliancegenome.org/<RELEASE>/HTP/TAGS/HTP_TAGS_0.json`
  - if new one is available, check to see that the `AGR_data_retrieval_datasets.py` script has right URL for the `process_agr_category_tags_file` function.
3. If running scripts locally:
  - Get updates and checkout appropriate branches for `alliance-flybase` and `agr_schemas` repos.  
  - Update the config file here: `/data/credentials/alliance/connection_info.cfg`.  
  - Move a python3 virtual environment and install requirements:  
```
cd <your_local_path>/alliance-flybase/
bash
workon <python_virtual_env_of_your_choice>
pip install -r requirements.txt
```  

### GoCDPipeline
Most of this process is now automated using a GoCD pipeline.  
See [Alliance_DQM_Submission](http://flysql22:8153/go/admin/pipelines/Alliance_DQM_Submission/general) in "Reporting_Build" pipeline group.  
See also general [GoCD documentation](https://github.com/FlyBase/harvdev-docs/blob/master/gocd/gocd.md).  

The pipeline automates these steps:  
1. Builds directory for data and log output in `/data/alliance/` folder.  
2. Gets HarvDev docker container.  
3. Builds docker image with FB `alliance-flybase` and Alliance `agr_schemas` repos.  
4. Fetches data using `alliance-flybase` scripts from a specified branch or release tag.  
5. Validates data files generated against Alliance schema (some specified branch) using `agr_validate.py`.   
6. Sends e-mail that files are ready for review (size, contents), then upload.

Some manual steps are still required before and after the pipeline runs:  
1. Download the most recent Panther file. See [LocalEnvironment](#LocalEnvironment) section above.  
2. Generate the appropriate `alliance-flybase` release git repo tag. See [ReleaseTag](#ReleaseTag) section above.  
3. Update `Materials`.  
  - For `alliance-flybase`, choose `master` in materials. During docker build, appropriate tag will be checked out.  
  - For `agr_schemas`, specify the appropriate release branch: e.g., release-1.0.1.5.  
4. Review the shared `Reporting_Build` environment.  
  - In the GoCD interface, look under `Admin`, then `Environments` for the `Reporting_Build` environment.  
  - This environment is shared by all reporting-related pipelines, including this `Alliance_DQM_Submission` pipeline.  
  - SERVER - flysql machine where the reporting db is located: e.g., `flysql25`  
  - DATABASE - the name of the reporting db to use: e.g., `fb_2022_02_reporting`  
  - ASSEMBLY - the Dmel reference genome assembly, currenty `R6`  
  - ALLIANCESCHEMA - e.g., `1.0.1.5`  
  - ALLIANCERELEASE - e.g., `5.2.0`  
  - TAG - the release tag for `flybase/alliance-flybase` Git repo that is to be used: e.g., `1.0.0.9-hotfix1`  
  - ITERATION - MMMDD label to use for output file folder in `data/alliance/`: e.g., `JUL22`  
  - RELEASE - e.g., `2021_02`. Only supply this if using an older release than is currently building for FB public release.  
  - ANNOTATIONRELEASE - e.g., `R6.46`.  
5. Update the specific `Environment variables` for this GoCD pipeline.  
  - Sometimes the FB release used for DQM submission will not be the latest release build.  
  - This happens when a release build has just been sent to FB-Indiana, but the public GFF file is not yet available for that release.  
  - In that case, specify pipeline-specific `Environment Variables` values as needed (these will override values for the `Reporting_Build` environment).  
6. Ensure files are of the correct size before upload (files lacking data can nonetheless validate OK).  
  - If code has changed, also worthwhile opening files to ensure no odd formatting.  
  - Particularly important for GFF files, for which there's no validation right now.  

## DataSubmission

### FileUpload
You will need the API access token. See [SubmissionAccessToken](#SubmissionAccessToken) section above.  
Each file must now be gzipped before upload (*.json, *.gff, *mitab.tsv).  
Upload the files, as per these [directions](https://github.com/alliance-genome/agr_chipmunk/blob/master/README.md).  
For each file, the upload endpoint is `ReleaseVersion_DataType_DataSubType/MOD`, linked by `=@` string to the path to the corresponding file to be uploaded.  
For example:  
```
curl \
-H "Authorization: Bearer <token>" \
-X POST "https://fms.alliancegenome.org/api/data/submit" \
-F "4.0.0_AGM_FB=@FB_1.0.1.4_agm.json.gz" \
-F "4.0.0_ALLELE_FB=@FB_1.0.1.4_allele.json.gz" \
-F "4.0.0_BGI_FB=@FB_1.0.1.4_BGI.json.gz" \
-F "4.0.0_CONSTRUCT_FB=@FB_1.0.1.4_construct.json.gz" \
-F "4.0.0_DAF_FB=@FB_1.0.1.4_disease.json.gz" \
-F "4.0.0_EXPRESSION_FB=@FB_1.0.1.4_expression.json.gz" \
-F "4.0.0_GFF_FB=@FB_1.0.1.4_GFF.gff.gz" \
-F "4.0.0_HTPDATASET_FB=@FB_1.0.1.4_dataset.json.gz" \
-F "4.0.0_HTPDATASAMPLE_FB=@FB_1.0.1.4_dataset_sample.json.gz" \
-F "4.0.0_HTVCF_FB=@dgrp.r6.vcf.gz" \
-F "4.0.0_INTERACTION-SOURCE_FB-GEN=@FB_1.0.1.4_genetic_interactions_mitab.tsv.gz" \
-F "4.0.0_INTERACTION-SOURCE_FB-MOL=@FB_1.0.1.4_physical_interactions_mitab.tsv.gz" \
-F "4.0.0_PHENOTYPE_FB=@FB_1.0.1.4_phenotype.json.gz" \
-F "4.0.0_REFERENCE_FB=@FB_1.0.1.4_reference.json.gz" \
-F "4.0.0_REF-EXCHANGE_FB=@FB_1.0.1.4_reference_exchange.json.gz" \
-F "4.0.0_RESOURCE_FB=@FB_1.0.1.4_resource.json.gz" \
-F "4.0.0_VARIATION_FB=@FB_1.0.1.4_variant.json.gz"
```

Success is indicated by this type of response:  
```
{"status":"success","message":null,"fileStatus": {"2.3.0_VARIATION_FB":"success"}}
```  

### FinalCheck
It's good to confirm that files landed in the s3 bucket as expected.  
One can use the AWS CLI to see recently uploaded files (within an appropriate python virtual environment).
Files are sorted first by release, then by datatype, then by mod/datasubtype. For example:  
```
aws s3 ls s3://mod-datadumps/4.0.0/VARIATION/FB/
```

## PersistentStore
DQMs have also started submitting files according to the `agr_curation_schema` LinkML model.  
Details on upload procedures are [here](https://github.com/alliance-genome/agr_curation#submitting-data).

Here's an overview:
Validate your JSON files with the [`util/validate_agr_schema.py`](https://github.com/alliance-genome/agr_curation_schema/blob/main/util/validate_agr_schema.py) in the `agr_curation_schema` repo.
- Be sure to validate against the appropriate release.

To get the API access key for file uploads, log into the requested A-Team curation site (one of the links below):
[ALPHA](https://alpha-curation.alliancegenome.org)
[BETA](https://beta-curation.alliancegenome.org)
[PROD](https://curation.alliancegenome.org)
- Then click on your profile in the upper right corner to find the key.

You will need to turn on your VPN access to the Alliance. instructions are [here](https://docs.google.com/document/d/1bp-8FP-4vcTOhExWPhlGZPfCa5eqkuIcpEvTT7lpCf4/edit).
- For me, this only works on my home computer, not on FlyBase workstations (Gil)

This is an example of the command to upload a GENE file.
```
curl -H "Authorization: Bearer ABC..." \
-X POST "https://curation.alliancegenome.org/api/data/submit" \
-F "GENE_FB=@gene_curation_fb_2022_02_reporting_audit.json"
```
