#!/bin/bash

# Define the directories and prefix


####################################################################################################################
####################################################################################################################
####################################################################################################################
source_dir="/storage/caishangLab/fangxiunan/project-rawdata/multiomic/singlecell/plate/20250313_enrichCRC/313-96CAT"
destination_dir="/storage/caishangLab/fangxiunan/project-ongoing/multiomic/singlecell/plate/20250313_enrichCRC/313-96CAT"
####################################################################################################################
####################################################################################################################
####################################################################################################################




exec > exec$(date +'%Y%m%d_%H%M').log 2>&1

scriptdir="/storage/caishangLab/fangxiunan/scripts/scplate"

mkdir -p $destination_dir
# Prepare the environment and generate the first job script

./s1.sh -s "$source_dir" -d "$destination_dir" -p "$scriptdir"

# Submit the first job and get the job ID

JOB_ID=$(sbatch job1.sh | awk '{print $4}')

# Check if the job ID was successfully obtained

if [ -z "$JOB_ID" ]; then

  echo "Failed to submit the first job."
  exit 1

fi

echo "First job submitted successfully with ID: $JOB_ID"

# Submit the second job with dependency on the first job

./s2.sh -d "$destination_dir" -p "$scriptdir" -j "$JOB_ID"

echo "Second job submitted successfully with dependency on job ID: $JOB_ID"
