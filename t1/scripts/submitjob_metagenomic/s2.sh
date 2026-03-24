#!/bin/bash

# Parse command-line arguments

while getopts d:p:j: flag; do

  case "${flag}" in

    d) destination_dir=${OPTARG};;
    p) scriptdir=${OPTARG};;
    j) job_id=${OPTARG};;
  esac

done

# Generate the Slurm submission script for the second job

cat <<EOL > job2.sh
#!/bin/bash
#SBATCH -p intel-sc3,amd-ep2,intel-sc3-32c
#SBATCH -q normal
#SBATCH --cpus-per-task=20
#SBATCH --output=/storage/caishangLab/fangxiunan/log/exec_%A_%a.log
#SBATCH --error=/storage/caishangLab/fangxiunan/log/exec_%A_%a.err
#SBATCH --mem=100G
#SBATCH --dependency=afterok:$job_id

source /home/caishangLab/fangxiunan/anaconda3/etc/profile.d/conda.sh

conda activate metagenomics

scriptdir="${scriptdir}"
destination_dir="${destination_dir}"

find "\$destination_dir" -mindepth 3 -type f \( -name "*.dedupe.1.fq" -o -name "*.dedupe.2.fq" \) -exec cp {} "\$destination_dir/" \;

\${scriptdir}/step3_bacquantification.sh -d "\$destination_dir"


\${scriptdir}/step4_baclast.sh -d "\$destination_dir"

EOL

# Make the generated script executable

chmod +x job2.sh

# Submit the second job

sbatch job2.sh

if [ $? -ne 0 ]; then

  echo "Failed to submit the second job."
  exit 1

fi
