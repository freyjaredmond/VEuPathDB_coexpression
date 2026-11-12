#!/bin/bash -l
#SBATCH -J alphafold3
#SBATCH -o alphafold3_%A_%a.out
#SBATCH -e alphafold3_%A_%a.err
#SBATCH -p gpu
#SBATCH --ntasks-per-node=1
#SBATCH -c 32
#SBATCH --gres=gpu:1
#SBATCH --mem=128G
#SBATCH -t 3-00:00:00
#SBATCH --array=1-N%3          

# AlphaFold3 pairwise predictions (SLURM array job)
# Input:
#   - JSON_DIR:   folder of AlphaFold3 input JSONs, one per protein pair
#   - JSON_LIST:  text file listing the JSON file names, one per line (array task i runs line i)
#   - MODEL_DIR:  AlphaFold3 model parameters 
#   - DB_DIR:     AlphaFold3 public databases
#   - AF3_SIF:    AlphaFold3 container image
#
# Output:
#   - OUTDIR/<job name>/   AlphaFold3 output for each pair (structures, confidence and PAE files)
#
# actifpTM is calculated from these outputs afterwards (see run_reactifptm.py)

# ---- paths: change these for your system ----
JSON_DIR=/path/to/af3_input/jsons
JSON_LIST=/path/to/AF3_interactors_jsons.txt
MODEL_DIR=/path/to/alphafold3/model_parameters
DB_DIR=/path/to/alphafold3_databases/public_databases
AF3_SIF=/path/to/alphafold3.sif
OUTDIR=/path/to/af3_output
# ----------------------------------------------

mkdir -p ${OUTDIR}

# separate tmp directory for each array task, so tasks don't overwrite each other's files
TMPDIR_TASK=${OUTDIR}/tmp/${SLURM_ARRAY_TASK_ID}
mkdir -p ${TMPDIR_TASK}

# JSON file for this array task
FILE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${JSON_LIST})

# memory settings let the model use unified (GPU + CPU) memory for large pairs
singularity exec \
     --nv \
     --bind ${JSON_DIR}:/root/af_input \
     --bind ${OUTDIR}:/root/af_output \
     --bind ${MODEL_DIR}:/root/models \
     --bind ${DB_DIR}:/root/public_databases \
     --bind ${TMPDIR_TASK}:/tmp \
     --env XLA_PYTHON_CLIENT_PREALLOCATE=false,TF_FORCE_UNIFIED_MEMORY=true,XLA_CLIENT_MEM_FRACTION=3.2 \
     ${AF3_SIF} \
     python /app/alphafold/run_alphafold.py \
     --json_path=/root/af_input/${FILE} \
     --model_dir=/root/models \
     --db_dir=/root/public_databases \
     --pdb_database_path=/root/public_databases/mmcif_files \
     --output_dir=/root/af_output
