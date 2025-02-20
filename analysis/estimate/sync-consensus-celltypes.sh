#!/bin/bash

# This script syncs consensus cell types from the OpenScPCA module
# Because of different credentials needed for getting these files, this
#  script should be run on its own and not from run-estimate-analysis.sh.
#
# Usage:
# bash sync-consensus-celltypes.sh <AWS profile name>

set -euo pipefail

aws_profile=$1

# Run script from its location
basedir=$(dirname "${BASH_SOURCE[0]}")
cd $basedir

consensus_dir="data/consensus-celltypes"
mkdir -p $consensus_dir

for project in SCPCP000001 SCPCP000002 SCPCP000006 SCPCP000009 SCPCP000017; do
  aws s3 sync s3://openscpca-nf-workflow-results/2024-11-25/cell-type-consensus/${project} ${consensus_dir}/${project} --profile ${aws_profile}
done
