#!/bin/bash
exec > hostalignment_parallel.log 2>&1

scriptdir="/mnt/resources/FXN/singlecell/sh_scripts/scplate"
source_dir="/mnt/resources/FXN/singlecell/plate_test"
destination_dir="/mnt/resources/FXN/singlecell/plate_test/0305_parallel"


$scriptdir/step1_hostalignment.sh -s $source_dir -d $destination_dir
$scriptdir/step2_bacquantification.sh -s $source_dir -d $destination_dir