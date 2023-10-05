#!/bin/bash

############################################################
# UPLOAD_FILE
#
# Upload objects to database. A workaround as I'm having
# issues with Python versions and installed modules.
#
############################################################

# Load Python version with necessary packages installed
module purge
ml Python/3.9.5-GCCcore-10.3.0

# Extract inputs
py_file=$1  # Python file to execute
sv_file=$2  # File name to save to
creds=$3    # Database credentials

# Call appropriate version of Python
python $py_file $sv_file $creds

