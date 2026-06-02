#!/bin/sh
# This script is used to create a environment.yml file for creating
# a conda environment (does not include pip)

conda env export --from-history  \
  | sed /prefix/d | sed 1d | sed 1i"name: OmegaCAnalysis" \
  > environment.yml

# This will optionally also add pip packages, but this does not seem to be very portable
# because the previous step can also add additional packages that would be listed here
# echo "  - pip:" >> environment.txt
# pip list | sed 1,2d | awk '{print "      - " $1 "==" $2}' >> environment.yml
#
#
# This will add the python dependencies that are specified in the custom requirements.txt file
if [ ! -f requirements.txt ]; then
  echo "Error: requirements.txt does not exist" 1>&2
else
  echo "  - pip:" >> environment.yml
  awk '{print "      - " $1}' requirements.txt >> environment.yml
fi


