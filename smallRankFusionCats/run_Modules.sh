#!/bin/bash

echo Dimension Time total_defect
ls moduleData/*.txt | xargs -I {} python3 run_jacobs_algebra.py moduleData/{}

date
