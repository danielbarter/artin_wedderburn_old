#!/bin/bash

echo Dimension Time total_defect
ls data | xargs -I {} python3 run_jacobs_algebra.py data/{}

date
