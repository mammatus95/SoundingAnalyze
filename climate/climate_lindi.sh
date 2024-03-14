#!/bin/bash
#SBATCH --job-name=climate             # Job name, will show up in squeue output
#SBATCH --time=0-12:59:59              # Runtime in DAYS-HH:MM:SS format
#SBATCH --partition=main              # On which Partition #main,calc,calclong, cip
#SBATCH --output=job1_%j.out           # File to which standard out will be written
#SBATCH --error=job1_%j.err            # File to which standard err will be written
#SBATCH --mail-type=ALL                # Type of email notification- BEGIN,END,FAIL,ALL


path_script=~/Dokumente/Programmierung/Radio/climate
path_zedat=public_html/temps/climate/image
python3=/home/mammatus95/Dokumente/miniconda3/bin/python3
cd ${path_script}

time ${python3} climateradio.py

scp -r *.png mammatus95@login.zedat.fu-berlin.de:${path_zedat}

#rm -rf *.png
