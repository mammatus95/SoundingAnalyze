#!/bin/bash

for X in 10035 10113 10184 10238 10304 10393 10410 10548 10618 10771 10739 10868 06610 16045 11520 11747 11035
do
    python3 readhtml_bufr.py $(date +%Y-%m-%d) 00 ${X}
    ./auswertung $(wc -l ${X}.txt) $(date +%m) $(date +%d) 00

    #rm -rf *.png *.txt *.html ${X}thermo.png
    echo "${X} done it"
done

#rm -f *.png *.txt *.html
