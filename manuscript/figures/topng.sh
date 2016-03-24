#!/bin/bash

# Convert the EPS files to PNG for the initial submission. The EPS are too high
# res (especially the Moho models) and don't render properly on the journal
# system.
# use: bash topng.sh DPI

for f in *.eps
do
    convert -density 400 $f "$(basename "$f" .eps).png"
done
convert -density 800 south-america-moho.eps south-america-moho.png
