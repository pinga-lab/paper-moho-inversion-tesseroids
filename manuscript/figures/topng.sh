#!/bin/bash

# Convert the EPS files to PNG for the initial submission. The EPS are too high
# res (especially the Moho models) and don't render properly on the journal
# system.
# use: bash topng.sh figure.eps
OUT=${1%.eps}.png
echo "Converting "$1 to $OUT
# Convert first to highres PNG because converting the Moho estimates to
# lowres is ugly. White lines appear between the model tiles.
convert -density 1200 $1 $OUT
# Now I can resample the highres PNG to a lower DPI without the ugly white
# lines.
convert -units PixelsPerInch $OUT -resample 300 $OUT
