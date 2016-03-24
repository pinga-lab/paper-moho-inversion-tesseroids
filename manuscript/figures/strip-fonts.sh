#!/bin/bash

# Remove the fonts embedded in the EPS files. The journal doesn't like them.
# Inspired by the solutions in
# http://stackoverflow.com/questions/28797418/replace-all-font-glyphs-in-a-pdf-by-converting-them-to-outline-shapes
# http://tex.stackexchange.com/questions/24355/convert-from-pdf-to-eps-while-keeping-outlining-fonts

for f in *.eps
do
    epstool --copy --bbox $f tmp.eps
    gs -dNOPAUSE -dNOCACHE -dBATCH -sDEVICE=pdfwrite -dNoOutputFonts -dEPSCrop \
        -sOutputFile=tmp.pdf tmp.eps
    pdftops tmp.pdf -eps $f
    rm tmp.eps tmp.pdf
done
