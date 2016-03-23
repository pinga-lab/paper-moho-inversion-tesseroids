#!/bin/bash

# Inspired by the solutions in 
# http://stackoverflow.com/questions/28797418/replace-all-font-glyphs-in-a-pdf-by-converting-them-to-outline-shapes
# http://tex.stackexchange.com/questions/24355/convert-from-pdf-to-eps-while-keeping-outlining-fonts

for f in *.eps
do
    epstool --copy --bbox $f tmp.eps
    gs -dNOPAUSE -dNOCACHE -dBATCH -sDEVICE=pdfwrite -dNoOutputFonts -dEPSCrop \
        -sOutputFile="$(basename "$f" .eps).pdf" tmp.eps    
    pdftops "$(basename "$f" .eps).pdf" -eps $f
    convert -density 300 $f "$(basename "$f" .eps).png"
    rm tmp.eps "$(basename "$f" .eps).pdf"
done
