#!/bin/bash
touch rmp.tex
texname=tmp.tex
echo "\documentclass{jsarticle}
      \pagestyle{empty}
      \usepackage{amsmath,amssymb,bm,braket,physics}
      \usepackage[dvipdfmx]{graphicx}
      \usepackage{gnuplot-lua-tikz}
      \begin{document}
        \include {hamildiff_elp}
      \end{document}" > ${texname}

platex ${texname}
dviname=${texname/.tex/.dvi}
dvipdfmx $dviname
pdfname=${dviname/.dvi/.pdf}

#pdfcrop tmp.pdf
#mv tmp-crop.pdf al${al}p2_z2.pdf
mv tmp.pdf hamildiff_elp.pdf
rm tmp.*
rm *.tex
rm *.sty
