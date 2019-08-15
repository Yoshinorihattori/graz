#!/bin/bash

for i in $(ls *.eps); do epstopdf $i; done
pdfjam T.pdf dT.pdf qw.pdf vel.pdf Nu.pdf alpha.pdf --nup 2x3 --no-landscape --outfile merged.pdf
