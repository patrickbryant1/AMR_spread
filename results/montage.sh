#!/usr/bin/env bash

montage 'Acinetobacter spp.|Aminoglycosides.png' 'Acinetobacter spp.|Carbapenems.png' 'Acinetobacter spp.|Combined resistance (fluoroquinolones, aminoglycosides and carbapenems).png' 'Acinetobacter spp.|Fluoroquinolones.png' -tile 2x2 -geometry +2+2 act.png

montage 'Enterococcus faecalis|Aminopenicillins.png' 'Enterococcus faecalis|High-level gentamicin.png' 'Enterococcus faecalis|Vancomycin.png' 'Enterococcus faecium|Aminopenicillins.png' 'Enterococcus faecium|High-level gentamicin.png' 'Enterococcus faecium|Vancomycin.png' -tile 2x3 -geometry +2+2 ent.png

montage 'Escherichia coli|Aminoglycosides.png' 'Escherichia coli|Aminopenicillins.png' 'Escherichia coli|Carbapenems.png' 'Escherichia coli|Combined resistance (third-generation cephalosporin, fluoroquinolones and aminoglycoside).png' 'Escherichia coli|Fluoroquinolones.png' 'Escherichia coli|Third-generation cephalosporins.png' -tile 2x3 -geometry +2+2 esch.png

montage 'Klebsiella pneumoniae|Aminoglycosides.png' 'Klebsiella pneumoniae|Carbapenems.png' 'Klebsiella pneumoniae|Combined resistance (third-generation cephalosporin, fluoroquinolones and aminoglycoside).png' 'Klebsiella pneumoniae|Fluoroquinolones.png' 'Klebsiella pneumoniae|Third-generation cephalosporins.png' -tile 2x3 -geometry +2+2 kleb.png

montage 'Pseudomonas aeruginosa|Aminoglycosides.png' 'Pseudomonas aeruginosa|Carbapenems.png' 'Pseudomonas aeruginosa|Ceftazidime.png' 'Pseudomonas aeruginosa|Combined resistance (at least three of piperac. and tazob., fluoroq., ceftaz., aminogl. and carbapenems).png' 'Pseudomonas aeruginosa|Fluoroquinolones.png' 'Pseudomonas aeruginosa|PiperacillinTazobactam.png' -tile 2x3 -geometry +2+2 pseudo.png

#'Staphylococcus aureus|Meticillin (MRSA).png'

montage 'Streptococcus pneumoniae|Macrolides.png' 'Streptococcus pneumoniae|Penicillins.png' -tile 1x2 -geometry +2+2 strep.png
