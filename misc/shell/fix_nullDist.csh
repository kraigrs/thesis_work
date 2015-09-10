#!/bin/csh -f

#foreach chr (2L 2LHet 2R 3R 3L 3LHet 4 U X XHet YHet dmel_mitochondrion_genome)

foreach chr (2L 2R 3R 3L 4 X)

set wd = /Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chr$chr

cd $wd

sed 's/,/_/g' zhr_z30.chr$chr\_k50.imprinting.txt > temp_gene.txt
sed 's/,/_/g' zhr_z30.chr$chr\_wind500000_step10000.imprinting.txt > temp_dist.txt

#cat zhr_z30.chr$chr\_k50.imprinting.txt | perl -pe 's/(coordinates)/\n$1/' > temp_gene.txt
#cat zhr_z30.chr$chr\_wind500000_step10000.imprinting.txt | perl -pe 's/(coordinates)/\n$1/' > temp_dist.txt

rm zhr_z30.chr$chr\_*

mv temp_dist.txt zhr_z30.chr$chr\_wind500000_step10000.imprinting.txt
mv temp_gene.txt zhr_z30.chr$chr\_k50.imprinting.txt

end
