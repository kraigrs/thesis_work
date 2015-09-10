# define expression values for each gene
open(EXP,"$expression") or die "\nError opening $expression\n";
while(<EXP>)    
{
  chomp;
  @elements = split("\t",$_); 
  @parts = split("_",$elements[0]);
  $gene = $parts[0];
  $reg = $parts[1];
  $ref_level = $elements[1];
  $alt_level = $elements[2];
  unless($ref_level + $alt_level < 20)
  {
    $expr{$gene}{$reg} = [$ref_level,$alt_level];
    #print "$gene\t$reg\n";
  }
}
close EXP;
