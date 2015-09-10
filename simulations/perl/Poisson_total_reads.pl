# define expression values for each gene
open(EXP,"$expression") or die "\nError opening $expression\n";
while(<EXP>)    
{
  chomp;
  @elements = split("\t",$_); 
  @parts = split("_",$elements[0]);
  $gene = $parts[0];
  $reg = $parts[1];
  $level = $elements[1];
  unless($level < 20){$expr{$gene}{$reg} = $level;}
}
close EXP;
