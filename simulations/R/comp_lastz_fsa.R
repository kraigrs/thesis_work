berlin_lastz <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin_lastz.txt",sep="\t");
berlin_fsa <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin_fsa.txt",sep="\t");
c1674_lastz <- read.table("/Users/kraigrs/Wittkopp/Graze/c1674_lastz.txt",sep="\t");
c1674_fsa <- read.table("/Users/kraigrs/Wittkopp/Graze/c1674_fsa.txt",sep="\t");

TUs <- merge(berlin_fsa,c1674_fsa,by.x="V1",by.y="V1");

nrow(subset(TUs,V3.x-V2.x > 0.05*(abs(V2.x-V2.y))));

nrow(subset(TUs,V3.y-V2.y > 0.05*(abs(V2.x-V2.y))));