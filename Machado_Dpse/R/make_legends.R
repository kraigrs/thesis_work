# create legends

pdf(file="/Users/kraigrs/Wittkopp/Machado_Dpse/reg_div_legend.pdf");

plot(0:10,rep(1,11),type='n',ylab="",xlab="",xaxt="n",yaxt="n",bty="n");
legend("center",
	legend=c("cis","trans","cis+trans","cisXtrans","compensatory","conserved","ambiguous"),
	fill=c("gray20","red","purple","green","orange","yellow","gray"),
	);

dev.off();

pdf(file="/Users/kraigrs/Wittkopp/Machado_Dpse/reg_div_legend.pdf");

plot(0:10,rep(1,11),type='n',ylab="",xlab="",xaxt="n",yaxt="n",bty="n");
legend("center",
	legend=c("cis","trans","cis & trans","conserved","ambiguous"),
	fill=c("blue","red","purple","black","gray"),
	);

dev.off();

pdf(file="/Users/kraigrs/Wittkopp/Machado_Dpse/MOI_legend.pdf");

plot(0:10,rep(1,11),type='n',ylab="",xlab="",xaxt="n",yaxt="n",bty="n");
legend("center",
	legend=c("similar","Dpse-dominant","Dbog-dominant","additive","underdominant","overdominant"),
	fill=c("cyan","green","magenta","orange","red","blue"),
	);
	
dev.off();