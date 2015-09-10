#!/bin/csh -f

R --slave < histGenes.R mel_sim_data/Hybrids/s_2_sequence/s_2_sequence.dmel.geneSummary.txt s_2_sequence.dmel >> ./trash &
R --slave < histGenes.R mel_sim_data/Hybrids/s_2_sequence/s_2_sequence.filtered.dmel.geneSummary.txt s_2_sequence.filtered.dmel >> ./trash &
#R --slave < histReads.R mel_sim_data/Hybrids/s_2_sequence/s_2_sequence.dmel.geneInfo.txt s_2_sequence.dmel >> ./trash &
#R --slave < histReads.R mel_sim_data/Hybrids/s_2_sequence/s_2_sequence.filtered.dmel.geneInfo.txt s_2_sequence.filtered.dmel >> ./trash &

R --slave < histGenes.R mel_sim_data/Hybrids/s_3_sequence/s_3_sequence.dmel.geneSummary.txt s_3_sequence.dmel >> ./trash &
R --slave < histGenes.R mel_sim_data/Hybrids/s_3_sequence/s_3_sequence.filtered.dmel.geneSummary.txt s_3_sequence.filtered.dmel >> ./trash &
#R --slave < histReads.R mel_sim_data/Hybrids/s_3_sequence/s_3_sequence.dmel.geneInfo.txt s_3_sequence.dmel >> ./trash &
#R --slave < histReads.R mel_sim_data/Hybrids/s_3_sequence/s_3_sequence.filtered.dmel.geneInfo.txt s_3_sequence.filtered.dmel >> ./trash &

R --slave < histGenes.R mel_sec_data/Hyb_nov.dmel.geneSummary.txt Hyb_nov.dmel >> ./trash &
R --slave < histGenes.R mel_sec_data/Hyb_nov.filtered.dmel.geneSummary.txt Hyb_nov.filtered.dmel >> ./trash &
#R --slave < histReads.R mel_sec_data/Hyb_nov.dmel.geneInfo.txt Hyb_nov.dmel >> ./trash &
#R --slave < histReads.R mel_sec_data/Hyb_nov.filtered.dmel.geneInfo.txt Hyb_nov.filtered.dmel >> ./trash &