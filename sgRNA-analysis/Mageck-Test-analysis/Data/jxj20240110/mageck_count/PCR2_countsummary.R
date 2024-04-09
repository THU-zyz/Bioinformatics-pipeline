Sweave("PCR2_countsummary.Rnw");
library(tools);

texi2dvi("PCR2_countsummary.tex",pdf=TRUE);

