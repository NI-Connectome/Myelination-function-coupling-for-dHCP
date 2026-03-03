This folder contains the code required to associate MFC and gene expression, as well as data and result. The gene expression data is sourced from BrainSpan Atlas of the Developing Human Brain (www.brainspan.org).

Specifically,

the program 'Transcriptomic_association_analysis.m' uses the PLS regression to select genes which have significantly correlated spatial pattern with MFC. Then, PLS1+ genes and PLS1- genes are used to GO analysis separately.

the subfolder 'result' includes selected genes (OrderedTerms_PLS_Loadings_Pfdr001_prenatal.csv), GO result for PLS+ genes, as well as GO result for PLS- genes.
