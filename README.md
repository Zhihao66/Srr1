# An evolutionarily ancient transcription factor drives spore morphogenesis in mushroom-forming fungi
[DOI: 10.1016/j.cub.2025.02.025](https://doi.org/10.1016/j.cub.2025.02.025)

### RBH_Phylostratygraphy_250210.R 
Code for phylostratygraphic analysis
### Srr1_RNA-Seq_analysis.R
Code for RNA-seq analysis
### all_counts.tsv
Raw count of RNA-seq data
### org.CABnewInterGO1.eg.db.zip
OrgDB for GO enrichment
org.CABnewInterGO1.eg.db is build base on the GO terms linked from interpro annotation(version5.57-90.0,database date:20240207) with '-goterms' parameter, see STAR method in the paper

see file Copci_IPR_anno_20240207.tsv

Installation:

unzip the file

``
install.packages("./org.CABnewInterGO1.eg.db", repos=NULL,type="source")
``

``
library(org.CABnewInterGO1.eg.db)
``
### Copci_IPR_anno_20240207.tsv: 
IPR annotatation of new protein model of Coprinopsis cinerea (AmutBmut strain), see more information of the new genome and protein model here:

[Unraveling Morphogenesis, Starvation, and Light Responses in a Mushroom-Forming Fungus, Coprinopsis cinerea, Using Long Read Sequencing and Extensive Expression Profiling ](https://doi.org/10.1101/2024.05.10.593147)
 and [mushroomDB](https://mushroomdb.brc.hu/)
