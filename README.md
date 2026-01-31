# Biopython_Project
In-silico Analysis and Functional Prediction of an Unknown DNA-Encoded Protein from Escherichia col

Functional Annotation

Functional annotation was performed to predict the biological role of the selected open reading frame (ORF) identified from a bacterial genome. The longest ORF (2387 amino acids) was subjected to homology-based and domain-based annotation using established bioinformatics tools.

 Homology-Based Annotation (BLASTp)

Tool used: NCBI BLASTp (nr database)

Top hit: Putative adhesin

Organism: Escherichia coli

Sequence identity: ~99.6%

E-value: 0.0

The extremely low E-value and high sequence identity indicate strong evolutionary conservation and suggest that the predicted ORF encodes a protein with a function similar to known bacterial adhesins.

 Conserved Domain Analysis (InterProScan)

InterProScan analysis identified multiple conserved immunoglobulin-like (Ig-like) domains throughout the protein sequence.

InterPro ID: IPR013783

Domain name: Immunoglobulin-like fold

Database source: Gene3D

Domain distribution: Multiple repeated regions across the protein length

Statistical significance: E-values ranging from 10⁻¹³ to 10⁻¹⁸

Immunoglobulin-like domains are well known for their role in protein–protein interactions and are commonly found in bacterial surface and adhesion proteins.






Biological Interpretation

Biological interpretation was performed by integrating sequence quality analysis, ORF prediction, homology search, and conserved domain annotation.

The predicted protein is unusually large and contains multiple immunoglobulin-like domains, suggesting a modular architecture optimized for surface interaction. Such features are characteristic of bacterial adhesins, which mediate attachment to host tissues or abiotic surfaces.

Adhesins play a crucial role in:

Bacterial colonization

Biofilm formation

Host–pathogen interactions

The high degree of sequence conservation observed across Escherichia coli strains indicates that this protein likely performs an important and conserved biological function. Based on the combined computational evidence, the analyzed gene is predicted to encode a surface-associated adhesin involved in bacterial attachment and interaction processes.





Conclusion

This study demonstrates how an unknown DNA sequence can be functionally characterized using a systematic bioinformatics pipeline. Through ORF prediction, homology-based similarity analysis, and conserved domain identification, the target gene was successfully annotated as a putative bacterial adhesin. The repository provides a reproducible framework for computational functional prediction of uncharacterized genes.
