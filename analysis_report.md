
Comparative Analysis of MSA Tools for Uniprot Sequences

I. Project Overview

This project demonstrates the development of a modular and automated bioinformatics pipeline written in Python. The goal is to perform Multiple Sequence Alignment (MSA) using industry-standard tools and quantitatively assess the quality of the resulting alignments using key metrics like **Shannon Entropy**.

The pipeline automatically integrates three different alignment methodologies: **MAFFT (Fast and Scalable)**, **Clustal Omega (General Purpose)**, and **T-Coffee (High Accuracy/Consistency)**.

II. Methodology
* **Input Data:** A set of protein sequences downloaded from Uniprot were used in FASTA format.
* **Pipeline:** The custom `msa_full_pipeline.py` script was executed using the 'Run All Aligners' option to automate sequential execution and scoring.
* **Scoring Metrics:**
    * **Gap Percentage (%):** Measures the sparsity/completeness of the alignment.
    * **Match Percentage (%):** Measures the prevalence of the most common residue at each column.
    * **Average Shannon Entropy:** A fundamental measure of sequence conservation. **Lower entropy indicates higher conservation** (less variability) across the alignment.

III. Results: Comparative Analysis Summary

COMPARATIVE ANALYSIS SUMMARY
======================================================================
Alignment File                 |    Gap (%) |  Match (%) |    Entropy | Tool
----------------------------------------------------------------------
uniprot_tcoffee.fasta          |     90.921 |     91.334 |     0.5036 | Tcoffee | Highest Accuracy/Consistency |
uniprot_mafft.fasta            |     84.212 |     85.320 |     0.7989 | Mafft | Fast and Balanced |
uniprot_clustal.fasta          |     84.721 |     85.448 |     0.7932 | Clustal | General Purpose / Fast |
----------------------------------------------------------------------


IV. Conclusion and AI/Tech Insight
The quantitative analysis clearly indicates a significant difference in conservation quality based on the algorithm used.

The **T-Coffee** alignment yielded the most conserved result, demonstrated by the lowest Average Shannon Entropy value of **0.5036**. This is consistent with its methodology, which prioritizes sequence consistency even at the cost of high gap insertion (90.921% Gap). MAFFT and Clustal Omega provided faster, yet less conserved alignments (Entropy approx. 0.79).

This project demonstrates proficiency in:
1.  **Tool Orchestration:** Successfully automating the execution and output handling of external bioinformatics software (T-Coffee, MAFFT, ClustalO) via a Python pipeline.
2.  **Quantitative Evaluation:** Applying information theory metrics (Shannon Entropy) for data-driven quality assessment.
3.  **Data-Driven Decision Making:** Providing a comparative framework to select the optimal alignment tool (T-Coffee for high-accuracy studies, MAFFT/Clustal for large, fast screening) based on the required balance between speed and quality.
