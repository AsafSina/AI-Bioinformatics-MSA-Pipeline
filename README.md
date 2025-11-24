# AI-Bioinformatics-MSA-Pipeline
Automated Multiple Sequence Alignment (MSA) pipeline using Python to orchestrate T-Coffee, MAFFT, and Clustal Omega, with quantitative quality assessment via Shannon Entropy.

**FOCUS:** Automated tool orchestration and quantitative quality assessment for **AI/Computational Biology**.

This repository features a modular Python pipeline for the rapid comparative analysis of Multiple Sequence Alignment (MSA) tools, demonstrating skills essential for data science in biotech.

## Pipeline Capabilities

* **Tool Orchestration:** Integrates and executes external MSA tools (`T-Coffee`, `MAFFT`, `Clustal Omega`) via Python.
* **Quantitative Evaluation:** Calculates key quality metrics like **Average Shannon Entropy** (measure of conservation) to objectively assess alignment quality.
* **Comparative Analysis:** Benchmarks tool performance for data-driven decisions (speed vs. accuracy trade-off).

## Results and Analysis

The pipeline successfully identified that for the tested dataset, the **T-Coffee** alignment yielded the highest conservation (**Entropy: 0.5036**).

The full analytical write-up, including the comparison table and conclusions, is available here: **[analysis_report.md](analysis_report.md)**

##  Installation & Usage

1.  **Prerequisites:** Python 3.x and the command-line tools: `t_coffee`, `mafft`, `clustalo`.
2.  **Execution:**
    ```bash
    git clone [https://github.com/YOUR_USERNAME/AI-Bioinformatics-MSA-Pipeline.git](https://github.com/YOUR_USERNAME/AI-Bioinformatics-MSA-Pipeline.git)
    cd AI-Bioinformatics-MSA-Pipeline
    python msa_full_pipeline.py
    ```


