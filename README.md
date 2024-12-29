# Alternative Splicing Analysis Pipeline

## Notes
- Ensure all required software and dependencies are installed before running the pipeline.
- Adjust paths and parameters based on your specific data and computational environment.

### Required Software and Dependencies
- **Operating System**: Linux Server (Ubuntu 20.04 server with more than 256 GB of RAM)
- **Genome and Annotation**: Version 106 of ENSEMBL for *Mus musculus*
- **Software Versions**:
  - Trim-galore: 0.6.4_dev
  - STAR: 2.7.3a
  - rMATS: v4_1_2
  - Sashimiplot: 2.0.4
  - vast-tools: v2.5.1

### Pipeline Summary

#### Folder and Data Preparation
1. **Create Directory Structure**: Organized folders for annotation, genome data, Salmon index, sample processing, trimming, and results.
2. **Download Genome Annotation and Sequence Files**:
   - GTF annotation file (*Mus musculus*, ENSEMBL version 106).
   - Primary assembly, cDNA, and ncRNA FASTA files.
3. **Prepare FASTA Files**:
   - Concatenate cDNA and ncRNA files.
4. **Transfer and Combine Raw FASTQ Files**:
   - Transfer raw sequencing files.
   - Combine replicates into unified files.
   - Clean up intermediate files.

#### Read Processing
5. **Trim Reads**: Use Trim-galore to preprocess paired-end FASTQ files and generate quality control reports.

#### Genome Indexing and Alignment
6. **Generate STAR Genome Index**:
   - Use STAR with the GTF annotation and genome FASTA.
   - Set `sjdbOverhang` to 149 (read length - 1).
7. **Adjust Open Files Limit**: Set `ulimit -n` to 3076 for handling multiple files.

#### Differential Splicing Analysis
8. **Run rMATS**:
   - Set up Python 2.7 environment.
   - Perform differential splicing analysis using rMATS.
   - Input paired-end RNA-Seq data, with 150 bp read length and a 5% statistical significance cutoff.
   - Specify paths for wild-type and condition samples.

#### Visualization
9. **Generate Sashimi Plots**: Use rmats2sashimiplot to visualize alternative splicing events (SE, MXE, A5SS, A3SS, RI).

#### Alternative Splicing Analysis with vast-tools
10. **Align Reads**: Align reads for each sample using vast-tools.
11. **Combine Results**: Merge alignment results across all samples.
12. **Differential Inclusion Levels**:
    - Compare inclusion levels between conditions.
    - Generate tidy outputs for better readability.
13. **Visualization**:
    - Create plots for differential alternative splicing results.
    - Compare gene expression levels between conditions.
14. **Differential Splicing Thresholds**:
    - Perform differential splicing analysis at multiple thresholds.
    - Tidy outputs for downstream analysis.

#### Final Outputs
- Comprehensive differential splicing results and tidy tables.
- Sashimi plots for visualization.
- Gene expression comparisons between conditions.

---
For questions or troubleshooting, contact the pipeline maintainer or consult the respective software documentation.
