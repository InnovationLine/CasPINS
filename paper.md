# CasPINS: An Integrated Web-Based Platform for CRISPR/TALEN gRNA Design, Primer Generation, and Indel Decomposition Analysis

Rinki Dasgupta¹,²,*, Kaushik Das³

¹ Department of Psychiatry and Neuroscience, Dell Medical Center at University of Texas at Austin, Austin, TX, 78701, USA

² Department of Biology, Texas Woman's University, Denton, TX, 76204, USA

³ Department of Computer Applications, National Institute of Technology Jharkhand, Jamshedpur, 831014, India

*Corresponding author

---

## ABSTRACT

Genome editing technologies, particularly CRISPR-Cas systems [[1]](https://doi.org/10.1126/science.1258096) [[2]](https://doi.org/10.1126/science.1225829) and TALENs [[3]](https://doi.org/10.1038/nbt.1755), have revolutionized biological research, yet researchers must navigate multiple disconnected software platforms to design guide RNAs [[4]](https://doi.org/10.1038/nbt.3437), generate validation primers [[5]](https://doi.org/10.1093/nar/gks596), and analyze editing outcomes [[6]](https://doi.org/10.1093/nar/gku936)—a fragmented workflow that impedes discovery and limits accessibility. Here, we present CasPINS (Cas-Primer-Indel Suite), an open-source integrated platform that unifies the complete genome editing workflow within a single interactive environment. CasPINS represents a significant advance over existing tools by providing exceptional capabilities that accelerate biological insights: support for 90+ species across diverse taxa; 14 CRISPR-Cas variants [[7]](https://doi.org/10.1016/j.cell.2014.05.010) alongside comprehensive TALEN design with RVD sequence generation and off-target analysis [[8]](https://doi.org/10.1038/nrm3486); direct NCBI/Ensembl database integration for automated sequence retrieval; and intelligent nested PCR primer design (PCR I for genomic amplification, PCR II for Sanger sequencing) positioned relative to predicted cut sites [[9]](https://doi.org/10.1534/g3.114.016527). Building upon established open-source tools including Primer3 [[5]](https://doi.org/10.1093/nar/gks596) and BioPython, CasPINS extends their capabilities through seamless workflow integration. The platform incorporates AB1 trace decomposition using Non-Negative Least Squares (NNLS) [[10]](https://doi.org/10.1137/1.9781611971217)—an open, reproducible algorithm similar to TIDE [[6]](https://doi.org/10.1093/nar/gku936)—to quantify indel spectra and editing efficiency, supporting both single-sample and unlimited batch analysis for high-throughput screening [[11]](https://doi.org/10.1101/gr.244293.118). By consolidating target design, primer engineering, and outcome analysis into an intuitive graphical interface requiring no programming expertise, CasPINS democratizes genome editing technology and substantially reduces technical barriers. CasPINS is freely available under the MIT license at https://github.com/raju1stnov/CasPINS.

---

## ABBREVIATIONS

AB1, Applied Biosystems chromatogram file format; bp, base pair; Cas, CRISPR-associated protein; CasPINS, Cas-Primer-Indel Suite; cDNA, complementary DNA; CFD, cutting frequency determination; CRISPR, Clustered Regularly Interspaced Short Palindromic Repeats; DSB, double-strand break; FokI, Flavobacterium okeanokoites restriction enzyme I; GC, guanine-cytosine; gRNA, guide RNA; GUI, graphical user interface; HDR, homology-directed repair; ICE, Inference of CRISPR Edits; indel, insertion and/or deletion; MM, mismatch; mRNA, messenger RNA; NCBI, National Center for Biotechnology Information; NGS, next-generation sequencing; NHEJ, non-homologous end joining; NNLS, Non-Negative Least Squares; nt, nucleotide; OT, off-target; PAM, protospacer adjacent motif; PCR, polymerase chain reaction; QC, quality control; RefSeq, Reference Sequence; RFLP, restriction fragment length polymorphism; RVD, Repeat Variable Diresidue; sgRNA, single guide RNA; SpCas9, Streptococcus pyogenes Cas9; SaCas9, Staphylococcus aureus Cas9; T7E1, T7 endonuclease I; TALE, Transcription Activator-Like Effector; TALEN, Transcription Activator-Like Effector Nuclease; TIDE, Tracking of Indels by Decomposition; Tm, melting temperature; WT, wild-type.

---

## INTRODUCTION

The advent of programmable nucleases has fundamentally transformed molecular biology, enabling precise genome modification across virtually all model organisms and cell types [[1]](https://doi.org/10.1126/science.1258096) [[7]](https://doi.org/10.1016/j.cell.2014.05.010) [[12]](https://doi.org/10.1146/annurev-biochem-060815-014607). CRISPR-Cas systems, derived from bacterial adaptive immunity, have emerged as the predominant technology due to their simplicity, efficiency, and versatility [[2]](https://doi.org/10.1126/science.1225829) [[13]](https://doi.org/10.1126/science.1231143). Concurrently, TALENs (Transcription Activator-Like Effector Nucleases) remain valuable for applications requiring high specificity or targeting PAM-restricted genomic regions [[3]](https://doi.org/10.1038/nbt.1755) [[8]](https://doi.org/10.1038/nrm3486). Both technologies generate targeted double-strand breaks that are repaired through error-prone non-homologous end joining (NHEJ), producing insertions and deletions (indels) that can disrupt gene function [[14]](https://doi.org/10.1146/annurev-genet-110410-132435). The widespread adoption of these technologies across biomedical research, agriculture, and biotechnology has created an urgent need for accessible, comprehensive computational tools that can support the entire experimental workflow without requiring specialized programming expertise.

A complete genome editing workflow encompasses several distinct computational and experimental phases. First, researchers must identify optimal guide RNA (gRNA) sequences for CRISPR-Cas or design appropriate TALE binding domains for TALENs, considering factors including on-target efficiency, off-target potential, and genomic context [[15]](https://doi.org/10.1038/nbt.3026) [[4]](https://doi.org/10.1038/nbt.3437) [[16]](https://doi.org/10.1038/nbt.2647). Second, PCR primers must be designed to amplify regions flanking the predicted cut site for downstream validation by Sanger sequencing or next-generation sequencing [[9]](https://doi.org/10.1534/g3.114.016527). Third, sequencing results must be analyzed to quantify editing efficiency and characterize the indel spectrum [[6]](https://doi.org/10.1093/nar/gku936) [[11]](https://doi.org/10.1101/gr.244293.118). Each phase requires specialized computational approaches, and the quality of execution at each step directly impacts experimental success and the reliability of biological conclusions.

Currently, this workflow requires researchers to transition between multiple independent tools, each with distinct interfaces, input requirements, and output formats. gRNA design platforms such as CHOPCHOP [[17]](https://doi.org/10.1093/nar/gkz365), CRISPOR [[18]](https://doi.org/10.1093/nar/gky354), and Benchling provide sophisticated target identification with extensive scoring metrics, yet lack integrated primer design capabilities. Primer design tools like Primer3 [[5]](https://doi.org/10.1093/nar/gks596) and NCBI Primer-BLAST operate independently of editing context and require manual coordination with nuclease target sites. Indel analysis software including TIDE [[6]](https://doi.org/10.1093/nar/gku936), ICE [[19]](https://doi.org/10.1101/251082), and ampliCan [[11]](https://doi.org/10.1101/gr.244293.118) require separately prepared input files and operate in isolation from the design phase. This fragmentation creates significant inefficiencies, increases error potential during data transfer between platforms, limits experimental reproducibility, and presents a substantial barrier for researchers without dedicated bioinformatics support. Moreover, the lack of workflow integration impedes the iterative optimization that is often essential for successful genome editing experiments.

The open-source software ecosystem has been instrumental in democratizing computational biology, with tools such as Primer3 for oligonucleotide design, BioPython for sequence analysis, and various CRISPR design utilities providing foundational capabilities to the research community. However, these tools operate independently, and their integration requires programming expertise that many experimental biologists lack. There exists a critical need for unified platforms that build upon these established resources while providing seamless workflow integration accessible to non-computational researchers.

Here, we present CasPINS (Cas-Primer-Indel Suite), an open-source unified platform that integrates all computational steps of the genome editing workflow into a single, accessible environment. CasPINS provides comprehensive gRNA design for 14 CRISPR-Cas variants with multiple scoring algorithms across 90+ species, full TALEN design with RVD sequence generation, off-target analysis, and restriction site mapping, automated primer design with direct NCBI and Ensembl database integration for both genomic DNA and cDNA/mRNA templates, and AB1 trace decomposition using Non-Negative Least Squares (NNLS) for indel quantification with single-sample and unlimited batch analysis modes. By extending the capabilities of established open-source tools within an intuitive graphical interface, CasPINS represents a significant advance that enables researchers across all experimental backgrounds to efficiently execute genome editing experiments from initial design through final validation, accelerating biological discovery while ensuring reproducibility and accessibility.

---

## METHODS

### Software Architecture and Implementation

CasPINS is implemented in Python 3.8+ and deployed as a web-based application using the Streamlit framework, enabling interactive visualization and real-time user feedback without requiring client-side installation. The application architecture follows a modular design pattern with four core modules: gRNA Design, TALEN Design, Primer Design, and Indel Analysis. Each module operates independently while sharing common data structures and utility functions to ensure workflow continuity.

### Database Integration

Sequence retrieval is accomplished through programmatic interfaces to Ensembl REST API and NCBI E-utilities. For Ensembl queries, the application constructs REST requests to retrieve genomic sequences, gene annotations, and transcript information using species-specific endpoints. NCBI integration utilizes the Entrez Programming Utilities for RefSeq sequence retrieval and gene annotation. All database queries implement error handling with retry logic to accommodate transient network failures.

### gRNA Design Algorithm

The gRNA identification algorithm scans input sequences for PAM motifs specific to the selected Cas variant. For each potential target site, the upstream 20-nucleotide sequence (or appropriate length for non-SpCas9 variants) is extracted and evaluated against multiple criteria. GC content is calculated as the proportion of guanine and cytosine bases. Homopolymer detection identifies consecutive identical nucleotides exceeding user-specified thresholds. Self-complementarity analysis evaluates potential secondary structure formation using thermodynamic calculations.

Efficiency scoring implements three established algorithms. The Doench 2016 score [[4]](https://doi.org/10.1038/nbt.3437) applies position-specific nucleotide weights derived from large-scale screening data. The Moreno-Mateos CRISPRscan score [[20]](https://doi.org/10.1038/nmeth.3543) emphasizes features predictive of in vivo activity in zebrafish. The Xu score [[21]](https://doi.org/10.1101/gr.191452.115) incorporates sequence features correlated with cutting efficiency in human cells. A weighted composite score combines individual algorithm outputs with user-adjustable weighting factors.

### TALEN Design Algorithm

TALEN pair identification scans genomic sequences for compatible binding sites meeting established design criteria [[3]](https://doi.org/10.1038/nbt.1755) [[8]](https://doi.org/10.1038/nrm3486). The algorithm requires 5' thymine residues for both TALE binding domains, evaluates arm lengths between 15-20 bp with preference for 17-18 bp, and constrains spacer regions to 12-21 bp with optimal range of 14-16 bp. RVD sequences are generated using standard codon assignments: NI for adenine, HD for cytosine, NG for thymine, and NN for guanine.

Off-target analysis estimates potential binding at genomic sites with 0-3 mismatches using sequence alignment algorithms. Restriction enzyme site identification within spacer regions facilitates RFLP-based screening of edited clones by scanning for recognition sequences of common restriction enzymes.

### Primer Design Algorithm

Primer design leverages Primer3 [[5]](https://doi.org/10.1093/nar/gks596) as the core oligonucleotide design engine, with CasPINS-specific wrappers that position primers relative to predicted nuclease cut sites. PCR I primers are designed with 500+ bp flanking regions to generate 800-2000 bp amplicons suitable for T7E1 assays and preparative applications. PCR II primers are positioned to place the cut site 150-300 bp from the sequencing primer binding site, optimizing Sanger trace quality in the analysis window. Design parameters include target Tm of 60°C ± 2°C, GC content of 40-60%, and filters for secondary structure avoidance.

### Indel Analysis Algorithm

The indel quantification algorithm implements Non-Negative Least Squares (NNLS) decomposition [[10]](https://doi.org/10.1137/1.9781611971217) to model edited chromatograms as linear combinations of shifted control traces. For control trace C(x) and edited trace E(x), the decomposition solves:

E(x) = Σ αᵢ × C(x - δᵢ)

where αᵢ represents the fractional contribution of indel size δᵢ (ranging from -10 to +10 bp), subject to constraints Σ αᵢ = 1 and αᵢ ≥ 0. The optimization is performed using scipy.optimize.nnls. Editing efficiency is calculated as (1 - α₀) × 100%, where α₀ represents the wild-type fraction.

AB1 trace file parsing extracts fluorescence intensity values for all four channels (A, C, G, T) along with base calls and quality scores. Trace alignment between control and edited samples uses the gRNA sequence to identify the expected cut position. Quality metrics include R² goodness-of-fit and confidence assessment based on trace quality in the decomposition window.

### Batch Processing

Batch analysis mode implements parallel processing of multiple edited samples against a common control. The system automatically detects edited AB1 files within organized project directories based on naming conventions. Progress tracking provides real-time status updates with estimated completion times. Results aggregation generates summary statistics including mean editing efficiency, standard deviation, and outlier identification.

---

## RESULTS

### Overview of CasPINS Architecture

CasPINS is implemented as a web-based application built on the Streamlit framework, providing an interactive graphical interface accessible through any modern web browser without requiring local software installation or programming expertise. The platform comprises four integrated modules, gRNA Design [[4]](https://doi.org/10.1038/nbt.3437), TALEN Design [[3]](https://doi.org/10.1038/nbt.1755), Primer Design [[5]](https://doi.org/10.1093/nar/gks596), and Indel Analysis [[6]](https://doi.org/10.1093/nar/gku936). Each module operates independently but shares a common data architecture, enabling seamless transition between workflow stages. User-specified data directories maintain organized project structures, and all intermediate and final outputs are preserved for reproducibility.

### gRNA Design Module, Multi-Species Multi-Cas Support

The gRNA Design module supports target identification for CRISPR-Cas genome editing across 90+ species spanning mammals, birds, fish, amphibians, reptiles, invertebrates, plants, fungi, and bacteria. Users provide target gene identifiers through multiple input formats, gene symbols, Ensembl gene IDs, RefSeq accession numbers, genomic coordinates, or direct sequence input. Species selection automatically configures appropriate genome assemblies, with support for both current and legacy reference assemblies such as GRCh38/hg38 and GRCh37/hg19 for human.

CasPINS implements 14 Cas nuclease variants with their respective PAM requirements. Cas9 variants include SpCas9 (NGG), SpCas9-NG (NG), SpCas9-VQR (NGAN), SpCas9-EQR (NGAG), SaCas9 (NNGRRT), SaCas9-KKH (NNNRRT), NmeCas9 (NNNNGATT), CjCas9 (NNNNRYAC), eSpCas9, SpCas9-HF1, and HypaCas9 as enhanced specificity variants. Cas12 variants include Cas12a/Cpf1 (TTTV), AsCas12a, LbCas12a, FnCas12a, and enAsCas12a (TTYN/VTTV).

For each potential gRNA, CasPINS calculates three established efficiency prediction scores, Doench et al. [[4]](https://doi.org/10.1038/nbt.3437), Moreno-Mateos et al. CRISPRscan [[20]](https://doi.org/10.1038/nmeth.3543), and Xu et al. [[21]](https://doi.org/10.1101/gr.191452.115). These scores incorporate position-specific nucleotide preferences, GC content, di-nucleotide features, and self-complementarity penalties. A weighted composite score combines individual metrics, with adjustable weighting for specific experimental contexts.

Results are presented in an interactive table displaying gRNA sequence, PAM, genomic position, strand orientation, GC content, and all efficiency scores. Users can filter results by score thresholds, GC content range with default 40-60%, maximum homopolymer length, and targeting of specific exon regions where early exons are recommended for knockout applications. Selected gRNAs are exported to project directories for downstream primer design and indel analysis.

### TALEN Design Module, Paired Nuclease Engineering

For researchers requiring TALEN-based editing, particularly in contexts where CRISPR PAM requirements are limiting or higher specificity is desired, CasPINS provides comprehensive TALEN pair design capabilities. The TALEN module accepts the same input formats as gRNA design and identifies optimal paired binding sites within target genes.

TALEN pairs are identified according to established design principles [[3]](https://doi.org/10.1038/nbt.1755) [[8]](https://doi.org/10.1038/nrm3486), each arm requires a 5' thymine (T) for optimal N-terminal domain binding, arm lengths of 15-20 bp with optimal 17-18 bp, and spacer regions of 12-21 bp with optimal 14-16 bp for effective FokI dimerization. For each identified target site, CasPINS generates TALE 1 (Left arm) consisting of target sequence, RVD (Repeat Variable Diresidue) sequence with standard code where NI=A, HD=C, NG=T, NN=G, length, and GC content; TALE 2 (Right arm) consisting of complementary strand target with corresponding RVD array; and Spacer region consisting of sequence and length of the intervening region containing the predicted cut site.

The module displays RVD sequences in both linear format and a grid visualization showing the tandem repeat array architecture. Off-target analysis estimates potential binding at genomic sites with 0-3 mismatches, presented in a found/total format such as 0/50 indicating zero off-target sites among 50 potential matching loci.

Restriction enzyme site analysis identifies recognition sequences within the spacer region, facilitating restriction fragment length polymorphism (RFLP)-based screening of edited clones. Results include enzyme names and cut positions for 40+ common restriction enzymes. A composite efficiency score ranks TALEN pairs based on arm length optimization, GC content, RVD composition, and specificity metrics.

### Primer Design Module, Database-Integrated PCR Engineering

The Primer Design module addresses a critical gap in existing workflows by generating CRISPR/TALEN-aware primer sets with direct integration to NCBI and Ensembl sequence databases [[5]](https://doi.org/10.1093/nar/gks596). This module designs two nested primer pairs optimized for distinct experimental objectives.

PCR I Primers for Genomic DNA Amplification target genomic DNA sequences flanking the predicted nuclease cut site with 500+ bp flanking regions. These primers generate 800-2000 bp amplicons suitable for T7 endonuclease I (T7E1) mismatch assays [[9]](https://doi.org/10.1534/g3.114.016527), RFLP screening, and preparative amplification for cloning or sequencing library preparation. Design parameters include optimal Tm of 60°C ± 2°C, GC content of 40-60%, and avoidance of secondary structures.

PCR II Primers for Sanger Sequencing are designed for direct Sanger sequencing of edited populations, positioned to place the cut site 150-300 bp from the sequencing primer binding site, which is optimal for trace quality in the critical analysis window. These primers can target either genomic DNA or mRNA/cDNA, enabling detection of both genomic modifications and expression-level effects. Amplicon sizes of 400-800 bp balance sequencing quality with sufficient flanking context.

Sequence acquisition occurs through three pathways, direct database query to Ensembl or NCBI RefSeq using gene identifiers, user-provided local sequence files, or previously stored project sequences. The module fetches genomic coordinates corresponding to user-specified gRNA sequences, automatically identifying cut site positions and orienting primer design accordingly.

Results are presented with comprehensive primer analytics including sequence, length, Tm, GC content, 3' stability, hairpin potential, and predicted amplicon size. A detailed design report documents all parameters and provides formatted sequences suitable for direct oligonucleotide ordering.

### Indel Analysis Module, Sanger Trace Decomposition

The Indel Analysis module implements trace decomposition for quantifying CRISPR/TALEN-mediated editing outcomes from Sanger sequencing chromatograms. This module processes AB1 trace files from control (unedited) and experimental (edited) samples to determine editing efficiency and characterize the indel spectrum.

The algorithmic approach employs Non-Negative Least Squares (NNLS) decomposition [[10]](https://doi.org/10.1137/1.9781611971217), a standard signal processing technique that models the edited chromatogram as a linear combination of shifted control traces. For a given control trace C(x), the edited trace E(x) is decomposed as E(x) = Σ αᵢ × C(x - δᵢ), where αᵢ represents the fractional contribution of each indel size δᵢ ranging from -10 to +10 bp, subject to the constraint Σ αᵢ = 1 and αᵢ ≥ 0. The wild-type fraction where δ = 0 directly indicates the unedited population, with editing efficiency calculated as (1 - α₀) × 100%.

For Single Sample Analysis, users upload paired AB1 files consisting of control.ab1 and edited.ab1 and specify the target gene from their project directory. CasPINS automatically retrieves the associated gRNA sequence to determine the expected cut position. The analysis generates sequence alignment visualization with base-by-base comparison of control and edited sequences with mismatch highlighting, chromatogram overlay showing normalized trace profiles with signal quality and the onset of trace divergence, indel spectrum histogram displaying distribution of insertion and deletion sizes with percentage contributions, and summary metrics including editing efficiency, dominant indel size, R² quality score, and confidence assessment.

For Batch Analysis Mode supporting high-throughput applications, CasPINS processes unlimited samples against a common control within organized project directories. Users select target genes from their data directory, and CasPINS automatically identifies all edited AB1 files matching the naming convention. Progress tracking displays real-time analysis status with estimated completion time.

Batch results include summary visualizations with pie charts showing edited versus wild-type proportions across all samples, comparative indel spectrum plots, and exportable tabular data. Statistical summaries report average editing efficiency, standard deviation, and identification of outlier samples requiring manual review.

---

## DISCUSSION

CasPINS represents a comprehensive solution for genome editing experimental design and analysis, consolidating functionalities that previously required multiple independent software platforms. By integrating gRNA/TALEN design, primer engineering, and indel quantification within a unified interface, CasPINS substantially reduces the technical overhead and potential for error inherent in multi-tool workflows.

### Unified Workflow Advantages

The principal advantage of CasPINS is workflow continuity. A researcher can progress from target gene identification through publication-ready editing efficiency data without exporting data files, reformatting inputs, or learning multiple software interfaces. This integration is particularly valuable for laboratories without dedicated bioinformatics support, enabling efficient experimental execution by researchers with diverse computational backgrounds.

The shared data architecture ensures consistency across workflow stages. gRNA sequences selected during target design are automatically available for primer positioning during the subsequent module. Sample organization established during project setup propagates to indel analysis, eliminating filename management errors. This coherence reduces the troubleshooting burden that often accompanies multi-platform analyses.

### Comprehensive Species and System Support

CasPINS supports an extensive organism range, over 90 species from diverse taxonomic groups, enabling genome editing projects across fundamental research, agricultural, and biomedical applications. The inclusion of 14 Cas variants addresses the expanding CRISPR toolkit, accommodating both standard SpCas9 workflows and specialized applications requiring alternative PAM specificities or enhanced fidelity variants.

The integration of TALEN design within the same platform provides flexibility for targets where CRISPR systems face limitations. Researchers can evaluate both nuclease classes for a given target, selecting the optimal approach based on PAM availability, predicted specificity, and prior laboratory experience.

### Database Integration for Error Reduction

Direct integration with NCBI and Ensembl databases addresses a common source of experimental failure, sequence discrepancies between design and target. By retrieving reference sequences programmatically, CasPINS ensures primer design against validated genomic coordinates rather than potentially outdated local sequence files. This integration also simplifies primer design for researchers unfamiliar with database query interfaces.

### Accessible Indel Quantification

The NNLS-based trace decomposition algorithm implemented in CasPINS provides robust indel quantification from standard Sanger sequencing data, the most widely accessible sequencing modality. While next-generation sequencing approaches offer higher resolution, Sanger sequencing remains the practical choice for initial screening and routine editing confirmation in most laboratories. The batch analysis capability enables systematic evaluation of multiple clones or experimental conditions, supporting both small-scale pilot experiments and larger screening efforts.

### Comparison with Existing Tools

Several excellent tools address individual components of the genome editing workflow. CHOPCHOP [[17]](https://doi.org/10.1093/nar/gkz365), CRISPOR [[18]](https://doi.org/10.1093/nar/gky354), and other gRNA design platforms provide sophisticated target identification with extensive scoring metrics and off-target analysis. TIDE [[6]](https://doi.org/10.1093/nar/gku936), ICE [[19]](https://doi.org/10.1101/251082), and ampliCan [[11]](https://doi.org/10.1101/gr.244293.118) offer specialized indel analysis with various algorithmic approaches. Primer3 [[5]](https://doi.org/10.1093/nar/gks596) remains the standard for general primer design.

CasPINS does not seek to replace these specialized tools but rather to provide an integrated alternative that prioritizes workflow efficiency over maximum feature depth in any single module. Researchers requiring advanced capabilities, such as comprehensive off-target profiling against whole-genome alignments or maximum-likelihood indel decomposition, may benefit from supplementing CasPINS with dedicated tools. However, for the majority of standard genome editing projects, CasPINS provides sufficient capability across all workflow stages within a single, accessible interface.

### Future Development

We anticipate several enhancements to CasPINS based on evolving genome editing practices. Integration of base editing and prime editing design rules would extend coverage to precision editing modalities beyond nuclease-induced indels. Incorporation of Cas13 for RNA targeting would support transcriptome-level interventions. Enhanced off-target analysis through integration of precomputed whole-genome alignments would strengthen specificity predictions. Finally, cloud deployment options would enable access without local installation requirements.

---

## SOFTWARE AVAILABILITY

CasPINS is freely available as open-source software at https://github.com/raju1stnov/CasPINS under the MIT license. The platform runs locally through Python installation or via containerized Docker deployment. Comprehensive documentation, tutorial workflows, and example datasets are provided within the repository.

---

## DATA AVAILABILITY

All data generated or analyzed during this study are included in this published article.

---

## CODE AVAILABILITY

The CasPINS software is open-source and is publicly available on GitHub at [link to be provided]. The specific version used for this publication has been archived on Zenodo [DOI to be provided].

---

## FUNDING

No funding received.

---

## AUTHOR CONTRIBUTIONS

Rinki Dasgupta conceived the project, conceptualized the development of the software, and wrote the manuscript. Kaushik Das contributed to the development of software, algorithmic design, and manuscript review.

---

## ACKNOWLEDGEMENTS

The authors thank the developers of open-source tools including Primer3, BioPython, and Streamlit, upon which CasPINS builds. We acknowledge the publicly available genome databases from NCBI and Ensembl that enable CasPINS functionality.

**AI Disclosure:** In accordance with COPE guidelines, we disclose that artificial intelligence tools (Claude, Anthropic) were used for manuscript drafting assistance, formatting references, and preparing figure legends. All AI-generated content was thoroughly reviewed, verified, and edited by the authors to ensure accuracy and scientific validity. The authors take full responsibility for the content of this manuscript. AI tools were not used for data generation, analysis, or interpretation of results.

---

## COMPETING INTERESTS

The authors declare no competing interests. As an independent, open-source initiative developed without dedicated institutional or grant funding, CasPINS leverages publicly available datasets from open scientific repositories for its validation and test cases.

---

## REFERENCES

1. Doudna, J. A., & Charpentier, E. (2014). The new frontier of genome engineering with CRISPR-Cas9. *Science*, 346(6213), 1258096. [https://doi.org/10.1126/science.1258096](https://doi.org/10.1126/science.1258096)
2. Jinek, M., Chylinski, K., Fonfara, I., Hauer, M., Doudna, J. A., & Charpentier, E. (2012). A programmable dual-RNA-guided DNA endonuclease in adaptive bacterial immunity. *Science*, 337(6096), 816-821. [https://doi.org/10.1126/science.1225829](https://doi.org/10.1126/science.1225829)
3. Miller, J. C., Tan, S., Qiao, G., Barber, K. A., Rebar, E. J., Babiarz, J. E., Isik, I., Minber, A. L., Piatt, D. P., Hamelin, A. M., Guschin, D. Y., & Gregory, P. D. (2011). A TALE nuclease architecture for efficient genome editing. *Nature Biotechnology*, 29(2), 143-148. [https://doi.org/10.1038/nbt.1755](https://doi.org/10.1038/nbt.1755)
4. Doench, J. G., Fusi, N., Sullender, M., Hegde, M., Vaimberg, E. W., Donovan, K. F., Smith, I., Tothova, Z., Wilen, C., Orchard, R., Virgin, H. W., Listgarten, J., & Root, D. E. (2016). Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. *Nature Biotechnology*, 34(2), 184-191. [https://doi.org/10.1038/nbt.3437](https://doi.org/10.1038/nbt.3437)
5. Untergasser, A., Cutcutache, I., Koressaar, T., Ye, J., Faircloth, B. C., Remm, M., & Rozen, S. G. (2012). Primer3—new capabilities and interfaces. *Nucleic Acids Research*, 40(15), e115. [https://doi.org/10.1093/nar/gks596](https://doi.org/10.1093/nar/gks596)
6. Brinkman, E. K., Chen, T., Amendola, M., & van Steensel, B. (2014). Easy quantitative assessment of genome editing by sequence trace decomposition. *Nucleic Acids Research*, 42(22), e168. [https://doi.org/10.1093/nar/gku936](https://doi.org/10.1093/nar/gku936)
7. Hsu, P. D., Lander, E. S., & Zhang, F. (2014). Development and applications of CRISPR-Cas9 for genome engineering. *Cell*, 157(6), 1262-1278. [https://doi.org/10.1016/j.cell.2014.05.010](https://doi.org/10.1016/j.cell.2014.05.010)
8. Joung, J. K., & Sander, J. D. (2013). TALENs: A widely applicable technology for targeted genome editing. *Nature Reviews Molecular Cell Biology*, 14(1), 49-55. [https://doi.org/10.1038/nrm3486](https://doi.org/10.1038/nrm3486)
9. Vouillot, L., Thélie, A., & Pollet, N. (2015). Comparison of T7E1 and surveyor mismatch cleavage assays to detect mutations triggered by engineered nucleases. *G3: Genes, Genomes, Genetics*, 5(3), 407-415. [https://doi.org/10.1534/g3.114.016527](https://doi.org/10.1534/g3.114.016527)
10. Lawson, C. L., & Hanson, R. J. (1995). *Solving Least Squares Problems*. Society for Industrial and Applied Mathematics. [https://doi.org/10.1137/1.9781611971217](https://doi.org/10.1137/1.9781611971217)
11. Labun, K., Guo, X., Chavez, A., Church, G., Gagnon, J. A., & Valen, E. (2019). Accurate analysis of genuine CRISPR editing events with ampliCan. *Genome Research*, 29(5), 843-847. [https://doi.org/10.1101/gr.244293.118](https://doi.org/10.1101/gr.244293.118)
12. Wang, H., La Russa, M., & Qi, L. S. (2016). CRISPR/Cas9 in genome editing and beyond. *Annual Review of Biochemistry*, 85, 227-264. [https://doi.org/10.1146/annurev-biochem-060815-014607](https://doi.org/10.1146/annurev-biochem-060815-014607)
13. Cong, L., Ran, F. A., Cox, D., Lin, S., Barretto, R., Habib, N., Hsu, P. D., Wu, X., Jiang, W., Marraffini, L. A., & Zhang, F. (2013). Multiplex genome engineering using CRISPR/Cas systems. *Science*, 339(6121), 819-823. [https://doi.org/10.1126/science.1231143](https://doi.org/10.1126/science.1231143)
14. Symington, L. S., & Gautier, J. (2011). Double-strand break end resection and repair pathway choice. *Annual Review of Genetics*, 45, 247-271. [https://doi.org/10.1146/annurev-genet-110410-132435](https://doi.org/10.1146/annurev-genet-110410-132435)
15. Doench, J. G., Hartenian, E., Graham, D. B., Tothova, Z., Hegde, M., Smith, I., Sullender, M., Ebert, B. L., Xavier, R. J., & Root, D. E. (2014). Rational design of highly active sgRNAs for CRISPR-Cas9-mediated gene inactivation. *Nature Biotechnology*, 32(12), 1262-1267. [https://doi.org/10.1038/nbt.3026](https://doi.org/10.1038/nbt.3026)
16. Hsu, P. D., Scott, D. A., Weinstein, J. A., Ran, F. A., Konermann, S., Agarwala, V., Li, Y., Fine, E. J., Wu, X., Shalem, O., Cradick, T. J., Marraffini, L. A., Bao, G., & Zhang, F. (2013). DNA targeting specificity of RNA-guided Cas9 nucleases. *Nature Biotechnology*, 31(9), 827-832. [https://doi.org/10.1038/nbt.2647](https://doi.org/10.1038/nbt.2647)
17. Labun, K., Montague, T. G., Krause, M., Torres Cleuren, Y. N., Tjeldnes, H., & Valen, E. (2019). CHOPCHOP v3: Expanding the CRISPR web toolbox beyond genome editing. *Nucleic Acids Research*, 47(W1), W171-W174. [https://doi.org/10.1093/nar/gkz365](https://doi.org/10.1093/nar/gkz365)
18. Concordet, J. P., & Haeussler, M. (2018). CRISPOR: Intuitive guide selection for CRISPR/Cas9 genome editing experiments and screens. *Nucleic Acids Research*, 46(W1), W242-W245. [https://doi.org/10.1093/nar/gky354](https://doi.org/10.1093/nar/gky354)
19. Hsiau, T., Conant, D., Rosber, N., Maures, T., Waite, K., Yang, J., Kelber, S., Sternberg, S., Nguyen, M., Stoner, R., Derr, A., Hinkley, S., Bhattacharyya, R., & Bhattacharya, S. (2019). Inference of CRISPR edits from Sanger trace data. *bioRxiv*. [https://doi.org/10.1101/251082](https://doi.org/10.1101/251082)
20. Moreno-Mateos, M. A., Vejnar, C. E., Beaudoin, J. D., Fernandez, J. P., Mis, E. K., Khokha, M. K., & Giraldez, A. J. (2015). CRISPRscan: Designing highly efficient sgRNAs for CRISPR-Cas9 targeting in vivo. *Nature Methods*, 12(10), 982-988. [https://doi.org/10.1038/nmeth.3543](https://doi.org/10.1038/nmeth.3543)
21. Xu, H., Xiao, T., Chen, C. H., Li, W., Meyer, C. A., Wu, Q., Wu, D., Cong, L., Zhang, F., Liu, J. S., Brown, M., & Liu, X. S. (2015). Sequence determinants of improved CRISPR sgRNA design. *Genome Research*, 25(8), 1147-1157. [https://doi.org/10.1101/gr.191452.115](https://doi.org/10.1101/gr.191452.115)

---

## FIGURE LEGENDS

Figure 1. CasPINS gRNA Design Module. (A) Input interface showing target gene selection including Gene Symbol, Ensembl ID, RefSeq ID, Genomic Coordinates, or Paste Sequence, species selection with 90+ organisms, CRISPR-Cas system selection with 14 Cas variants, and comprehensive sequence filters including GC content range, maximum homopolymer length, and target region preferences. (B) gRNA Design Results displaying ranked guide sequences with comprehensive scoring using Doench 2016, Moreno-Mateos, and Xu algorithms with genomic locations. Selected gRNA details panel shows sequence information, PAM, genomic coordinates, and individual efficiency scores.

Figure 2. CasPINS TALEN Design Module. (A) TALEN input interface with gene symbol entry and TALEN-specific parameters. (B) TALEN Pair Results showing designed pairs with TALE1 (left arm) and TALE2 (right arm) sequences, RVD arrays displayed in grid format, spacer region details, off-target analysis in found/total format covering MM0-MM3, restriction enzyme sites within the spacer, and composite efficiency scores. Color-coded target sequences distinguish TALE1 binding, spacer, and TALE2 binding regions.

Figure 3. CasPINS Primer Design Module. (A) Primer Design interface with database integration options for NCBI and Ensembl, sequence type selection for Genomic DNA or cDNA/mRNA, and gRNA sequence input for cut site positioning. (B) Enhanced Primer Design Results showing PCR I for genomic amplification and PCR II for Sanger sequencing primer pairs with comprehensive analytics including Tm, GC content, amplicon size, and formatted ordering sequences. Sequence information panel displays fetched genomic coordinates and gene metadata.

Figure 4. CasPINS Single Sample Indel Analysis. (A) Indel Spectrum histogram showing distribution of insertions as positive values in blue and deletions as negative values in red with wild-type fraction in green. (B) Complete analysis visualization including sequence alignment with mismatch highlighting, chromatogram traces for control and edited samples with cut site indication, indel spectrum, and pie chart showing edited versus wild-type proportions. Summary metrics display editing efficiency, dominant indel size, and confidence level.

Figure 5. CasPINS Batch Analysis Mode. (A) Batch Analysis interface showing gene selection with detected samples, processing options, and progress tracking. (B) Summary results with pie chart array comparing editing efficiency across multiple clonal cell lines. Each pie displays edited in red versus wild-type in green proportions with percentage labels. Statistical summary provides average efficiency and sample-by-sample metrics.

---

## SUPPLEMENTARY DATA

Supplementary Data are available at NAR Online.
