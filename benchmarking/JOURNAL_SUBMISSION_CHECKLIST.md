# Journal Submission Checklist

## Journal Recommendation: Bioinformatics Advances (Primary) | BIOMAP (Backup)

---

## Part 1: Journal Comparison Decision

### Bioinformatics Advances (RECOMMENDED)

| Criterion | Assessment |
|-----------|-----------|
| **Scope match** | EXCELLENT - "bioinformatics methods including algorithms, statistics, databases, and software" |
| **Paper type** | Application Notes or Original Article - both fit CasPINS |
| **Publisher** | Oxford University Press + ISCB (International Society for Computational Biology) |
| **Review speed** | 27 days to first decision, 67 days to final |
| **Indexing** | PubMed Central, Scopus, Web of Science, Google Scholar |
| **APC** | ~$1,900 (20% discount for ISCB members) |
| **Open access** | Fully open access |
| **Audience** | Computational biologists and bioinformaticians - your primary users |
| **Submission** | mc.manuscriptcentral.com/bioadv |
| **Format** | Format-free at first submission |

**Why this is the best choice:**
1. NAR editor explicitly said the contribution is "bioinformatics-side" - this journal IS bioinformatics
2. Application Notes category is designed exactly for software tools like CasPINS
3. ISCB backing gives credibility in the computational biology community
4. Fast review with quality feedback
5. Format-free first submission = less reformatting work

### Biology Methods & Protocols (BIOMAP) - Backup

| Criterion | Assessment |
|-----------|-----------|
| **Scope match** | GOOD - accepts "Computational Methods" but broader biology focus |
| **Paper type** | Methods paper |
| **Publisher** | Oxford University Press |
| **Review speed** | ~75% acceptance rate after review |
| **Indexing** | PubMed Central, Scopus, Web of Science, Google Scholar |
| **APC** | ~$1,800 (estimated) |
| **Transfer path** | NAR editor offered direct manuscript transfer |
| **Audience** | Experimental biologists broadly |

**Advantages:** Direct transfer from NAR (easiest path), high acceptance rate
**Disadvantages:** Less prestige for bioinformatics tools, smaller computational biology audience

### Decision: Submit to **Bioinformatics Advances**

Reasons:
- Perfect scope alignment
- Better visibility in the computational biology community
- Application Notes format matches CasPINS
- The benchmarking additions we're making will strengthen the paper significantly
- BIOMAP remains a fallback if needed

---

## Part 2: Pre-Submission Checklist

### Manuscript Content Revisions

- [ ] **Add Feature Comparison Table (Table 1)**
  - Source: `benchmarking/feature_comparison_table.md`
  - Include as main-text table or supplementary table
  - Verify all feature claims are current (see MANUAL_BENCHMARKING_GUIDE.md Task C)

- [ ] **Add Benchmarking Results Section**
  - gRNA design concordance: CasPINS recovers 68.8% of CHOPCHOP, 67.2% of CRISPOR gRNAs (DONE - data ready)
  - Indel analysis: CasPINS analyzed 8/8 DDC samples; ICE 0/8, TIDE 0/8 due to low quality (DONE - data ready)
  - Synthetic validation results (from `benchmark_indel_analysis.py`) (DONE - data ready)

- [ ] **Add Worked Example**
  - TP53 workflow demonstration (from `worked_example_tp53.py` output)
  - Include as Figure or Supplementary Figure showing complete workflow

- [ ] **Strengthen Introduction**
  - Add paragraph on workflow fragmentation problem with quantified impacts
  - Reference: `benchmarking/integration_novelty_narrative.md`

- [ ] **Strengthen Discussion**
  - Replace mild "workflow continuity" language with quantified claims
  - Add: error reduction, time savings (80-90%), reproducibility
  - Add: explicit statement as first open-source integrated platform
  - Reference: `benchmarking/integration_novelty_narrative.md`

- [ ] **Update "Comparison with Existing Tools" section**
  - Reference the feature comparison table
  - Include benchmarking concordance data
  - Explicitly state what CasPINS provides that no other single tool offers

### Manuscript Formatting for Bioinformatics Advances

- [ ] **Remove any NAR-specific formatting**
  - Check DOCX for "Supplementary Data are available at NAR Online" (remove if present)
  - Update any NAR-specific header/footer formatting

- [ ] **Verify DOI**
  - Zenodo DOI already present: 10.5281/zenodo.18370068
  - Ensure it resolves correctly

- [ ] **Format for Bioinformatics Advances**
  - Check author guidelines: https://academic.oup.com/bioinformaticsadvances/pages/instructions-for-authors
  - Format-free at first submission (less reformatting needed!)
  - Decide on article type: Application Note (~2-4 pages) vs Original Article (longer)
  - **Recommendation**: Submit as **Original Article** since you have extensive content + benchmarking

- [ ] **Update cover letter**
  - Address to Bioinformatics Advances editors
  - Highlight: integrated workflow, benchmarking vs existing tools, open-source
  - Mention the NAR feedback constructively (shows you've addressed concerns)

### Code/Repository Requirements

- [ ] **Ensure GitHub repo is public and clean**
  - Repository: https://github.com/InnovationLine/CasPINS
  - Verify README is comprehensive
  - Verify installation instructions work

- [ ] **Verify Zenodo archive**
  - DOI: 10.5281/zenodo.18370068
  - Ensure it points to correct version

- [ ] **Add benchmarking materials to repository**
  - Include `benchmarking/` directory
  - Include feature comparison table
  - Include benchmark scripts

- [ ] **Test installation from scratch**
  - Clone repo → pip install → verify basic functionality

### Supplementary Materials

- [ ] **Supplementary Table S1**: Full feature comparison (if not in main text)
- [ ] **Supplementary Table S2**: gRNA benchmarking results (detailed)
- [ ] **Supplementary Table S3**: Indel analysis benchmarking results
- [ ] **Supplementary Figure S1**: Worked example workflow (if not in main text)
- [ ] **Supplementary Video** (optional): Screen recording of CasPINS workflow vs traditional approach

---

## Part 3: Submission Steps

1. Go to: https://mc.manuscriptcentral.com/bioadv
2. Create account if needed
3. Start new submission
4. Select article type: **Original Article**
5. Upload manuscript (DOCX)
6. Upload figures (high-resolution, individual files)
7. Upload supplementary materials
8. Enter metadata (title, authors, abstract, keywords)
9. Write cover letter
10. Submit

### Suggested Keywords
- CRISPR
- genome editing
- gRNA design
- indel analysis
- primer design
- bioinformatics workflow
- open-source software
- Sanger trace decomposition

### Suggested Cover Letter Structure

```
Dear Editors,

We submit "CasPINS: An Integrated Web-Based Platform for CRISPR/TALEN gRNA Design, 
Primer Generation, and Indel Decomposition Analysis" for consideration as an 
Original Article in Bioinformatics Advances.

CasPINS addresses a critical gap in genome editing workflows: no existing 
open-source tool integrates gRNA design, CRISPR-aware primer engineering, and 
Sanger-based indel quantification into a single platform. This integration 
eliminates error-prone manual data transfers between 3-5 separate tools, 
reduces computational workflow time by approximately 80%, and improves 
experimental reproducibility.

Key contributions:
- First open-source platform integrating the complete genome editing computational workflow
- Support for 90+ species and 14 CRISPR-Cas variants
- Benchmarking demonstrates concordance with CRISPOR (67.2%) and CHOPCHOP (68.8%) for gRNA design
- CasPINS indel analysis succeeds on traces that both TIDE and ICE reject (robustness advantage)
- Feature comparison shows unique capabilities not available in any single existing tool

The manuscript includes comprehensive benchmarking against existing tools, 
a detailed worked example, and a feature comparison table addressing the 
practical utility of the integrated approach.

We believe this work is well-suited for Bioinformatics Advances given its 
focus on bioinformatics software methodology and practical utility for the 
genome editing research community.

Sincerely,
[Authors]
```

---

## Part 4: Timeline

| Week | Task |
|------|------|
| Week 1 | ~~Run automated benchmarks, collect CRISPOR/CHOPCHOP/TIDE data manually~~ **DONE** |
| Week 2 | Revise manuscript (add tables, benchmarking, strengthen narrative) |
| Week 3 | Format, prepare supplementary materials, final review |
| Week 4 | Submit to Bioinformatics Advances |

---

## Part 5: If Rejected from Bioinformatics Advances

Fallback options in order of preference:
1. **BIOMAP** (Biology Methods & Protocols) - NAR offered transfer
2. **JOSS** (Journal of Open Source Software) - shorter review, software-focused
3. **BMC Bioinformatics** - established journal, good for software tools
4. **PeerJ Computer Science** - fast review, open access
5. **F1000Research** - open peer review, fast publication
