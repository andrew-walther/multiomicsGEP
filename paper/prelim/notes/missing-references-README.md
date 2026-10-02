# Missing references for `project3-ssbmf.qmd`

Draft entries are in `missing-references.bib`, which follows the house style of
`bios-dissertation/prelim/references.bib`. They have not been merged into that
file or into `paper/ssbmf-refs.bib`. Checked 2026-10-02.

Status codes:
- **VERIFIED**: every field matches the named source.
- **PARTIAL**: the listed fields could not be checked and were left out.
- **UNCERTAIN**: the citable record could not be confirmed.

## Already in the master bib (no new entry)

| Reference | Existing key | Note |
|---|---|---|
| MOFA (Argelaguet et al. 2018, Mol Syst Biol) | `Argelaguet2018` | Present, PubMed-verified per its own note |
| iCluster (Shen, Olshen, Ladanyi 2009, Bioinformatics) | `Shen2009` | Present |
| EBMF (Wang & Stephens 2021, JMLR) | `Wang2021` | Present; flashier's CRAN page cites this as its methods paper |
| ebnm (Willwerscheid et al. 2025, JSS) | `Willwerscheid2025` | Present. This is the ebnm paper, not flashier |
| PACA-AU | `Bailey2016` | Present |

## New entries

| # | Reference | Key | Source used | Status |
|---|---|---|---|---|
| 1 | flashier R package | `Willwerscheid2023` | CRAN package page (authors, v1.0.7, published 2023-10-17, DOI 10.32614/CRAN.package.flashier); Crossref record for that DOI; arXiv search (no flashier preprint); local `packageVersion("flashier")` = 1.0.7 | VERIFIED as `@misc`. No JSS article or arXiv/bioRxiv preprint for flashier was found, so the CRAN package is the citable record. CRAN lists Wei Wang as an author, so he is included |
| 2 | Cao et al. 2021, Cell (CPTAC PDAC) | `Cao2021` | Crossref, DOI 10.1016/j.cell.2021.08.023 | VERIFIED. Author list taken from Crossref; the printed byline also credits the CPTAC consortium, which Crossref omits |
| 3 | Puleo et al. 2018, Gastroenterology | `Puleo2018` | Crossref, DOI 10.1053/j.gastro.2018.08.033 | VERIFIED |
| 4 | Dijk et al. 2020, Sci Rep | `Dijk2020` | PubMed PMID 31941932 and Crossref, DOI 10.1038/s41598-019-56826-9 | VERIFIED. Cohort match: the abstract reports RNA-seq of 90 resected PDAC specimens, which agrees with the repo's Dijk cohort (`docs/PDAC_data_audit.qmd`: RNA-seq, 90 samples). This is the only 2020 Dijk/Bijlsma PDAC RNA-seq paper found. The repo does not record the data's GEO/ArrayExpress accession, so confirm with the data provider if a formal accession is needed |
| 5 | Troyanskaya et al. 2001, Bioinformatics | `Troyanskaya2001` | Crossref, DOI 10.1093/bioinformatics/17.6.520 | VERIFIED |
| 6 | Owen & Perry 2009, Ann Appl Stat | `OwenPerry2009` | Crossref/doi.org BibTeX (vol, issue, DOI); Project Euclid article page (pages 564–594) | VERIFIED |
| 7 | Schwarz 1978, Ann Stat | `Schwarz1978` | Crossref/doi.org BibTeX (vol, issue, DOI); Project Euclid article page (pages 461–464) | VERIFIED |
| 8a | Breiman, Friedman, Olshen, Stone 1984, CART | `Breiman1984` | OpenLibrary ISBN 0534980538 (Wadsworth International Group, Belmont, CA, 1984); Crossref for author order | VERIFIED. DOI left out on purpose: 10.1201/9781315139470 is the 2017 Routledge reprint |
| 8b | Hastie, Tibshirani, Friedman 2009, ESL 2nd ed. | `Hastie2009` | Crossref, DOI 10.1007/978-0-387-84858-7; OpenLibrary ISBN 9780387848570 | VERIFIED. The 1-SE rule is in Sec. 7.10; cite one or both of 8a/8b |
| 9 | Efron & Tibshirani 1993, An Introduction to the Bootstrap | `EfronTibshirani1993` | OpenLibrary ISBN 0412042312 (Chapman & Hall, New York, 1993) | VERIFIED. DOI left out on purpose: 10.1201/9780429246593 is the CRC Press re-issue, dated 1994 in Crossref |
| 10 | Bender, Augustin, Blettner 2005, Stat Med | `Bender2005` | Crossref, DOI 10.1002/sim.2059 | VERIFIED |
| 12 | Young et al. 2026 (DeSurv) | `Young2026` (replacement) | Crossref (title, author and PNAS-prefix searches), PubMed (title and author searches), DeSurv GitHub README | UNCERTAIN. No published PNAS record exists in any of these sources as of 2026-10-02, so volume, number, pages and DOI are unknown. The draft changes the entry to `@unpublished` with a neutral note so it no longer renders as "In *PNAS*." with nothing after it. Replace with a full `@article` when the paper is published. Whether to describe it as "submitted" is the authors' call |

## Merging

- New keys do not clash with existing keys in the master bib.
- The two two-author keys (`OwenPerry2009`, `EfronTibshirani1993`) follow the
  master bib's existing pattern (`BhattacharyaDunson2011`, `KalbfleischPrentice2002`).
- Add entries to `bios-dissertation/prelim/references.bib` first, then copy
  them verbatim to `paper/ssbmf-refs.bib`, to keep the two files consistent.
