# Data and software provenance

## Scientific inputs

| Analysis | Public source | Use in this repository |
| --- | --- | --- |
| §5.2 | [Tsunemoto et al. (2018)](https://doi.org/10.1038/s41586-018-0103-5) | Screen calls, corrected yield, physical pair slots and clone-identity evidence. The repository includes the required derived tables and source locators. |
| §5.2 | [GSE53951](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53951), [GSE114809](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114809), [GSE43916](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43916), [GSE101397](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101397) | Expression contrasts and chromatin accessibility used for the reported pair rankings and sensitivity analyses. |
| §5.4 | [Hammelman et al. (2022)](https://doi.org/10.1038/s41592-022-01522-2), [ReprogrammingRecovery](https://github.com/gifford-lab/ReprogrammingRecovery) | Candidate universe and historical references. The three endpoints are P0 liver, P1 PCM1-positive cardiomyocyte nuclei, and the pancreatic β-TC6 cell line. |
| §5.4 | [GSE66534](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66534), [GSE72744](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72744), [GSM2525591](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2525591), [GSM2666437](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2666437), [ENCSR096STK](https://www.encodeproject.org/experiments/ENCSR096STK/), [ENCSR609OHJ](https://www.encodeproject.org/experiments/ENCSR609OHJ/) | RNA and ATAC inputs. Consult the bundled sample and retrieval tables for exact accessions, roles, and files. |

The subsystem manifests record accession-level locations, resource versions and hashes where available. The manifest `INPUT_SHA256.tsv` hashes the packaged copies. If a portable adaptation changed bytes, its packaged checksum and historical source checksum serve different purposes.

## Upstream code and terms

The upstream terms below were checked against the local source snapshots. Obtain upstream packages from their release locations. A citation or public download URL does not carry the license of this repository.

| Upstream component | Source or frozen release | License evidence |
| --- | --- | --- |
| ReprogrammingRecovery | [Repository](https://github.com/gifford-lab/ReprogrammingRecovery) | CC0 1.0 in `LICENSE.txt`. |
| TRIAGE | [e57492c6eb67f4ac96a7ccb121496f12d3225414](https://github.com/woojunshim/TRIAGE/tree/e57492c6eb67f4ac96a7ccb121496f12d3225414) | GNU GPL v3 in upstream `LICENSE`. |
| VIPER | [Bioconductor viper](https://bioconductor.org/packages/viper/) | Columbia software evaluation license in the inspected 1.40.0 source. The package and its code are obtained separately. |
| Taiji | [Taiji v1.3.0](https://github.com/Taiji-pipeline/Taiji/releases/tag/v1.3.0) | BSD-style three-clause license in the Taiji source snapshot. |
| Taiji-reprogram | [0b6337aae2bebdcb5f1db35031ecc38e8b49179a](https://github.com/Wang-lab-UCSD/Taiji-reprogram/tree/0b6337aae2bebdcb5f1db35031ecc38e8b49179a) | MIT in the local upstream snapshot. |
| ANANSE | [580907581bda42e7bae172bcd346076cb87dcd41](https://github.com/vanheeringen-lab/ANANSE/tree/580907581bda42e7bae172bcd346076cb87dcd41) | MIT in upstream `LICENSE`. |
| IMCell | [v0.1.0 archive](https://doi.org/10.5281/zenodo.20563020) | MIT in upstream `LICENSE` and package metadata. |
| CompForce | [v0.1.0 archive](https://doi.org/10.5281/zenodo.20563029) | MIT in upstream `LICENSE` and package metadata. |
| SFA | [Repository](https://github.com/dwgoon/sfa) | BSD 3-Clause in upstream `LICENSE`. |
| TFcomb | [d80164bef403d567206c0127d802f9dce90aaba8](https://github.com/Chen-Li-17/TFcomb/tree/d80164bef403d567206c0127d802f9dce90aaba8) | MIT in upstream `LICENSE`. The publication discloses patent application 202410023244.5. |
| CollecTRI | [Zenodo resource](https://doi.org/10.5281/zenodo.8192729) | Obtain the frozen network and any required mapping resources using the corresponding source manifest. |
| TRRUST v2 | [Database license](https://www.grnpedia.org/trrust//?homeType=newUser), [Han et al.](https://doi.org/10.1093/nar/gkx1013), [frozen source table](data/pairs/provenance/SOURCE_TRRUST.tsv) | The network release is CC BY-SA 4.0. The article text has separate CC BY-NC terms. See the pair workflow for attribution and the derived files that use the network. |

The archive contains manuscript-analysis scripts and lightweight derived data. The repository owners should select the license for their original code and data, complete citation metadata and supply their release DOI at publication.
