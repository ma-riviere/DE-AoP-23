<div align="center">
 
 <h1>Apnea of Prematurity and transcriptional cerebellar changes</h1>
 <h3><i>Main branch</i></h3>

 [![License](https://img.shields.io/badge/license-CCBY-blue.svg)](/LICENSE)
 [![DOI](https://zenodo.org/badge/665567093.svg)](https://zenodo.org/badge/latestdoi/665567093)
 
 <a href = "https://ma-riviere.github.io/DE-AoP-23/">Link to the documentation</a>

</div>

> **Note**  
> This repository contains the data and R code for the *"Apnea of Prematurity induces short and long-term development-related transcriptional changes in the murine cerebellum"* paper

## ❔ Requirements:

- R version 4.3 or newer
- R Studio version 2022.07 or newer


## 💻 Repository structure:

- `DE-AoP-23.RProj`: R Studio project (open this first).
- `analysis`: R Markdown files for the data analysis, split by data type (PCR & IHC). 
  - The first code chunk of any of the .Rmd files will install and load all the packages required for the project, based on the `renv.lock` file.
- `data`: The PCR and IHC data, both raw and processed (when applicable).
- `src`: R scripts declaring the functions called within the analysis files (e.g. `viz.R` for the figures, `data.R` for the data loading).
- `_configs.yml`: Lists the paths to various external files used within the code (e.g. data).
- `_dependencies.yml`: Lists the packages required for this project (which will be auto-installed based on the existing `renv.lock` file).


## 📜 Licence:

[CC-BY](LICENSE)


## 💬 Citation:

- **Paper:** Rodriguez-Duboc, A., Basille-Dugay, M., Debonne, A., Rivière, M.-A., Vaudry, D., & Burel, D. (2023). Apnea of prematurity induces short and long-term development-related transcriptional changes in the murine cerebellum. *Current Research in Neurobiology, 5*, 100113. https://doi.org/10.1016/j.crneur.2023.100113

- **Code:** Marc-Aurèle Rivière, & Agalic Rodriguez-Duboc. (2023). ma-riviere/DE-AoP-23: public release (v1.0). Zenodo. https://doi.org/10.5281/zenodo.8139284


## ✨ Contributors:

- **Marc-Aurèle Rivière**:  
[![ORCID](https://img.shields.io/badge/ORCID-A6CE39?style=flat-square&labelColor=white&logo=orcid&logoColor=A6CE39)][ORCID_MAR]
[![Research Gate](https://img.shields.io/badge/ResearchGate-00CCBB?style=flat-square&labelColor=white&logo=researchgate&logoColor=00CCBB)][RG_MAR]

- **Agalic Rodriguez-Duboc**:  
[![ORCID](https://img.shields.io/badge/ORCID-A6CE39?style=flat-square&labelColor=white&logo=orcid&logoColor=A6CE39)][ORCID_ARD]
[![Research Gate](https://img.shields.io/badge/ResearchGate-00CCBB?style=flat-square&labelColor=white&logo=researchgate&logoColor=00CCBB)][RG_ARD]


## 📫 Contact:

For any questions, please contact the primary author of the paper, **Agalic Rodriguez-Duboc**:  
<a href="mailto:agalic.rd@gmail.com?subject=Apnea%20of%20Prematurity%20and%20transcriptional%20cerebellar%20changes">![Gmail](https://img.shields.io/badge/Gmail-C71610?style=flat-square&labelColor=white&logo=Gmail&logoColor=C71610)</a>


<!----------------------------------->

[RG_MAR]: https://www.researchgate.net/profile/Marc_Aurele_Riviere2
[ORCID_MAR]: https://orcid.org/0000-0002-5108-3382
[RG_ARD]: https://www.researchgate.net/profile/Agalic-Rodriguez-Duboc
[ORCID_ARD]: https://orcid.org/0000-0002-2084-3780
[RG_ARD]: https://www.researchgate.net/profile/Agalic-Rodriguez-Duboc
[ORCID_ARD]: https://orcid.org/0000-0002-2084-3780
