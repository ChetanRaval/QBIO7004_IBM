# QBIO7004 - Computational Methods in Biology

**Chetan Raval | S4768622 | c.raval@uqconnect.edu.au**


## Aims
This repository contains a simplified stochastic event-based model that utilises tau-leaping to simulate the dynamics of protein folding within a generic biological system. This model attempts to better understand the implications of misfolding in neurodegenerative diseases.

## Folder Structure

- `code` contains `protein_model.R` which houses the entire simulation script. `lhs.R` and `lhs_bash.sh` were used solely for simulations running on the HPC
- `docs` contains the `index.html` file for GitHub Pages comptability
- `figs` contains the model diagram and plotting outputs from the script
- `references` contains `sessionInfo.txt` and two BiBTeX files for easier referencing in Zotero
- `report` contains the .rmd file and .html outputs of the markdown report

## Packages Required

    Bookdown (https://bookdown.org/home/)
    Tidyverse (https://www.tidyverse.org/)
    adaptivetau (https://cran.r-project.org/web/packages/adaptivetau/index.html)
    lhs (https://cran.r-project.org/web/packages/lhs/index.html)
    
Please see `references/sessionInfo.txt` for version information about R, the OS and attached or loaded packages at the time of writing the script.