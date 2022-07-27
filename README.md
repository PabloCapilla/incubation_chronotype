

# Reproductive fitness is associated with female chronotype in a songbird

This repository contains materials used for a project investigating incubation in female great tits, incubation chronotype and the association between chronotype and reproductive success.

---

Robyn J. Womack^#, Pablo Capilla-Lasheras^#, Ciara L. O. McGlade,  Davide M. Dominoni, Barbara Helm. **Reproductive fitness is associated with female chronotype in a songbird**. *bioRxiv*. DOI: 10.1101/2022.07.01.498449v1
^# joint first authors

---

For any further information, please contact: [Pablo Capilla-Lasheras](https://scholar.google.com/citations?hl=en&user=5JMTO-kAAAAJ&view_op=list_works&sortby=pubdate), email: pacapilla@gmail.com

## Code:

All R code is available in [`scripts`](https://github.com/PabloCapilla/incubation_chronotype/tree/main/scripts). Scripts are divided into [`data prep scripts`](https://github.com/PabloCapilla/incubation_chronotype/tree/main/scripts/01_data_prep) and [`statistical models`](https://github.com/PabloCapilla/incubation_chronotype/tree/main/scripts/02_models). within each sub-folder, scripts are numbered. The paths provided to import datasets into R assume your location is the main folder of the repository (i.e., the general folder `incubation_chronotype`). The R version used for this project was 4.2.0.

## Folders:

[`data`](https://github.com/PabloCapilla/incubation_chronotype/tree/main/data): contains two files, one for [`raw data`] (https://github.com/PabloCapilla/incubation_chronotype/tree/main/data/01_raw_data) and one with a dataset ready to be analysed ([`data_incubation.RDS`](https://github.com/PabloCapilla/incubation_chronotype/tree/main/data)).

[`scripts`](https://github.com/PabloCapilla/incubation_chronotype/tree/main/scripts): R code to produce results and figures included in the manuscript. It contains three subfolders:
* [`01_data_prep`](https://github.com/PabloCapilla/incubation_chronotype/tree/main/scripts/01_data_prep): script to clean and generate a table for further analysis (saved as [`data_incubation.RDS`](https://github.com/PabloCapilla/incubation_chronotype/tree/main/data) in [`data`](https://github.com/PabloCapilla/incubation_chronotype/tree/main/data)).
* [`02_models`](https://github.com/PabloCapilla/incubation_chronotype/tree/main/scripts/02_models): scripts for statistical models and visualisation of results.

[`plots`](https://github.com/PabloCapilla/incubation_chronotype/tree/main/plots): contains all figures (main and supplementary ones) created for this study, all of which were created using the R package [ggplot2 v.3.3.6](https://cran.r-project.org/web/packages/ggplot2/index.html).

## Notes

See details of the licence of this repository in [`LICENSE`](https://github.com/PabloCapilla/incubation_chronotype/blob/main/LICENSE).
