# FDEP

This package contains methods for investigating dependence in functional data. See the unpublished work

Hörmann, S., & Strenger, D. (2024), "Measuring Dependence between a Scalar  Response and a Functional Covariate".

for further details.

## Installation

The R-package `FDEP` can be installed and loaded via the commands

```
devtools::install_github("danielstrenger/FDEP")
library(FDEP)
```

## Help

The function for the test for independence is `soFun.test`. Type

```
help(soFun.test)
```

for further details. Data of the real data application are contained in `vaccine_data`.

## License

The package `FDEP` is licensed under GNU General Public License v3.0.

## References

Mona Azadkia and Sourav Chatterjee. A simple measure of conditional dependence. The Annals of Statistics, 49(6):3070–3102, 2021.

Statistik Austria. Bevölkerung zu Jahresbeginn 2021, 2021. URL https://statcube.at. Accessed: 23.10.2021.

Bundesministerium für Soziales, Gesundheit, Pflege und Konsumentenschutz (BMSGPK). Covid-19 Schutzimpfungen – Impfungen in Gemeinden, 2021. URL https://www.data.gv.at. Accessed: 23.10.2021.