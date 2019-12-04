# dplfit

The self-similar discrete power law distribution — R distribution, density,
quantile and random deviate functions, and the maximum likelihood estimator.
This repository also includes OCaml code for generating self-similar random
trees.

dplfit was written by Mitchell Newberry <mitchell@localpost.io> and is (c)
Mitchell Newberry 2019.  Bug reports, comments, and feature requests are
welcome.

## Data Availability and Reproducibility

Data to reproduce the results of Newberry and Savage (2019) are in `data/`.

`figures.R` reproduces figures and numerical results from these files. 

`yeh1976tracheobronchial.D6.tsv` was hand-transcribed from a US Government
report as cited in Newberry and Savage (2019).

`mouse-vasculature-tekin2016.tsv` is derived from Angicart scans of MicroCT
data from Tekin et al. (2016) as used in Newberry and Savage (2019).

`earthcat/mags.all` is assembled by `make_earthcat.sh` from mirrors of the
Southern California Seismographic Network data catalog on internet archive to
reproduce the data source, time interval and results of Bak et al. (2002).

## References

Bak, P., Christensen, K., Danon, L., & Scanlon, T. (2002). Unified scaling law for earthquakes. Physical Review Letters, 88(17), 178501.

Newberry, M. G., & Savage, V. M. (2019). Self-similar processes follow a power law in discrete logarithmic space. Physical review letters, 122(15), 158303.

Tekin, E., Hunt, D., Newberry, M. G., & Savage, V. M. (2016). Do vascular networks branch optimally or randomly across spatial scales? PLoS computational biology, 12(11), e1005223.
