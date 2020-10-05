# mrgsolve_replace_NMsim
Purpose: use mrgsolve to replace `$SIM` and `vpc` with NONMEM/psn

The main motivation is to avoid using NONMEM for any simulation tasks (and associated dataset generation & CSV handling), so that the only NONMEM run you would do is $EST.  
You still need to develop `mrgsolve` model file, but from my experience that effort will pays off quite immediately.  


## Quick look

See https://genentech.github.io/mrgsolve_replace_NMsim/01_mrgsolve_replace_NMsim.html for output from the RMarkdown template


## How to use

You can download this repository from the green `Clone or download` button on top right.  
Once cloned or unarchived the ZIP file, you can double-click "mrgsolve_replace_NMsim.Rproj" to open this repository in RStudio. 

The main RMarkdown template file is `script/01_mrgsolve_replace_NMsim.Rmd`, which should be directly executable once you install all the necessary R packages from CRAN.


## Author

This tool was initiated by Kenta Yoshida (yoshida.kenta@gene.com) within Genentech Clinical Pharmacology Department.


## License 

The source code is distributed under a MIT license. Full description can be found in [LICENSE.txt](LICENSE.txt)

