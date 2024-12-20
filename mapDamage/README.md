## mapDamage

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/mapdamage2/README.html) [![Conda](https://img.shields.io/conda/dn/bioconda/mapdamage2.svg)](https://anaconda.org/bioconda/mapdamage2/files) 

![Conda](https://anaconda.org/bioconda/mapdamage2/badges/latest_release_date.svg) ![Conda](https://anaconda.org/bioconda/mapdamage2/badges/version.svg) [![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](https://www.repostatus.org/badges/latest/inactive.svg)](https://www.repostatus.org/#inactive)

#### `bioconda` installation

* python3 version **2.2.1**

``` 
conda install -c bioconda mapdamage2=2.2.1
```

* python3 version **2.2.1** **with** R and 4 mandatory packages for the Bayesian inference:

``` 
conda install -c bioconda mapdamage2=2.2.1=pyr36_1
```

---

### Important

* From version `2.2.1` the `master` branch is requiring **python3** as `python2` is not supported from 2020-01-01.

* Users with versions dating prior to June 12 2013 please update. A nasty bug that caused the statistical part of `mapDamage` to use half of the data for estimation of the damage parameters, sorry for the inconvenience.

### Introduction

Complete documentation, instructions, examples, screenshots and FAQ are available at [this address](http://ginolhac.github.io/mapDamage/).

[mapDamage2](https://geogenetics.ku.dk/publications/mapdamage2.0/) is a computational framework written in **Python3** and **R**, which tracks and quantifies DNA damage patterns
among ancient DNA sequencing reads generated by Next-Generation Sequencing platforms.

`mapDamage` was developed at the [Centre for GeoGenetics](https://geogenetics.ku.dk/) by the [Orlando Group ](https://geogenetics.ku.dk/research_groups/palaeomix_group/).

### Citation

If you use this program, please cite the following publication:  
Jónsson H, Ginolhac A, Schubert M, Johnson P, Orlando L.
[mapDamage2.0: fast approximate Bayesian estimates of ancient DNA damage parameters.
_Bioinformatics_ 23rd April 2013. doi: 10.1093/bioinformatics/btt193](http://bioinformatics.oxfordjournals.org/content/early/2013/05/17/bioinformatics.btt193)

The original `mapDamage1` was published in the following article:  
Ginolhac A, Rasmussen M, Gilbert MT, Willerslev E, Orlando L.
[mapDamage: testing for damage patterns in ancient DNA sequences. _Bioinformatics_ 2011 **27**(15):2153-5
http://bioinformatics.oxfordjournals.org/content/27/15/2153](http://bioinformatics.oxfordjournals.org/content/27/15/2153)

### Test the no-stats part and rescaling

you can test `mapDamage` by running:

``` 
cd mapDamage/mapdamage/
python3 mp_test.py
```

should return

``` 
Started with the command: /usr/local/bin/mapDamage -i tests/test.bam -r tests/fake1.fasta -d tests/results --no-stats
	Reading from 'tests/test.bam'
	Writing results to 'tests/results/'
pdf tests/results/Fragmisincorporation_plot.pdf generated
additional tests/results/Length_plot.pdf generated
Successful run
.
----------------------------------------------------------------------
Ran 2 tests in 3.357s

OK
```

### Contact

Please report bugs and suggest possible improvements to Aurélien Ginolhac, Mikkel Schubert or Hákon Jónsson by email:
ginolhac at gmail.com, mikkelsch at gmail.com or jonsson.hakon at gmail.com.
