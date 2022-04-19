# Prefactor

## IMPORTANT NOTICE
This repository contains the old version of prefactor. The new version of the pipeline is called LOFAR Initial Calibration (LINC) pipeline and can be found [here](https://git.astron.nl/RD/LINC).

## The (old) LOFAR pre-facet calibration pipeline.

**prefactor** is a pipeline to correct for various instrumental and ionospheric effects in both **LOFAR HBA** and **LOFAR LBA** observations.
It will prepare your data so that you will be able to use any direction-dependent calibration software, like [Rapthor](https://github.com/darafferty/rapthor), [factor](https://github.com/lofar-astron/factor), or [killMS](https://github.com/saopicc/killMS/).

It includes:
* removal of clock offsets between core and remote stations (using clock-TEC separation)
* correction of the polarization alignment between XX and YY
* robust time-independent bandpass correction
* ionospheric RM corrections with [RMextract](https://github.com/lofar-astron/RMextract/)
* removal of the element beam
* advanced flagging and removal of bad data
* mitigation of broad-band RFI and bad stations
* direction-independent phase correction of the target, using a global sky model from [TGSS ADR](https://http://tgssadr.strw.leidenuniv.nl/)  or the new Global Sky Model [GSM](http://172.104.228.177/)
* detailled diagnostics

The full documentation can be found at the [prefactor webpage](https://www.astron.nl/citt/prefactor/).

**WARNING**: The current skymodel used for 3C295 is not using the Scaife&Heald flux density scale.

### Software requirements
* the full "offline" LOFAR software installation (version >= 3.1)
* [DPPP](https://github.com/lofar-astron/DP3) (version v5.1)
* [LoSoTo](https://github.com/revoltek/losoto) (version of Jun 15, 2020, commit [c8fbd61](https://github.com/revoltek/losoto/tree/c8fbd6194074bef4009cb66dd7ecd59e98664d63))
* [LSMTool](https://github.com/darafferty/LSMTool) (version>= 1.4.2)
* [EveryBeam](https://git.astron.nl/RD/EveryBeam) (commit 62ed313d or later)
* [RMextract](https://github.com/maaijke/RMextract) (commit a618fff or later)
* Python (including matplotlib, scipy, and astropy)
* [AOFlagger](https://gitlab.com/aroffringa/aoflagger) (version 2.14)
* [WSClean](https://gitlab.com/aroffringa/wsclean) (for Initial-Subtract; version >= 2.5)
* for Initial-Subtract-IDG(-LowMemory).parset and Pre-Facet-Image.parset: WSClean must be compiled with [IDG](https://gitlab.com/astron-idg/idg)
* APLpy (for Initial-Subtract and Pre-Facet-Image)

### Installation
The recommended way to install prefactor is to download it from github with:

```
git clone https://github.com/lofar-astron/prefactor.git
```

This allows for easy updating of the code to include bugfixes or new features.
It is also possible to download tar files of releases from the [release page](https://github.com/lofar-astron/prefactor/releases).

Once downloaded, the installation is complete; to set up a run, see the detailed setup information at the [prefactor webpage](https://www.astron.nl/citt/prefactor/).

### Directory structure
prefactor contains the following sub-directories:
* **bin**: scripts for your convenience
* **docs**: contains the documentation for prefactor
* **plugins**: scripts for manipulating mapfiles
* **rfistrategies**: strategies for statistical RFI mitigation using [AOFlagger](https://sourceforge.net/p/aoflagger/wiki/Home/)
* **scripts**: scripts that the pipeline calls to process data, generate plots, etc.
* **skymodels**: skymodels that are used by the pipeline (e.g. for demixing or calibrating the calibrator)
* **solutions**: template solutions for the use of non-supported calibrator sources


The main directory contains the different parsets for the genericpipeline:
* Pre-Facet-Calibrator.parset : The calibrator part of the "standard" pre-facet calibration pipeline.
* Pre-Facet-Target.parset : The target part of the "standard" pre-facet calibration pipeline.
* Concatenate.parset : A pipeline that concatenates single-subband target data to produce concatenated bands suitable for the initial-subtract pipeline.
* Initial-Subtract.parset : A pipeline that generates full-FoV images and subtracts the sky-models from the visibilities. (Needed for facet-calibration.)
* Initial-Subtract-IDG.parset : Same as Initial-Subtract-Fast.parset, but uses the image domain gridder (IDG) in WSClean.
* Initial-Subtract-IDG-LowMemory.parset : Same as Initial-Subtract-Fast.parset, but uses the image domain gridder (IDG) in WSClean for high-res imaging.
* Pre-Facet-Image.parset : A pipeline that generates a full-bandwidth, full-FoV image.
* make\_calibrator/target\_plots.losoto\_parset : Losoto parsets for making diagnostic plots from the output h5parms.

The Pre-Facet-Calibration pipeline and its scripts where developed by:
* Alexander Drabent
* David Rafferty
* Andreas Horneffer
* Francesco de Gasperin
* Marco Iacobelli
* Emanuela Orru
* Björn Adebahr
* Martin Hardcastle
* George Heald
* Soumyajit Mandal
* Carole Roskowinski
* Jose Sabater Montes
* Timothy Shimwell
* Sarrvesh Sridhar
* Reinout van Weeren
* Wendy Williams

With special thanks to Stefan Fröhlich for developing the genericpipeline.

### Acknowledgements
The Prefactor v3 procedure is described in this paper:
* de Gasperin, F.; Dijkema, T. J.; Drabent, A.; Mevius, M.; Rafferty, van Weeren, R., et al. 2019, [A&A, 662, A5](http://adsabs.harvard.edu/abs/2018arXiv181107954D)

The Factor procedure is described in these papers:
* van Weeren, R. J., Williams, W. L., Hardcastle, M. J., et al. 2016, [ApJS, 223, 2](http://adsabs.harvard.edu/abs/2016ApJS..223....2V)
* Williams, W. L., van Weeren, R. J., Röttgering, H. J. A., et al. 2016, [MNRAS, 460, 2385W](http://adsabs.harvard.edu/abs/2016MNRAS.460.2385W)


