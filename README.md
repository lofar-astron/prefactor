# prefactor
## The LOFAR pre-facet calibration pipeline.

Parsets for the genericpipeline that do the first calibration of LOFAR data. Originally in order 
to prepare said data for the Factor facet calibration (https://github.com/lofar-astron/factor), but 
also useful if you don't plan to run Factor.

It includes:
* applying Ionospheric RM corrections
* clock-TEC separation with transfer of clock from the calibrator to the target
* some flagging and averaging of amplitude solutions
* grouping of subbands by actual frequency
* speed and disk usage improvements by optimized usage of NDPPP
* (optional) wide-band cleaning in Initial-Subtract 
* diagnostic plots
* at least some documentation on the [wiki pages](https://github.com/lofar-astron/prefactor/wiki)

The documentation can be found on the GitHub wiki pages: https://github.com/lofar-astron/prefactor/wiki

There are several pipeline parsets in this repository:
* Pre-Facet-Calibrator.parset : The calibrator part of the "standard" pre-facet calibration pipeline. 
* Pre-Facet-Target.parset : The target part of the "standard" pre-facet calibration pipeline. 
* Pre-Facet-Cal.parset : One parset calling first the calibrator and then the target pipelines. This is deprecated, please have a look at the [pipeline description](https://github.com/lofar-astron/prefactor/wiki/Documentation%3A-Pipelines#pre-facet-cal)
* Initial-Subtract.parset : A pipeline that generates full FoV images and subtracts the sky-models from the visibilities. (Needed for facet-calibration.)
* Initial-Subtract-Deep.parset : Same as Initial-Subtract.parset, but it does only one image of the full bandwidth instead of imaging the bands separately.

Experimental and thus deprecated for "normal" use are:
* Pre-Facet-Calibrator-RawSingle.parset : A version of a pre-facet calibrator pipeline to work on raw (non NDPPP'ed) data
* Pre-Facet-Calibrator-RawCombine.parset : A version of a pre-facet calibrator pipeline to work on raw (non NDPPP'ed) data that does the subband concatenating in the first NDPPP step. (To reduce the number of files on systems where this is a problem, e.g. JURECA)
* Pre-Facet-Target-RawSingle.parset : A version of a pre-facet target pipeline to work on raw (non NDPPP'ed) data
* Pre-Facet-Target-RawCombine.parset : A version of a pre-facet target pipeline to work on raw (non NDPPP'ed) data that does the subband concatenating in the first NDPPP step.
* Simple-Selfcal.parset : As the name says, an experimental selfcal pipeline.

Software requirements:
* the full "offline" LOFAR software installation version >= 2.17
* LoSoTo (version >=0.3 -- see https://github.com/revoltek/losoto)
* LSMTool (see https://github.com/darafferty/LSMTool)
* Python-PP (see http://www.parallelpython.com/ or https://pypi.python.org/pypi/pp )
* RMextract (see https://github.com/maaijke/RMextract)
* Python matplotlib
* WSClean 
  * for Initial-Subtract.parset : version >=1.12
  * for Initial-Subtract-Deep.parset : version >=2.0 (not yet released)
  * see https://sourceforge.net/projects/wsclean/
* APLpy (for Initial-Subtract)

The Pre-Facet-Calibration pipeline and its scripts where developed by:
* Martin Hardcastle <mjh somewhere extragalactic.info>
* George Heald <heald somewhere astron.nl>
* Andreas Horneffer <ahorneffer somewhere mpifr-bonn.mpg.de>
* Soumyajit Mandal <mandal somewhere strw.leidenuniv.nl>
* David Rafferty <drafferty somewhere hs.uni-hamburg.de>
* Carole Roskowinski <carosko gmail.com>
* Jose Sabater Montes <jsm somewhere iaa.es>
* Timothy Shimwell <shimwell somewhere strw.leidenuniv.nl>
* Sarrvesh Sridhar <sarrvesh somewhere astro.rug.nl>
* Reinout van Weeren <rvanweeren somewhere cfa.harvard.edu>
* Wendy Williams <wwilliams somewhere strw.leidenuniv.nl>

With special thanks to Stefan Froehlich <s.froehlich somewhere fz-juelich.de> for developing the 
genericpipeline.

The procedure is also mostly described in these papers:
* van Weeren et al. ApJ submitted
* Williams et al. MNRAS submitted


