# prefactor
## The LOFAR pre-facet calibration pipeline.

Parsets for the genericpipeline that do the first calibration of LOFAR data. Originally in order 
to prepare said data for the Factor facet calibration (https://github.com/lofar-astron/factor), but 
also useful if you don't plan to run Factor.

It includes:
* clock-TEC separation with transfer of clock from the calibrator to the target
* some flagging and averaging of amplitude solutions
* diagnostic plots
* at least some documentation

The Pre-Facet-Calibration pipeline and its scripts where developed by:
* Reinout van Weeren <rvanweeren somewhere cfa.harvard.edu>
* Wendy Williams <wwilliams somewhere strw.leidenuniv.nl>
* Martin Hardcastle <mjh somewhere extragalactic.info>
* Timothy Shimwell <shimwell somewhere strw.leidenuniv.nl>
* David Rafferty <drafferty somewhere hs.uni-hamburg.de>
* Jose Sabater Montes <jsm somewhere iaa.es>
* George Heald <heald somewhere astron.nl>
* Sarrvesh Sridhar <sarrvesh somewhere astro.rug.nl>
* Andreas Horneffer <ahorneffer somewhere mpifr-bonn.mpg.de>

With special thanks to Stefan Froehlich <s.froehlich somewhere fz-juelich.de> for developing the 
genericpipeline.

This procedure is also fully described in these papers:
* van Weeren et al. ApJ submitted
* Williams et al. MNRAS submitted


