.. Prefactor documentation master file, created by
   sphinx-quickstart on Tue Nov 27 11:30:20 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Prefactor: Preprocessing for Facet Calibration for LOFAR
========================================================

**prefactor** is a pipeline to correct for various instrumental and ionospheric effects in both **LOFAR HBA** and **LOFAR LBA** observations.

It includes:

- removal of clock offsets between core and remote stations (using clock-TEC separation)
- correction of the polarization alignment between XX and YY
- robust time-independent bandpass correction
- ionospheric RM corrections with `RMextract`_
- removal of the element beam
- advanced flagging and interpolation of bad data
- mitigation of broad-band RFI and bad stations
- direction-independent phase correction of the target, using a global sky model from `TGSS ADR`_  or the new Global Sky Model (GSM)
- detailed diagnostics

It will prepare your data so that you will be able to use any direction-dependent calibration software, like `factor`_ or `killMS`_.

.. note::

    This documentation refers to the user version, intended to be run
    by users. For the production version, intended for automated
    processing on CEP4/LTA sites, see the documentation on the production branch.

    **NOTE:** The use of **prefactor** for the processing of long baselines (international stations) is not yet fully tested and thus **experimental**. Do not expect decent results if you use **prefactor** including long baselines.

Introduction
------------

.. toctree::
   :maxdepth: 2

   acknowledgements


Obtaining Prefactor
-------------------

.. toctree::
   :maxdepth: 2

   installation
   changelog


Setting Up and Running Prefactor
--------------------------------

.. toctree::
   :maxdepth: 2

   preparation
   parset
   running
   help


The Prefactor Pipelines
-----------------------

.. toctree::
   :maxdepth: 2

   pipelineoverview
   calibrator
   target
   concatenate
   initsubtract


.. _TGSS ADR: https://http://tgssadr.strw.leidenuniv.nl/
.. _RMextract: https://github.com/lofar-astron/RMextract/
.. _factor: https://github.com/lofar-astron/factor/
.. _killMS: https://github.com/saopicc/killMS/
