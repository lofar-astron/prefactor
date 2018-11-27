.. _pipeline_overview:

The pipelines of prefactor
==========================

Prefactor contains pipeline parsets for the tasks needed to prepare the data for facet calibration:

``Pre-Facet-Calibrator``
    Processes the (amplitude-)calibrator to derive direction-independent corrections.
``Pre-Facet-Target``
    Transfers the direction-independent corrections to the target and does direction-independent calibration of the target.
``Image``
    Images the full FoV.
``Initial-Subtract``
    Images the full FoV (and 1st side-lobe), generating a sky-model and subtracting it from the visibilities.
