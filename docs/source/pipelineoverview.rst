.. _pipeline_overview:

Pipeline overview
=================

Prefactor contains pipeline parsets for the tasks needed to prepare the data for facet calibration:

``Pre-Facet-Calibrator``
    Processes the (amplitude-)calibrator to derive direction-independent corrections. See :ref:`calibrator_pipeline` for details.
``Pre-Facet-Target``
    Transfers the direction-independent corrections to the target and does direction-independent calibration of the target. See :ref:`target_pipeline` for details.
``Initial-Subtract``
    Images the full FoV (and 1st side-lobe), generating a sky-model and subtracting it from the visibilities. See :ref:`initsubtract_pipeline` for details.
