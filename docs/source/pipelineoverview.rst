.. _pipeline_overview:

Pipeline overview
=================

**Prefactor** is organized in three major parts to process **LOFAR** data:

    .. image:: prefactor_workflow_sketch.png

``Pre-Facet-Calibrator``
    Processes the (amplitude-)calibrator to derive direction-independent corrections. See :ref:`calibrator_pipeline` for details.
``Pre-Facet-Target``
    Transfers the direction-independent corrections to the target and does direction-independent calibration of the target. See :ref:`target_pipeline` for details.
``Concatenate``
    Concatenates the single-subband target data retrieved from the LTA into bands suitable for further processing with the initial-subtract pipeline. See :ref:`concatenate_pipeline` for details.
``Initial-Subtract``
    Images the full FoV (and 1st side-lobe), generating a sky-model and subtracting it from the visibilities. See :ref:`initsubtract_pipeline` for details.
``Pre-Facet-Image``
    Images the full FoV using the full bandwidth. See :ref:`image_pipeline` for details.

