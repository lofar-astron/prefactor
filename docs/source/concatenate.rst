.. _concatenate_pipeline:

Concatenate pipeline
=========================

This pipeline concatenates single-subband target data produced by production
runs and retrieved through the LTA. The resulting concatenated files are
required for further processing with the initial-subtract pipeline.

.. note::

    If you processed the target data yourself instead of retrieving them from the
    LTA, this pipeline is not required as the concatenated datasets are already
    produced by the user version of the target pipeline.


Prepare data
------------
This part of the pipeline prepares the target data in order to be concatenated. The steps are
as follows:

``createmap_target``
    Generate a mapfile of all the target data (the single-subband datasets retrieved
    from the LTA, with the direction-independent phase-only calibration applied).
``combine_target_map``
    Generate a mapfile with all files in a single entry. This mapfile is used as
    input to the next step.
``sortmap_target``
    Compute frequency groupings
``do_magic_maps``, ``do_sortmap_maps``
    Convert the output of do_magic into usable mapfiles.


Concatenation
-------------
Subbands are concatenated into "bands".

``dpppconcat``
    Concatenate the data, averaging to the specified frequency and time resolution.
``make_results_mapfile``, ``move_results``
    Move the concatenated files to the results directory.



User-defined parameter configuration
------------------------------------

**Parameters you will need to adjust**

*Information about the input data*

``! target_input_path``
    Directory where your single-subband target data are stored.
``! target_input_pattern``
    Regular expression pattern of all your target files.
    .. note::

        These files should have the direction-independent calibration applied to the DATA
        column.

*Location of the software*

``! prefactor_directory``
    Path to your prefactor copy

**Parameters you may need to adjust**

*Interpolation options*

- ``interp_windowsize``: size of the window over which a value is interpolated. Should be odd. (default: 15)

*Averaging options*

- ``avg_timeresolution_concat``: final time resolution of the data in seconds after averaging and concatenation (default: 8)
- ``avg_freqresolution_concat``: final frequency resolution of the data after averaging and concatenation (default: 97.64kHz, which translates to 2 channels per subband)

*Concatenation options*

- ``num_SBs_per_group``: make concatenated measurement-sets with that many subbands (default: 10)
- ``reference_stationSB``: station-subband number to use as reference for grouping, (default: ``None`` -> use lowest frequency input data as reference)
