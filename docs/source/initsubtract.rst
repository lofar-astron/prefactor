.. _initsubtract_pipeline:

Intial-subtract pipeline
========================

This pipeline images the full FoV (and 1st side-lobe) at two resolutions and at
multiple frequencies, generating a sky-model and subtracting it from the
visibilities. The parset is named one of ``Initial-Subtract.parset``,
``Initial-Subtract-IDG.parset``, or ``Initial-Subtract-IDG-LowMemory.parset``,
depending on whether one wants to use IDG with WSClean. IDG is generally much
faster if you have GPUs.

.. note::

    At this time, only HBA data are supported.


Options
-------

``! data_input_path``
    Directory where your concatenated target data are stored.
``! data_input_pattern``
    Regular expression pattern of all your calibrator files.
``! data_input_pattern``
    Full path to the direction-independent calibration solutions.
``! cellsize_highres_deg``
    Cellsize in degrees for high-resolution images.
``! cellsize_lowres_deg``
    Cellsize in degrees for low-resolution images.
``! fieldsize_highres``
    Size of the high-resolution image is this value times the FWHM of mean semi-major axis of
    the station beam.
``! fieldsize_lowres``
    Size of the low-resolution image is this value times the FWHM of mean semi-major axis of
    the station beam.
``! maxlambda_highres``
    Maximum uv-distance in lambda that will be used for the high-resolution imaging.
``! maxlambda_lowres``
    Maximum uv-distance in lambda that will be used for the low-resolution imaging.
``! image_padding``
    How much padding shall we add during the imaging?
``! nbands_image``
    Number of bands to image (spread over the full bandwidth). Larger values
    result in better subtraction but longer runtimes.
``! min_flux_jy``
    Minimum flux density in Jy of clean components from the high-resolution
    imaging to include in subtract_high step.
``! idg_mode``
    IDG mode to use: cpu or hybrid (= CPU + GPU).
``! local_scratch_dir``
    Scratch directory for wsclean (can be local to the processing nodes!).


Steps
-----

``create_ms_map``
    Generate a mapfile of all the target data. The files must be supplied as a
    list of the full paths to the files.
``combine_mapfile``
    Generate a mapfile with all files in a single entry. This mapfile is used as
    input to the next step.
``do_magic``
    Compute frequency groupings, image sizes, and averaging values using the MS
    files from the previous step.
``do_magic_maps``
    Convert the output of do_magic into usable mapfiles.
``create_h5parm_map``
    Create a mapfile with the direction independent h5parm.
``expand_h5parm_mapfile``
    Expand the h5parm mapfile so that there is one entry for every file.
``select_imaging_bands``
    Select bands spread over the full bandwidth for imaging.
``select_apply_high_files``
    Select files spread over the full bandwidth for imaging.
``select_apply_h5parm``
    Adjust the dir-indep h5parm mapfile to match the selected bands.
``gsmcal_apply``
    Apply direction-independent solutions.
``select_high_size``
    Adjust the high_size mapfile to match the selected bands.
``select_high_nwavelengths``
    Adjust the nwavelengths mapfile to match the selected bands.
``wsclean_high``
    Image the data with WSClean to make the high-resolution images. The images will
    automatically be stretched along the y-axis to account for the elongation of the
    primary beam as a function of average elevation. A typical image at
    lower Declination (+7 degrees) looks like the one below.

    .. image:: initsub_high_image.png

``wsclean_low``
    Image the data with WSClean to make the low-resolution images. The images will
    automatically be stretched along the y-axis to account for the elongation of the
    primary beam as a function of average elevation. A typical image at
    lower Declination (+7 degrees) looks like the one below.

    .. image:: initsub_low_image.png

