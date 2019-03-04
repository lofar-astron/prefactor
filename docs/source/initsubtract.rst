.. _initsubtract_pipeline:

Initial-subtract pipeline
=========================

This pipeline images the full FoV (and first side lobe) at two resolutions and at
multiple frequencies, generating a sky model and subtracting it from the
visibilities. This pipeline need only be run if you want to use Factor to do the
direction-dependent imaging. The parset is named one of ``Initial-Subtract.parset``,
``Initial-Subtract-IDG.parset``, or ``Initial-Subtract-IDG-LowMemory.parset``,
depending on whether one wants to use IDG with WSClean. IDG is generally much
faster than the normal WSClean if you have GPUs.

.. note::

    At this time, only HBA data are supported.


Prepare data
------------
This part of the pipeline prepares the target data in order to be imaged. The steps are
as follows:

``create_ms_map``
    Generate a mapfile of all the target data (the concatenated datasets output by the
    target pipeline, with the direction-independent phase-only calibration applied).
``combine_mapfile``
    Generate a mapfile with all files in a single entry. This mapfile is used as
    input to the next step.
``do_magic``
    Compute frequency groupings, image sizes, and averaging values using the MS
    files from the previous step. The image size is calculated from the FWHM of the
    primary beam at the lowest frequency at the mean elevation of the observation.
``do_magic_maps``
    Convert the output of do_magic into usable mapfiles.
``create_h5parm_map``
    Create a mapfile with the direction independent h5parm.
``expand_h5parm_mapfile``
    Expand the h5parm mapfile so that there is one entry for every file.
``select_imaging_bands``
    Select bands spread over the full bandwidth for imaging.
``select_high_size``
    Adjust the high_size mapfile to match the selected bands.
``select_high_nwavelengths``
    Adjust the nwavelengths mapfile to match the selected bands.


Imaging and subtraction
-----------------------
Imaging is done at two resolutions to fully cover the expected range of source structure.
WSClean is used to produce the images. See the parset and the do_magic step above
for details of the parameters used. They are chosen to produce good results for
most standard observations.

``wsclean_high``
    Image the data with WSClean to make the high-resolution images. The images will
    automatically be stretched along the y-axis to account for the elongation of the
    primary beam as a function of average elevation. A typical image at
    lower Declination (+7 degrees) looks like the one below.

    .. image:: initsub_high_image.png

``mask_high``
    Make masks for the high-res images. Masks are used to exclude artifacts from
    being included in the subtract steps.
``mk_inspect_dir``
    Create the inspection_directory if needed.
``copy_mask``
    Copy the mask images to where we want them.
``plot_im_high``
    Plot the high-res image and mask as png files. Such an image is show above.
``move_high``
    Move the high-res images to where we want them.
``create_maxsize_high_map``
    Make a mapfile with maximum image size.
``pad_model_high``
    Pad the model images to a uniform size.
``pad_mask_high``
    Pad the mask images to a uniform size.
``combine_model_high_mapfile``
    Compress the model_high mapfile.
``expand_model_high``
    Expand the model_high mapfile so that there is one entry for every band.
``combine_mask_high_mapfile``
    Compress the mask_high mapfile.
``expand_mask_high``
    Expand the mask high mapfile so that there is one entry for every band.
``fits_to_bbs_high``
    Convert high-res model images to sky models that are understood by DPPP.
``make_sourcedb_high``
    Make sourcedbs from the high-res sky models.
``expand_sourcedb_high``
    Expand the sourcedb mapfile so that there is one entry for every file.
``subtract_high``
    Predict, corrupt, and subtract the high-resolution model. The subtraction is
    done from the DATA column to the SUBTRACTED_DATA_HIGH column. The SUBTRACTED_DATA_HIGH
    column is imaged later in the ``wsclean_low`` step to pick up any emission missed in
    the high-resolution image.
``select_low_size``
    Adjust the low size mapfile to match the selected bands.
``select_low_nwavelengths``
    Adjust the low nwavelengths mapfile to match the selected bands.
``wsclean_low``
    Image the data (after subtraction of the high-resolution model) with WSClean
    to make the low-resolution images. The images will automatically be
    stretched along the y-axis to account for the elongation of the primary beam
    as a function of average elevation. A typical image at lower Declination (+7
    degrees) looks like the one below.

    .. image:: initsub_low_image.png

``mask_low``
    Make masks for the low-res images. Masks are used to exclude artifacts from
    being included in the subtract steps.
``plot_im_low``
    Plot the low-res image and mask as png files. Such an image is show above.
``move_low``
    Move the low-res images to where we want them.
``create_maxsize_low_map``
    Make a mapfile with maximum image size.
``pad_model_low``
    Pad the model images to a uniform size.
``pad_mask_low``
    Pad the mask images to a uniform size.
``combine_model_low_mapfile``
    Compress the model_low mapfile.
``expand_model_low``
    Expand the model_low mapfile so that there is one entry for every band.
``combine_mask_low_mapfile``
    Compress the mask_low mapfile.
``expand_mask_low``
    Expand the mask low mapfile so that there is one entry for every band.
``fits_to_bbs_low``
    Convert low-res model images to sky models.
``make_sourcedb_low``
    Make sourcedbs from the low-res sky models.
``expand_sourcedb_low``
    Expand the sourcedb mapfile so that there is one entry for every file.
``subtract_low``
    Predict, corrupt, and subtract the low-resolution model. The subtraction is
    done from the SUBTRACTED_DATA_HIGH column to the SUBTRACTED_DATA_ALL column.
    Therefore, the SUBTRACTED_DATA_ALL column contains the final residual data needed
    for Factor.
``merge``
    Merge the high-res and low-res sky models together. These sky models are used
    by Factor to add sources back before calibration.
``copy_skymodels``
    Copy the merged sky models to the directory with the input data.
``createmap_plots``
    Create a map with the generated plots.
``move_plots``
    Move the plots to the inpection directory.



User-defined parameter configuration
------------------------------------

**Parameters you will need to adjust**

*Information about the input data*

``! data_input_path``
    Directory where your concatenated target data are stored.
``! data_input_pattern``
    Regular expression pattern of all your target files.
    .. note::

        These files should have the direction-independent calibration applied to the DATA
        column (usually the ``*.pre-cal.ms`` files from the target pipeline).

*Location of the software*

``! prefactor_directory``
    Path to your prefactor copy
``! wsclean_executable``
    Path to your local WSClean executable

**Parameters you may need to adjust**

*Imaging and subtraction options*

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


Parameters for **HBA** and **LBA** observations
-----------------------------------------------

At this time, only HBA data are supported.
