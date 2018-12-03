.. _image_pipeline:

Image pipeline
==============

This pipeline produces an image of the full FOV, using the full bandwidth. The
parset is named ``Pre-Facet-Image.parset``. It is intended primarily for use in
production.

Options
-------

``! data_input_filenames``
    List of input MS filenames (full path).
``! image_output_filenames``
    List of output image filenames for feedback (full path).
``! cellsize_highres_deg``
    Cellsize in degrees.
``! fieldsize_highres``
    Size of the image is this value times the FWHM of mean semi-major axis of
    the station beam at the lowest observed frequency.
``! maxlambda_highres``
    Maximum uv-distance in lambda that will be used for imaging.
``! image_padding``
    How much padding shall we add during the imaging?
``! idg_mode``
    IDG mode to use: cpu or hybrid (= CPU + GPU).
``! local_scratch_dir``
    Scratch directory for WSClean (can be local to the processing nodes!).
``! images_metadata_file``
    Feedback metadata file (full path).
``! parset_prefix``
    Feedback parset prefix.
``! image_rootname``
    Output image root name. The image will be named ``image_rootname-MFS-image.fits``.


Steps
-----

``create_ms_map``
    Generate a mapfile of all the target data. The files must be supplied as a
    list of the full paths to the files.
``combine_mapfile``
    Generate a mapfile with all files in a single entry. This mapfile is used as
    input to the next step.
``do_magic``
    Compute image sizes and the number of channels to use during imaging from the MS
    files from the previous step. The image size is calculated from the FWHM of the
    primary beam at the lowest frequency at the mean elevation of the observation. The
    number of channels is set simply as the number of subbands / 40, to result in
    enough channels to allow multi-frequency synthesis (MFS), but not so many that
    performance is impacted. A minimum of 2 channels is used.
``do_magic_maps``
    Convert the output of do_magic into usable mapfiles.
``average``
    Average the data as appropriate for imaging of the FOV. The amount of averaging
    depends on the size of the image (to limit bandwidth and time smearing). The
    averaging currently adopted is 16 s per time slot and 0.2 MHz per channel. These
    values result in low levels of bandwidth and time smearing for the target image
    sizes and resolutions.
``combine_mapfile_deep``
    Generate a mapfile with all files in a single entry. This mapfile is used as
    input to the next step.
``dpppconcat``
    Run DPPP to concatenate the data. Concatenating the data speeds up gridding
    and degridding with IDG by factors of several.
``wsclean_high_deep``
    Image the data with WSClean+IDG. Imaging is done in MFS mode, resulting in a
    single image for the full bandwidth. A typical HBA image looks like the one below.

    .. image:: image_pipeline_example.png

``copy_output_images``
    Copy the image to the output filenames expected for feedback.
``make_image_metadata``
    Make feedback metadata for the output image. This step is needed to make the
    info needed for ingest into the LTA. This info includes the RA, Dec of the image
    center, an estimate of the image rms noise, the frequency, etc.

