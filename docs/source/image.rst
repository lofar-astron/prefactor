.. _image_pipeline:

Image pipeline
==============

This pipeline produces an image of the full FOV, using the full bandwidth. It is intended primarily
for use in production.

The options for this pipeline are as follows:

``! data_input_filenames``
    List of input MS filenames (full path).
``! image_output_filenames``
    List of output image filenames for feedback (full path).
``! cellsize_highres_deg``
    Cellsize in degrees.
``! fieldsize_highres``
    Size of the image is this value times the FWHM of mean semi-major axis of
    the station beam.
``! maxlambda_highres``
    Maximum uv-distance in lambda that will be used for imaging.
``! image_padding``
    How much padding shall we add during the imaging?
``! axis_stretch``
    How much shall the y-axis be stretched or compressed?
``! idg_mode``
    IDG mode to use: cpu or hybrid (= CPU + GPU).
``! local_scratch_dir``
    Scratch directory for wsclean (can be local to the processing nodes!).
``! images_metadata_file``
    Feedback metadata file (full path).
``! parset_prefix``
    Feedback parset prefix.
``! image_rootname``
    Output image root name. The image will be named ``image_rootname-MFS-image.fits``.


The pipeline consists of the following steps:

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
``average``
    Average the data with the values determined above. The amount of averaging
    depends on the size of the image (to limit bandwidth and time smearing).
``combine_mapfile_deep``
    Generate a mapfile with all files in a single entry. This mapfile is used as
    input to the next step.
``dpppconcat``
    Run DPPP to concatenate the data. Concatenating the data speeds up gridding
    and degridding with IDG by factors of several.
``wsclean_high_deep``
    Image the data with WSClean+IDG.
``copy_output_images``
    Copy the image to the output filenames expected for feedback.
``make_image_metadata``
    Make feedback metadata for the output image. This step is needed to make the
    info needed for ingest into the LTA.

