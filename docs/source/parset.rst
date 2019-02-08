.. _parset:

Configuring prefactor
=====================

Preparing the configuration file
--------------------------------
To set up the genericpipeline for prefactor you need to customize the ``pipeline.cfg`` configuration file:

- Copy ``$LOFARROOT/share/pipeline/pipeline.cfg`` to someplace in your ``$HOME`` and open it in an editor.

- It starts with a section ``[DEFAULT]``, in there you need to edit three entries:

  - ``runtime_directory``: This is the directory where the pipeline puts
    logfiles, parsets, the status of successful steps etc. This can be set to a
    directory in your ``$HOME``, but it is recommended to set this to the same
    value as the `working_directory`

  - ``working_directory``: This is the directory where the processed data of the
    intermediate and final steps is put. Set this to a directory on your data
    disk, e.g. ``/data/scratch/<username>/PipelineExample``

  - ``recipe_directories``: This is a list of directories where the pipeline
    searches for recipes and plugins. There should already be an entry there, but
    another needs to be added so that the plugin scripts of prefactor are found.
    You need to add the pre-facet calibration directory to this list (so that the
    ``plugins`` directory is a subdirectory of one of the
    ``recipe\_directories``). E.g.:
    ``recipe_directories = [%(pythonpath)s/lofarpipe/recipes,/home/<username>/software/prefactor]``.

- There may be empty entries in the ``[DEFAULT]`` section of the
  ``pipeline.cfg`` file. This was set during the installation of the LOFAR
  software and usually is no need to worry.

- In case you do not run it on a cluster, you need to tell the pipeline to start the processes on the local machine:

  - Search for the section ``[cluster]``, set the entry clusterdesc to ``(lofarroot)s/share/local.clusterdesc``.

  - Add a section ``[remote]`` by adding the following lines::

      [remote]
      method = local
      max_per_node = <NumCPUCores>

  If there is already another section ``[remote]``, then remove that.

- The pipeline framework contains a number of parallel processing schemes for
  working on multi-node clusters. Ask you local sysadmin for advice.


Preparing the pipeline parset
-----------------------------

The pipeline parsets available in the prefactor directory (those files ending
with ``.parset``) are templates of genericpipeline parsets and need some small
editing before they can be run. To avoid confusion you should make a copy of the
parset that you want to run and give it a descriptive name. All parameters that
need to be changed are defined at the top of the parset with lines beginning
with ``!``. See the comments in the pipeline parsets and the notes below for some
hints on setting these parameters.

Some parameters depend on the observation to be processed and need to be
modified for each new observation, others are more "machine-dependent" so they
are the same for different observations that are processed on the same
machine(s).

See :ref:`pipeline_overview` for an overview of the pipeline parsets, and their
respective pages for a more in-depth description. Below are some general guidelines
for preparing the parsets:

- Don't edit the original parset files directly. Make a copy with a descriptive
  name (e.g. ``Pre-Facet-Cal-calibrator-3c295.parset``) and edit that copy.

- The pipeline framework will use the filename of the pipeline parset as the
  job-name if the latter is not explicitly given. That way there is no need to
  change the ``runtime_directory`` and ``working_directory`` entries in the
  ``pipeline.cfg`` for different pipeline runs. The pipeline framework will
  generate sub-directories with the job-name in there.

- The ``reference_station`` should be a station that is present in the
  observation, and didn't do anything strange. It is worth checking the plots and
  possibly changing the reference station if problems are found.

- The ``num_proc_per_node`` parameter is the number of processes of "small"
  programs (with little memory usage and which are not multi-threaded) that are
  run in parallel on one node.

- The ``num_proc_per_node_limit`` parameter is the number of processes of "big"
  programs (with large memory usage and/or multi-threaded programs) that are run
  in parallel on one node.

- In addition to setting how many DPPP processes (DPPP is a "big" program) are
  run in parallel, you can set how many threads each DPPP process may use with the
  ``max_dppp_threads`` parameter.

- So both: ``num_proc_per_node`` and (``num_proc_per_node_limit`` *
  ``max_dppp_threads``) should be equal or smaller than the number of cores that
  you have on your processing computers.

- Similarly (``max_imagers_per_node`` * ``max_percent_mem_per_img``) should be
  less that 100% and (``max_imagers_per_node`` * ``max_cpus_per_img``) should be
  equal or smaller than the number of cores.

- If your pipeline runs out of memory, then you can also lower these parameters
  to make the pipeline use less memory.

- Most of the actual processing is now done in DPPP, so the parameters that
  control its behavior are the important ones.
