.. _runfactor:

Starting a pipeline
-------------------

Once you have the data and the parsets ready, you can run the pipeline using the
genericpipeline script, e.g.::

    $ genericpipeline.py -d -c pipeline.cfg My_prefactor_calibrator.parset

.. note::

    The -d option is recommended: it does make the log-files extremely large
    (many megabytes), but without it, often the important information as to why a
    pipeline run fails is not included.

While the pipeline runs, in the specified ``runtime_directory`` (see previous
section) new files are generated in a directory named after the parset name (e.g.,
if you are running ``My_prefactor_calibrator.parset`` a directory named
``My_prefactor_calibrator`` will appear in your ``runtime_directory``)::

    $ ls My_prefactor_calibrator/
    logs mapfiles parsets statefile

The logs dir contains all the logs of the pipeline runs, identified by the date
and time of execution, e.g.::

    $ ls My_prefactor_calibrator/logs/2016-06-30T15:07:21/pipeline.log

These contain all the output from the processes the pipeline called and
diagnostic information about the pipeline. So they are useful to follow the
status of the process and possibly identify reasons why a process crashed.

While running the pipeline writes a statefile in the ``runtime_directory``, with all
the step which were successfully executed. If the pipeline stops for whatever
reason, you can re-run the same command and it will skip all the steps that are
already done and only work on those which are still missing.

The intermediate data files of the pipeline are written in the ``working_directory``
specified in the ``pipeline.cfg``.


Stopping and restarting the pipeline
------------------------------------

You can stop a pipeline run anytime by terminating the genericpipeline process
(typically by pressing CRTL-C in the terminal where you started it). Sometimes some of
the processes that the pipeline started don't get properly terminated, so if the
genericpipeline process doesn't terminate you should look for its child
processes and terminate them too.

.. note::

    If you stop and re-start pipelines a number of time then you should also
    check occasionally if there are orphaned children that are eating up
    resources on you computer.

As mentioned earlier, you can re-start the pipeline by running the same command
with which you started it.


Pipeline crashes
----------------

It can happen that the pipeline stops with a message like this::

     ERROR   genericpipeline: LOFAR Pipeline finished unsuccesfully.
     WARNING genericpipeline: recipe genericpipeline completed with errors

You need to read the log of that run to identify the reason why it stopped, e.g.::

    > less My_prefactor_calibrator/logs/2016-06-30T15:07:21/pipeline.log

It is usually best to first check at the end of the file for what ended the
pipeline and then search from the beginning of the file for error or diagnostic
messages that tell you what exactly went wrong. See :ref:`help` for tips on
interpreting the error messages.

If you identify the problem and it does not affect the products that have been
already produced, you can launch the pipeline again, after correcting the issue
causing the process to stop.


Rerunning parts of the pipeline
--------------------------------

You can fully rerun a pipeline by deleting the runtime and working directories and restarting the pipeline.

To rerun parts of the pipeline that were (allegedly) already executed
successfully, you need to modify the ``statefile`` of the pipeline. To do this
there is a ``statefile_manipulation.py`` script as part of prefactor::

    python prefactor/bin/statefile_manipulation.py My_Workdir/My_calibrator_job/statefile

If you then run the pipeline again, it will start at the step that you removed with the statefile manipulation tool.
