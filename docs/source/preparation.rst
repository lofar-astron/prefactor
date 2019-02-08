.. _data_preparation:

Preparing the data
------------------

.. note::

    Processing of interleaved datasets is not currently supported.

Prefactor requires LOFAR LBA or HBA raw or pre-processed data. These data are
typically obtained from the LOFAR Long-Term Archive at https://lta.lofar.eu. All
input measurement-sets for one pipeline run need to be in the same directory.

- The calibrator and target data have to match, i.e., be observed close enough
  in time that calibration values can be transferred.

- For each observation you should process all the calibrator data at once
  together. Clock/TEC separation and flagging of bad amplitudes work better with
  the full bandwidth.

- The target data can be processed one time-/frequency- block at a time. This
  strategy can be useful if you have limited disk space (the processed data is a
  lot smaller than the input data). You do need to process all data that should be
  concatenated into one file in one go. (Otherwise the pipeline will fail, or
  you'll get two output files with overlapping frequency coverage.)


.. _genericpipeline: http://www.astron.nl/citt/genericpipeline/
