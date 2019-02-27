.. _data_preparation:

Preparing the data
------------------

.. note::

    Processing of interleaved datasets is not currently supported.

**Prefactor** requires **LOFAR LBA** or **HBA** raw or pre-processed data. These data are
typically obtained from the LOFAR Long-Term Archive at https://lta.lofar.eu. All
input measurement-sets for one pipeline run need to be in the same directory.

- The calibrator and target data have to match, i.e., be observed close enough
  in time that calibration values can be transferred.

- For each observation you should process all the calibrator data at once
  together. Clock/TEC separation and flagging of bad amplitudes work better with
  the full bandwidth.

- For the target pipeline you will need to have internet access from the machine you are running **prefactor**.
  It is required in order to retrieve RM values from `CODE`_ and a global sky model (`TGSS`_ or `GSM`_). Both are hosted as online services.
  
  
.. _CODE: ftp://ftp.aiub.unibe.ch/CODE/
.. _TGSS: http://tgssadr.strw.leidenuniv.nl/doku.php
.. _GSM:  http://172.104.228.177/