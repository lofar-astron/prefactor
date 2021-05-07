#! /usr/bin/env python3
"""
Script to sort a list of MSs into frequency-bands, and compute additional values needed for initsubtract
"""
import casacore.tables as pt
import sys, os
import numpy as np
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct
import argparse
from argparse import RawTextHelpFormatter

class Band(object):
    """
    The Band object contains parameters needed for each band
    """
    def __init__(self, MSfiles):
        self.files = MSfiles
        self.msnames = [ MS.split('/')[-1] for MS in self.files ]
        self.numMS = len(self.files)
        # Get the frequency info and set name
        sw = pt.table(self.files[0]+'::SPECTRAL_WINDOW', ack=False)
        self.freq = sw.col('REF_FREQUENCY')[0]
        self.nchan = sw.col('NUM_CHAN')[0]
        self.chan_freqs_hz = sw.col('CHAN_FREQ')[0]
        self.chan_width_hz = sw.col('CHAN_WIDTH')[0][0]
        sw.close()
        self.name = str(int(self.freq/1e6))
        # Get the station diameter
        ant = pt.table(self.files[0]+'::ANTENNA', ack=False)
        self.diam = float(ant.col('DISH_DIAMETER')[0])
        ant.close()

    def get_image_sizes(self, cellsize_highres_deg=None, cellsize_lowres_deg=None,
                        fieldsize_highres=2.5, fieldsize_lowres=6.5):
        """
        Sets sizes for initsubtract images

        The image sizes are scaled from the mean primary-beam FWHM. For
        the high-res image, we use 2.5 * FWHM; for low-res, we use 6.5 * FHWM.

        Parameters
        ----------
        cellsize_highres_deg : float, optional
            cellsize for the high-res images in deg
        cellsize_lowres_deg : float, optional
            cellsize for the low-res images in deg
        fieldsize_highres : float, optional
            How many FWHM's shall the high-res images be.
        fieldsize_lowres : float, optional
            How many FWHM's shall the low-res images be.
        """
        if cellsize_highres_deg:
            self.cellsize_highres_deg = cellsize_highres_deg
        if cellsize_lowres_deg:
            self.cellsize_lowres_deg = cellsize_lowres_deg
        if not hasattr(self, 'mean_el_rad'):
            for MS_id in range(self.numMS):
                # calculate mean elevation
                tab = pt.table(self.files[MS_id], ack=False)
                el_values = pt.taql("SELECT mscal.azel1()[1] AS el from "
                                    + self.files[MS_id] + " limit ::10000").getcol("el")
                if MS_id == 0:
                    global_el_values = el_values
                else:
                    global_el_values = np.hstack( (global_el_values, el_values ) )
            self.mean_el_rad = np.mean(global_el_values)

        # Calculate mean FOV
        sec_el = 1.0 / np.sin(self.mean_el_rad)
        self.fwhm_deg = 1.1 * ((3.0e8 / self.freq) / self.diam) * 180. / np.pi * sec_el
        self.imsize_high_res = self.get_optimum_size(self.fwhm_deg
                                                     /self.cellsize_highres_deg * fieldsize_highres)
        self.imsize_low_res = self.get_optimum_size(self.fwhm_deg
                                                    /self.cellsize_lowres_deg * fieldsize_lowres)
        return (self.imsize_high_res, self.imsize_low_res)

    def get_optimum_size(self, size):
        """
        Gets the nearest optimum image size

        Taken from the casa source code (cleanhelper.py)

        Parameters
        ----------
        size : int
            Target image size in pixels

        Returns
        -------
        optimum_size : int
            Optimum image size nearest to target size

        """
        import numpy

        def prime_factors(n, douniq=True):
            """ Return the prime factors of the given number. """
            factors = []
            lastresult = n
            sqlast=int(numpy.sqrt(n))+1
            if n == 1:
                return [1]
            c=2
            while 1:
                 if (lastresult == 1) or (c > sqlast):
                     break
                 sqlast=int(numpy.sqrt(lastresult))+1
                 while 1:
                     if(c > sqlast):
                         c=lastresult
                         break
                     if lastresult % c == 0:
                         break
                     c += 1
                 factors.append(c)
                 lastresult = lastresult // c
            if (factors==[]): factors=[n]
            return  numpy.unique(factors).tolist() if douniq else factors

        n = int(size)
        if (n%2 != 0):
            n+=1
        fac=prime_factors(n, False)
        for k in range(len(fac)):
            if (fac[k] > 7):
                val=fac[k]
                while (numpy.max(prime_factors(val)) > 7):
                    val +=1
                fac[k]=val
        newlarge=numpy.product(fac)
        for k in range(n, newlarge, 2):
            if ((numpy.max(prime_factors(k)) < 8)):
                return k
        return newlarge

    def get_averaging_steps(self):
        """
        Sets the averaging step sizes

        Note: the frequency step must be an even divisor of the number of
        channels

        """
        # Get time per sample and number of samples
        t = pt.table(self.files[0], readonly=True, ack=False)
        for t2 in t.iter(["ANTENNA1","ANTENNA2"]):
            if (t2.getcell('ANTENNA1',0)) < (t2.getcell('ANTENNA2',0)):
                self.timestep_sec = t2[1]['TIME']-t2[0]['TIME'] # sec
                break
        t.close()
        # generate a (numpy-)array with the divisors of nchan
        tmp_divisors = []
        for step in range(self.nchan,0,-1):
            if (self.nchan % step) == 0:
                tmp_divisors.append(step)
        freq_divisors = np.array(tmp_divisors)

        # Average to 0.2 MHz per channel and 16 sec per time slot
        initsubtract_freqstep = max(1, min(int(round(0.2 * 1e6 / self.chan_width_hz)), self.nchan))
        initsubtract_freqstep = freq_divisors[np.argmin(np.abs(freq_divisors-initsubtract_freqstep))]
        initsubtract_timestep = max(1, int(round(16.0 / self.timestep_sec)))

        return (initsubtract_freqstep, initsubtract_timestep)

    def get_nwavelengths(self, cellsize_deg, timestep_sec,):
        """
        Returns nwavelengths for WSClean BL-based averaging
        The value depends on the integration time given the specified maximum
        allowed smearing. We scale it from the imaging cell size assuming normal
        sampling as:
        max baseline in nwavelengths = 1 / theta_rad ~= 1 / (cellsize_deg * 3 * pi / 180)
        nwavelengths = max baseline in nwavelengths * 2 * pi *
            integration time in seconds / (24 * 60 * 60) / 4
        Parameters
        ----------
        cellsize_deg : float
            Pixel size of image in degrees
        timestep_sec : float
            Length of one timestep in seconds
        """
        max_baseline = 1 / (3 * cellsize_deg * np.pi / 180)
        wsclean_nwavelengths_time = int(max_baseline * 2*np.pi * timestep_sec /
            (24 * 60 * 60) / 4)
        return wsclean_nwavelengths_time


def main(ms_input, outmapname=None, mapfile_dir=None, cellsize_highres_deg=0.00208, cellsize_lowres_deg=0.00694,
         fieldsize_highres=2.5, fieldsize_lowres=6.5, image_padding=1., y_axis_stretch=1.):
    """
    Check a list of MS files for missing frequencies

    Parameters
    ----------
    ms_input : list or str
        List of MS filenames, or string with list, or path to a mapfile
    outmapname: str
        Name of output mapfile
    mapfile_dir : str
        Directory for output mapfile
    cellsize_highres_deg : float, optional
        cellsize for the high-res images in deg
    cellsize_lowres_deg : float, optional
        cellsize for the low-res images in deg
    fieldsize_highres : float, optional
        How many FWHM's shall the high-res images be.
    fieldsize_lowres : float, optional
        How many FWHM's shall the low-res images be.
    image_padding : float, optional
        How much padding shall we add to the padded image sizes.
    y_axis_stretch : float, optional
        How much shall the y-axis be stretched or compressed.

    Returns
    -------
    result : dict
        Dict with the name of the generated mapfiles

    """

    if not outmapname or not mapfile_dir:
        raise ValueError('sort_times_into_freqGroups: outmapname and mapfile_dir are needed!')
    if type(ms_input) is str:
        if ms_input.startswith('[') and ms_input.endswith(']'):
            ms_list = [f.strip(' \'\"') for f in ms_input.strip('[]').split(',')]
        else:
            map_in = DataMap.load(ms_input)
            map_in.iterator = DataMap.SkipIterator
            ms_list = []
            for fname in map_in:
                if fname.startswith('[') and fname.endswith(']'):
                    for f in fname.strip('[]').split(','):
                        ms_list.append(f.strip(' \'\"'))
                else:
                    ms_list.append(fname.strip(' \'\"'))
    elif type(ms_input) is list:
        ms_list = [str(f).strip(' \'\"') for f in ms_input]
    else:
        raise TypeError('sort_into_freqBands: type of "ms_input" unknown!')

    cellsize_highres_deg = float(cellsize_highres_deg)
    cellsize_lowres_deg = float(cellsize_lowres_deg)
    fieldsize_highres = float(fieldsize_highres)
    fieldsize_lowres = float(fieldsize_lowres)
    image_padding = float(image_padding)
    y_axis_stretch = float(y_axis_stretch)

    msdict = {}
    for ms in ms_list:
        # group all MSs by frequency
        sw = pt.table(ms+'::SPECTRAL_WINDOW', ack=False)
        msfreq = int(sw.col('REF_FREQUENCY')[0])
        sw.close()
        if msfreq in msdict:
            msdict[msfreq].append(ms)
        else:
            msdict[msfreq] = [ms]
    bands = []
    bandfreqs = []
    print("InitSubtract_deep_sort_and_compute.py: Putting files into bands.")
    for MSkey in msdict.keys():
        bands.append( Band(msdict[MSkey]) )
        bandfreqs.append( Band(msdict[MSkey]).freq )

    ## min freq gives largest image size for deep image
    bandfreqs = np.array(bandfreqs)
    minfreq = np.min(bandfreqs)
    bandmin = np.argmin(bandfreqs)
    ## need to map the output from wsclean channels to the right frequencies
    ## just put the bands in the right freq order
    wsclean_channum = np.argsort(bandfreqs)
    bands = np.array(bands)
    bands = bands[wsclean_channum]

    #minfreq = 1e9
    #for ib, band in enumerate(bands):
        #if band.freq < minfreq:
            #minfreq = band.freq
            #bandmin = ib

    group_map = MultiDataMap()
    file_single_map = DataMap([])
    high_size_map = DataMap([])
    low_size_map = DataMap([])
    high_paddedsize_map = DataMap([])
    low_paddedsize_map = DataMap([])
    numfiles = 0
    nbands = np.int(len(bands)/10)
    if nbands > 8:
        nchansout_clean1 = np.int(nbands/4)
    elif nbands > 4:
        nchansout_clean1 = np.int(nbands/2)
    else:
        nchansout_clean1 = np.int(nbands)
    if nchansout_clean1 < 2:
        nchansout_clean1 = 2

    (freqstep, timestep) = bands[0].get_averaging_steps()
    int_time_sec = bands[0].timestep_sec * timestep   # timestep_sec gets added to band object in get_averaging_steps()
    nwavelengths_high = bands[0].get_nwavelengths(cellsize_highres_deg, int_time_sec)
    nwavelengths_low = bands[0].get_nwavelengths(cellsize_lowres_deg, int_time_sec)

    print("InitSubtract_deep_sort_and_compute.py: analyzing data...")
    for band in bands:
        group_map.append(MultiDataProduct('localhost', band.files, False))
        numfiles += len(band.files)
        for filename in band.files:
            file_single_map.append(DataProduct('localhost', filename, False))
        (imsize_high_res, imsize_low_res) = band.get_image_sizes(cellsize_highres_deg, cellsize_lowres_deg,
                                                                 fieldsize_highres, fieldsize_lowres)
        imsize_high_res_stretch = band.get_optimum_size(int(imsize_high_res*y_axis_stretch))
        high_size_map.append(DataProduct('localhost', str(imsize_high_res)+" "+str(imsize_high_res_stretch), False))
        imsize_low_res_stretch = band.get_optimum_size(int(imsize_low_res*y_axis_stretch))
        low_size_map.append(DataProduct('localhost', str(imsize_low_res)+" "+str(imsize_low_res_stretch), False))
        imsize_high_pad = band.get_optimum_size(int(imsize_high_res*image_padding))
        imsize_high_pad_stretch = band.get_optimum_size(int(imsize_high_res*image_padding*y_axis_stretch))
        high_paddedsize_map.append(DataProduct('localhost', str(imsize_high_pad)+" "+str(imsize_high_pad_stretch), False))
        imsize_low_pad = band.get_optimum_size(int(imsize_low_res*image_padding))
        imsize_low_pad_stretch = band.get_optimum_size(int(imsize_low_res*image_padding*y_axis_stretch))
        low_paddedsize_map.append(DataProduct('localhost', str(imsize_low_pad)+" "+str(imsize_low_pad_stretch), False))

        if band.freq == minfreq:
            deep_imsize_high_res = imsize_high_res
            deep_imsize_high_res_stretch = imsize_high_res_stretch
            deep_imsize_high_pad = imsize_high_pad
            deep_imsize_high_pad_stretch = imsize_high_pad_stretch
            deep_imsize_low_res = imsize_low_res
            deep_imsize_low_res_stretch = imsize_low_res_stretch
            deep_imsize_low_pad = imsize_low_pad
            deep_imsize_low_pad_stretch = imsize_low_pad_stretch

    deep_high_size_map = DataMap([DataProduct('localhost', str(deep_imsize_high_res)+" "+str(deep_imsize_high_res_stretch), False)])
    deep_high_paddedsize_map = DataMap([DataProduct('localhost', str(deep_imsize_high_pad)+" "+str(deep_imsize_high_pad_stretch), False)])
    deep_low_size_map = DataMap([DataProduct('localhost', str(deep_imsize_low_res)+" "+str(deep_imsize_low_res_stretch), False)])
    deep_low_paddedsize_map = DataMap([DataProduct('localhost', str(deep_imsize_low_pad)+" "+str(deep_imsize_low_pad_stretch), False)])
    nbands_map = DataMap([DataProduct('localhost', str(nbands), False)])
    nchansout_clean1_map = DataMap([DataProduct('localhost', str(nchansout_clean1), False)])

    # get mapfiles for freqstep and timestep with the length of single_map
    freqstep_map = DataMap([])
    timestep_map = DataMap([])
    nwavelengths_high_map        = DataMap([])
    nwavelengths_low_map         = DataMap([])

    for index in range(numfiles):
        # set time and frequency averaging values (in sec and Hz)
        freqstep_map.append(DataProduct('localhost', str(freqstep*bands[0].chan_width_hz), False))
        timestep_map.append(DataProduct('localhost', str(timestep*bands[0].timestep_sec), False))
    nwavelengths_high_map.append(DataProduct('localhost', str(nwavelengths_high), False))
    nwavelengths_low_map.append(DataProduct('localhost', str(nwavelengths_low), False))

    groupmapname = os.path.join(mapfile_dir, outmapname)
    group_map.save(groupmapname)
    file_single_mapname = os.path.join(mapfile_dir, outmapname+'_single')
    file_single_map.save(file_single_mapname)

    high_sizename = os.path.join(mapfile_dir, outmapname+'_high_size')
    high_size_map.save(high_sizename)
    low_sizename = os.path.join(mapfile_dir, outmapname+'_low_size')
    low_size_map.save(low_sizename)
    high_padsize_name = os.path.join(mapfile_dir, outmapname+'_high_padded_size')
    high_paddedsize_map.save(high_padsize_name)
    low_padsize_name = os.path.join(mapfile_dir, outmapname+'_low_padded_size')
    low_paddedsize_map.save(low_padsize_name)

    deep_high_sizename = os.path.join(mapfile_dir, outmapname+'_deep_high_size')
    deep_high_size_map.save(deep_high_sizename)
    deep_low_sizename = os.path.join(mapfile_dir, outmapname+'_deep_low_size')
    deep_low_size_map.save(deep_low_sizename)
    deep_high_padsize_name = os.path.join(mapfile_dir, outmapname+'_deep_high_padded_size')
    deep_high_paddedsize_map.save(deep_high_padsize_name)
    deep_low_padsize_name = os.path.join(mapfile_dir, outmapname+'_deep_low_padded_size')
    deep_low_paddedsize_map.save(deep_low_padsize_name)

    nbands_mapname = os.path.join(mapfile_dir, outmapname+'_nbands')
    nbands_map.save(nbands_mapname)
    nchansout_clean1_mapname = os.path.join(mapfile_dir, outmapname+'_nchansout_clean1')
    nchansout_clean1_map.save(nchansout_clean1_mapname)

    freqstepname = os.path.join(mapfile_dir, outmapname+'_freqstep')
    freqstep_map.save(freqstepname)
    timestepname = os.path.join(mapfile_dir, outmapname+'_timestep')
    timestep_map.save(timestepname)
    nwavelengths_high_name = os.path.join(mapfile_dir, outmapname+'_nwavelengths_high')
    nwavelengths_high_map.save(nwavelengths_high_name)
    nwavelengths_low_name = os.path.join(mapfile_dir, outmapname+'_nwavelengths_low')
    nwavelengths_low_map.save(nwavelengths_low_name)

    result = {'groupmap': groupmapname, 'single_mapfile' : file_single_mapname,
              'high_size_mapfile' : high_sizename, 'low_size_mapfile' : low_sizename,
              'high_padsize_mapfile' : high_padsize_name, 'low_padsize_mapfile' : low_padsize_name,
              'deep_high_size_mapfile' : deep_high_sizename, 'deep_low_size_mapfile' : deep_low_sizename,
              'deep_high_padsize_mapfile' : deep_high_padsize_name, 'deep_low_padsize_mapfile' : deep_low_padsize_name,
              'nbands' : nbands_mapname, 'nchansout_clean1' : nchansout_clean1_mapname,
              'freqstep' : freqstepname, 'timestep' : timestepname, 'nwavelengths_high_mapfile': nwavelengths_high_name,
              'nwavelengths_low_mapfile': nwavelengths_low_name}
    return result


class MultiDataProduct(DataProduct):
    """
    Class representing multiple files in a DataProduct.
    """
    def __init__(self, host=None, file=None, skip=True):
        super(MultiDataProduct, self).__init__(host, file, skip)
        if not file:
            self.file = list()
        else:
            self._set_file(file)

    def __repr__(self):
        """Represent an instance as a Python dict"""
        return (
            "{{'host': '{0}', 'file': {1}, 'skip': {2}}}".format(self.host, self.file, str(self.skip))
        )

    def __str__(self):
        """Print an instance as 'host:[filelist]'"""
        return ':'.join((self.host, str(self.file)))

    def _set_file(self, data):
        try:
            # Try parsing as a list
            if isinstance(data, list):
                self.file = data
            if isinstance(data, DataProduct):
                self._from_dataproduct(data)
            if isinstance(data, DataMap):
                self._from_datamap(data)

        except TypeError:
            raise DataProduct("No known method to set a filelist from %s" % str(file))

    def _from_dataproduct(self, prod):
        self.host = prod.host
        self.file = prod.file
        self.skip = prod.skip

    def _from_datamap(self, inmap):
        filelist = {}
        for item in inmap:
            if not item.host in filelist:
                filelist[item.host] = []
            filelist[item.host].append(item.file)
        self.file = filelist['i am']

    def append(self, item):
        self.file.append(item)


class MultiDataMap(DataMap):
    """
    Class representing a specialization of data-map, a collection of data
    products located on the same node, skippable as a set and individually
    """
    @DataMap.data.setter
    def data(self, data):
        if isinstance(data, DataMap):
            mdpdict = {}
            data.iterator = DataMap.SkipIterator
            for item in data:
                if not item.host in mdpdict:
                    mdpdict[item.host] = []
                mdpdict[item.host].append(item.file)
            mdplist = []
            for k, v in mdpdict.iteritems():
                mdplist.append(MultiDataProduct(k, v, False))
            self._set_data(mdplist, dtype=MultiDataProduct)
        elif isinstance(data, MultiDataProduct):
            self._set_data(data, dtype=MultiDataProduct)
        elif not data:
            pass
        else:
            self._set_data(data, dtype=MultiDataProduct)

    def split_list(self, number):
        mdplist = []
        for item in self.data:
            for i in range(0, len(item.file), number):
                chunk = item.file[i:i+number]
                mdplist.append(MultiDataProduct(item.host, chunk, item.skip))
        self._set_data(mdplist)
