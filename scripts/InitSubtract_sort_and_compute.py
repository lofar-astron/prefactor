#! /usr/bin/env python3
"""
Script to sort a list of MSs into frequency-bands, and compute additional values needed for initsubtract
"""
import pyrap.tables as pt
import os
import numpy as np
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct

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
        return(self.imsize_high_res, self.imsize_low_res)

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
                return([1])
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
                 lastresult /= c
            if (factors==[]): factors=[n]
            return(numpy.unique(factors).tolist() if douniq else factors)

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
                return(k)
        return(newlarge)

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

        # For initsubtract, average to 0.5 MHz per channel and 20 sec per time
        # slot. Since each band is imaged separately and the smearing and image
        # sizes both scale linearly with frequency, a single frequency and time
        # step is valid for all bands
        initsubtract_freqstep = max(1, min(int(round(0.5 * 1e6 / self.chan_width_hz)), self.nchan))
        initsubtract_freqstep = freq_divisors[np.argmin(np.abs(freq_divisors-initsubtract_freqstep))]
        initsubtract_timestep = max(1, int(round(20.0 / self.timestep_sec)))

        return(initsubtract_freqstep, initsubtract_timestep)

    def nwavelengths(self, cellsize_highres_deg, cellsize_lowres_deg, initsubtract_timestep):
        int_time_sec = self.timestep_sec * initsubtract_timestep
        max_baseline_in_nwavelenghts_h = 1.0/(cellsize_highres_deg*3.0*np.pi/180.0)
        max_baseline_in_nwavelenghts_l = 1.0/(cellsize_lowres_deg*3.0*np.pi/180.0)
        self.nwavelengths_high	=	max_baseline_in_nwavelenghts_h*2.0*np.pi*int_time_sec/(24.0*60.0*60.0)/2
        self.nwavelengths_low	=	max_baseline_in_nwavelenghts_l*2.0*np.pi*int_time_sec/(24.0*60.0*60.0)/2
        return(self.nwavelengths_high, self.nwavelengths_low)


def input2bool(invar):
    if invar is None:
        return(None)
    if isinstance(invar, bool):
        return(invar)
    elif isinstance(invar, str):
        if invar.upper() == 'TRUE' or invar == '1':
            return(True)
        elif invar.upper() == 'FALSE' or invar == '0':
            return(False)
        else:
            raise(ValueError('input2bool: Cannot convert string "'+invar+'" to boolean!'))
    elif isinstance(invar, int) or isinstance(invar, float):
        return(bool(invar))
    else:
        raise(TypeError('input2bool: Unsupported data type:'+str(type(invar))))


def main(ms_input, outmapname=None, mapfile_dir=None, cellsize_highres_deg=0.00208, cellsize_lowres_deg=0.00694,
         fieldsize_highres=2.5, fieldsize_lowres=6.5, image_padding=1., y_axis_stretch=1.,
         calc_y_axis_stretch=False, apply_y_axis_stretch_highres=True, apply_y_axis_stretch_lowres=True):
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
    calc_y_axis_stretch : bool, optional
        Adjust the image sizes returned by this script for the mean elevation.
        If True, the value of y_axis_stretch above is ignored
    apply_y_axis_stretch_highres : bool, optional
        Apply the y-axis stretch to the high-res image sizes
    apply_y_axis_stretch_lowres : bool, optional
        Apply the y-axis stretch to the low-res image sizes

    Returns
    -------
    result : dict
        Dict with the name of the generated mapfiles

    """
    if not outmapname or not mapfile_dir:
        raise(ValueError('sort_times_into_freqGroups: outmapname and mapfile_dir are needed!'))
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
        raise(TypeError('sort_into_freqBands: type of "ms_input" unknown!'))

    cellsize_highres_deg = float(cellsize_highres_deg)
    cellsize_lowres_deg = float(cellsize_lowres_deg)
    fieldsize_highres = float(fieldsize_highres)
    fieldsize_lowres = float(fieldsize_lowres)
    image_padding = float(image_padding)
    y_axis_stretch = float(y_axis_stretch)
    calc_y_axis_stretch = input2bool(calc_y_axis_stretch)
    apply_y_axis_stretch_highres = input2bool(apply_y_axis_stretch_highres)
    apply_y_axis_stretch_lowres = input2bool(apply_y_axis_stretch_lowres)

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
    print("InitSubtract_sort_and_compute.py: Putting files into bands.")
    for MSkey in msdict.keys():
        bands.append( Band(msdict[MSkey]) )

    group_map = MultiDataMap()
    file_single_map = DataMap([])
    high_size_map = DataMap([])
    low_size_map = DataMap([])
    high_paddedsize_map = DataMap([])
    low_paddedsize_map = DataMap([])
    numfiles = 0
    for i, band in enumerate(bands):
        print("InitSubtract_sort_and_compute.py: Working on Band:",band.name)
        group_map.append(MultiDataProduct('localhost', band.files, False))
        numfiles += len(band.files)
        for filename in band.files:
            file_single_map.append(DataProduct('localhost', filename, False))
        (imsize_high_res, imsize_low_res) = band.get_image_sizes(cellsize_highres_deg, cellsize_lowres_deg,
                                                                 fieldsize_highres, fieldsize_lowres)

        # Calculate y_axis_stretch if desired
        if calc_y_axis_stretch:
            if i == 0:
                y_axis_stretch = 1.0 / np.sin(band.mean_el_rad)
                print("InitSubtract_sort_and_compute.py: Using y-axis stretch of:",y_axis_stretch)
                y_axis_stretch_lowres = y_axis_stretch
                y_axis_stretch_highres = y_axis_stretch
        else:
            y_axis_stretch_lowres = 1.0
            y_axis_stretch_highres = 1.0

        # Adjust sizes so that we get the correct ones below
        if not apply_y_axis_stretch_highres:
            y_axis_stretch_highres = 1.0
        if not apply_y_axis_stretch_lowres:
            y_axis_stretch_lowres = 1.0
        imsize_low_res /= y_axis_stretch_lowres
        imsize_high_res /= y_axis_stretch_highres

        imsize_high_res = band.get_optimum_size(int(imsize_high_res))
        imsize_high_res_stretch = band.get_optimum_size(int(imsize_high_res*y_axis_stretch_highres))
        high_size_map.append(DataProduct('localhost', str(imsize_high_res)+" "+str(imsize_high_res_stretch), False))

        imsize_low_res = band.get_optimum_size(int(imsize_low_res))
        imsize_low_res_stretch = band.get_optimum_size(int(imsize_low_res*y_axis_stretch_lowres))
        low_size_map.append(DataProduct('localhost', str(imsize_low_res)+" "+str(imsize_low_res_stretch), False))

        imsize_high_pad = band.get_optimum_size(int(imsize_high_res*image_padding))
        imsize_high_pad_stretch = band.get_optimum_size(int(imsize_high_res*image_padding*y_axis_stretch_highres))
        high_paddedsize_map.append(DataProduct('localhost', str(imsize_high_pad)+" "+str(imsize_high_pad_stretch), False))

        imsize_low_pad = band.get_optimum_size(int(imsize_low_res*image_padding))
        imsize_low_pad_stretch = band.get_optimum_size(int(imsize_low_res*image_padding*y_axis_stretch_lowres))
        low_paddedsize_map.append(DataProduct('localhost', str(imsize_low_pad)+" "+str(imsize_low_pad_stretch), False))

    print("InitSubtract_sort_and_compute.py: Computing averaging steps.")
    (freqstep, timestep) = bands[0].get_averaging_steps()
    (nwavelengths_high, nwavelengths_low) = bands[0].nwavelengths(cellsize_highres_deg,
        cellsize_lowres_deg, timestep)

    # get mapfiles for freqstep and timestep with the length of single_map
    freqstep_map = DataMap([])
    timestep_map = DataMap([])
    nwavelengths_high_map = DataMap([])
    nwavelengths_low_map = DataMap([])
    for index in range(numfiles):
        freqstep_map.append(DataProduct('localhost', str(freqstep), False))
        timestep_map.append(DataProduct('localhost', str(timestep), False))
        freqstep_map.append(DataProduct('localhost', str(freqstep), False))
        timestep_map.append(DataProduct('localhost', str(timestep), False))
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
              'freqstep' : freqstepname, 'timestep' : timestepname, 'nwavelengths_high_mapfile': nwavelengths_high_name,
              'nwavelengths_low_mapfile': nwavelengths_low_name}
    return(result)


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
            "{{'host': '{0}', 'file': '{1}', 'skip': {2}}}".format(self.host,
                '[{}]'.format(','.join(self.file)), str(self.skip))
        )

    def __str__(self):
        """Print an instance as 'host:[filelist]'"""
        return(':'.join((self.host, str(self.file))))

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
            raise(DataProduct("No known method to set a filelist from %s" % str(file)))

    def _from_dataproduct(self, prod):
        print('setting filelist from DataProduct')
        self.host = prod.host
        self.file = prod.file
        self.skip = prod.skip

    def _from_datamap(self, inmap):
        print('setting filelist from DataMap')
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
    def __init__(self, data=list(), iterator=iter):
        super(MultiDataMap, self).__init__(data, iterator)

    @classmethod
    def load(cls, filename):
        """Load a data map from file `filename`. Return a DataMap instance."""
        with open(filename) as f:
            datamap = eval(f.read())
            for i, d in enumerate(datamap):
                file_entry = d['file']
                if file_entry.startswith('[') and file_entry.endswith(']'):
                    file_list = [e.strip(' \'\"') for e in file_entry.strip('[]').split(',')]
                    datamap[i] = {'host': d['host'], 'file': file_list, 'skip': d['skip']}
            return(cls(datamap))

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
