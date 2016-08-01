import os
#from lofarpipe.support.data_map import DataMap
#from lofarpipe.support.data_map import DataProduct


def plugin_main(args, **kwargs):
    """
    Creates a directory on the disk.

    Parameters
    ----------
    directory : str
        Filename of datamap to trim

    Returns
    -------
    result : dict
        Empty driectory

    """
    dir_to_create = kwargs['directory']

    if os.path.exists(dir_to_create):
        if not os.path.isdir(dir_to_create):
            print "Plugin makeDirectory: Path {} exists, but is no directory!".format(dir_to_create)
            raise ValueError, "Path {} exists, but is no directory!".format(dir_to_create)
    else:
        os.makedirs(dir_to_create,mode=0755)

    result = {}
    return result
