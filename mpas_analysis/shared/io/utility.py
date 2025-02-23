# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
"""
IO utility functions

Phillip J. Wolfram, Xylar Asay-Davis
"""

import glob
import os
import random
import string
from datetime import datetime
import numpy
import shutil
import warnings


def paths(*args):
    """
    Returns glob'd paths in list for arbitrary number of function arguments.
    Note, each expanded set of paths is sorted.

    Parameters
    ----------
    *args : list
        A list of arguments to pass to ``glob.glob``

    Returns
    -------
    paths : list of str
        A list of file paths
    """
    # Authors
    # -------
    # Phillip J. Wolfram

    paths = []
    for aargs in args:
        paths += sorted(glob.glob(aargs))
    return paths


def fingerprint_generator(size=12,
                          chars=string.ascii_uppercase + string.digits):
    """
    Returns a random string that can be used as a unique fingerprint

    Parameters
    ----------
    size : int, optional
        The number of characters in the fingerprint

    chars : list of char, optional
        The fingerprint

    Returns
    -------
    fingerprint : str
        A random string

    Reference
    ---------
    http://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits-in-python
    """
    # Authors
    # -------
    # Phillip J. Wolfram

    return ''.join(random.choice(chars) for _ in range(size))


def make_directories(path):
    """
    Make the given path if it does not already exist.

    Parameters
    ----------
    path : str
        the path to make

    Returns
    -------
    path : str
        the path unchanged
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    try:
        os.makedirs(path)
    except OSError:
        pass
    return path


def build_config_full_path(config, section, relativePathOption,
                           relativePathSection=None,
                           defaultPath=None,
                           baseDirectoryOption='baseDirectory'):
    """
    Get a full path from a base directory and a relative path

    Parameters
    ----------
    config : mpas_tools.config.MpasConfigParser
        configuration from which to read the path

    section : str
        the name of a section in `config`, which must have an option
        ``baseDirectory``

    relativePathOption : str
        the name of an option in ``section`` of the relative path within
        ``baseDirectory`` (or possibly an absolute path)

    relativePathSection : str, optional
        the name of a section for ``relativePathOption`` if not ``section``

    defaultPath : str, optional
        the name of a path to return if the resulting path doesn't exist.

    baseDirectoryOption : str, optional
        the name of the option in ``section`` for the base directorys

    Returns
    -------
    fullPath : str
        The full path to the given relative path within the given
        ``baseDirectory``
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    if relativePathSection is None:
        relativePathSection = section

    subDirectory = config.get(relativePathSection, relativePathOption)
    if os.path.isabs(subDirectory):
        fullPath = subDirectory
    else:
        fullPath = '{}/{}'.format(config.get(section, baseDirectoryOption),
                                  subDirectory)

    if defaultPath is not None and not os.path.exists(fullPath):
        fullPath = defaultPath
    return fullPath


def get_region_mask(config, regionMaskFile):
    """
    Get the full path for a region mask with a given file name

    Parameters
    ----------
    config : mpas_tools.config.MpasConfigParser
        configuration from which to read the path

    regionMaskFile : str
        the file name of the region mask, typically a relative path

    Returns
    -------
    fullFileName : str
        The absolute path to the given fileName within the custom or base
        diagnostics directories
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    if os.path.isabs(regionMaskFile):
        fullFileName = regionMaskFile
    else:
        tryCustom = config.get('diagnostics', 'customDirectory') != 'none'
        found = False
        fullFileName = None
        if tryCustom:
            # first see if region mask file is in the custom directory
            regionMaskDirectory = build_config_full_path(
                config, 'diagnostics', 'regionMaskSubdirectory',
                baseDirectoryOption='customDirectory')

            fullFileName = '{}/{}'.format(regionMaskDirectory,
                                          regionMaskFile)
            found = os.path.exists(fullFileName)

        if not found:
            # no, so second see if mapping files are in the base directory
            regionMaskDirectory = build_config_full_path(
                config, 'diagnostics', 'regionMaskSubdirectory',
                baseDirectoryOption='base_path')

            fullFileName = '{}/{}'.format(regionMaskDirectory,
                                          regionMaskFile)
            found = os.path.exists(fullFileName)

        if not found:
            # still not found, point to a local mask directory
            maskSubdirectory = build_config_full_path(config, 'output',
                                                      'maskSubdirectory')
            make_directories(maskSubdirectory)

            fullFileName = '{}/{}'.format(maskSubdirectory,
                                          regionMaskFile)

    return fullFileName


def build_obs_path(config, component, relativePathOption=None,
                   relativePathSection=None, relativePath=None):
    """

    Parameters
    ----------
    config : mpas_tools.config.MpasConfigParser
        configuration from which to read the path

    component : {'ocean', 'seaIce', 'iceberg'}
        the prefix on the ``*Observations`` section in ``config``, which must
        have an option ``obsSubdirectory``

    relativePathOption : str, optional
        the name of an option in `section` of the relative path within
        ``obsSubdirectory`` (or possibly an absolute path)

    relativePathSection : str, optional
        the name of a section for ``relativePathOption`` if not
        ``<component>Observations``

    relativePath : str, optional
        As an alternative to giving the option (and possibly section) of the
        relative path, it can be supplied directly

    Returns
    -------
    fullPath : str
        The full path to the given relative path within the observations
        directory for the given component
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    obsSection = '{}Observations'.format(component)

    if relativePath is None:
        if relativePathSection is None:
            relativePathSection = obsSection
        relativePath = config.get(relativePathSection, relativePathOption)

    if os.path.isabs(relativePath):
        fullPath = relativePath
    else:
        obsSubdirectory = config.get(obsSection, 'obsSubdirectory')
        if os.path.isabs(obsSubdirectory):
            fullPath = '{}/{}'.format(obsSubdirectory, relativePath)
        else:
            basePath = config.get('diagnostics', 'customDirectory')
            fullPath = '{}/{}/{}'.format(basePath, obsSubdirectory,
                                         relativePath)
            if basePath == 'none' or not os.path.exists(fullPath):
                basePath = config.get('diagnostics', 'base_path')
                fullPath = '{}/{}/{}'.format(basePath, obsSubdirectory,
                                             relativePath)

    return fullPath


def check_path_exists(path):
    """
    Raise an exception if the given path does not exist.

    Parameters
    ----------
    path : str
        Absolute path

    Raises
    ------
    OSError
        If the path does not exist
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    if not (os.path.isdir(path) or os.path.isfile(path)):
        raise OSError('Path {} not found'.format(path))


def get_files_year_month(fileNames, streamsFile, streamName):
    """
    Extract the year and month from file names associated with a stream

    Parameters
    ----------
    fileNames : list of str
        The names of files with a year and month in their names.

    streamsFile : ``StreamsFile``
        The parsed streams file, used to get a template for the

    streamName : str
        The name of the stream with a file-name template for ``fileNames``

    Returns
    -------
    years, months : list of int
        The years and months for each file in ``fileNames``
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    template = streamsFile.read_datetime_template(streamName)
    template = os.path.basename(template)
    dts = [datetime.strptime(os.path.basename(fileName), template) for
           fileName in fileNames]

    years = [dt.year for dt in dts]
    months = [dt.month for dt in dts]

    return years, months


def decode_strings(da):
    """
    Decode to unicode strings an array that might either be char or string type
    in the NetCDF file.

    Parameters
    ----------
    da : ``xarray.DataArray``
        the data array of strings to decode

    Returns
    -------
    strings : list
        The data array as a list of unicode strings
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    if da.dtype.type is numpy.string_:
        strings = [bytes.decode(name) for name in da.values]
    else:
        strings = [name for name in da.values]

    return strings


def copyfile(src, dst):
    """ Copy a file, retrying if temporarily unavailable """

    try:
        shutil.copyfile(src, dst)
    except BlockingIOError:
        # this is an occasional problem on Chrysalis.  Try a slow copy
        warnings.warn('Making a slow copy from {} to {}'.format(src, dst))
        with open(src, 'rb') as fsrc, open(dst, 'wb') as fdst:
            shutil.copyfileobj(fsrc, fdst)
