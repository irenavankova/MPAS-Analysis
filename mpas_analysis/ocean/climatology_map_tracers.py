# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/main/LICENSE
"""
Analysis tasks for comparing Antarctic climatology maps against observations
and reanalysis data.
"""
# Authors
# -------
# Xylar Asay-Davis

import numpy

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask

from mpas_analysis.ocean.remap_depth_slices_subtask import \
    RemapDepthSlicesSubtask
from mpas_analysis.shared.plot import PlotClimatologyMapSubtask
from mpas_analysis.ocean.remap_sose_climatology import RemapSoseClimatology

from mpas_analysis.shared.io.utility import build_obs_path


class ClimatologyMapTracers(AnalysisTask):
    """
    An analysis task for comparison of antarctic field against the Southern
    Ocean State Estimate
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, mpasClimatologyTask,
                 controlConfig=None):
        """
        Construct the analysis task.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        controlconfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        fields = \
            [{'prefix': 'temperature',
              'mpas': 'timeMonthly_avg_activeTracers_temperature',
              'units': r'$\degree$C',
              'titleName': 'Potential Temperature',
              '3D': True,
              'obsFilePrefix': 'pot_temp',
              'obsFieldName': 'theta',
              'obsBotFieldName': 'botTheta'},
             {'prefix': 'clifw',
              'mpas': 'timeMonthly_avg_freshwaterTracers_landIceFreshwaterConcentration',
              'units': r'',
              'titleName': 'LIFW Concentration',
              '3D': True,
              'obsFilePrefix': 'clifw',
              'obsFieldName': 'clifw',
              'obsBotFieldName': 'botClifw'}]

        tags = ['climatology', 'horizontalMap', 'sose', 'publicObs',
                'antarctic'] + [field['prefix'] for field in fields]

        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapTracers, self).__init__(
            config=config, taskName='climatologyMapTracers',
            componentName='ocean',
            tags=tags)

        sectionName = self.taskName

        fileSuffix = config.get(sectionName, 'fileSuffix')
        if fileSuffix.endswith('.nc'):
            fileSuffix = fileSuffix.strip('.nc')

        fieldList = config.getexpression(sectionName, 'fieldList')
        fields = [field for field in fields if field['prefix'] in fieldList]

        # read in what seasons we want to plot
        seasons = config.getexpression(sectionName, 'seasons')

        if len(seasons) == 0:
            raise ValueError('config section {} does not contain valid '
                             'list of seasons'.format(sectionName))

        comparisonGridNames = config.getexpression(sectionName,
                                                   'comparisonGrids')

        if len(comparisonGridNames) == 0:
            raise ValueError('config section {} does not contain valid '
                             'list of comparison grids'.format(
                                 sectionName))

        if not numpy.any([field['3D'] for field in fields]):
            depths = None
        else:
            depths = config.getexpression(sectionName, 'depths')

            if len(depths) == 0:
                raise ValueError('config section {} does not contain valid '
                                 'list of depths'.format(sectionName))

        variableList = [field['mpas'] for field in fields
                        if field['mpas'] != 'velMag']

        shallowVsDeepColormapDepth = config.getfloat(
            sectionName, 'shallowVsDeepColormapDepth')

        shallow = []
        for depth in depths:
            if depth == 'top':
                shallow.append(True)
            elif depth == 'bot':
                shallow.append(False)
            else:
                shallow.append(depth >= shallowVsDeepColormapDepth)

        if depths is None:
            remapMpasSubtask = RemapMpasClimatologySubtask(
                mpasClimatologyTask=mpasClimatologyTask,
                parentTask=self,
                climatologyName='SOSE',
                variableList=variableList,
                seasons=seasons,
                comparisonGridNames=comparisonGridNames,
                iselValues=None)
        else:
            remapMpasSubtask = RemapMpasVelMagClimatology(
                mpasClimatologyTask=mpasClimatologyTask,
                parentTask=self,
                climatologyName='SOSE',
                variableList=variableList,
                seasons=seasons,
                depths=depths,
                comparisonGridNames=comparisonGridNames,
                iselValues=None)

        for field in fields:
            fieldPrefix = field['prefix']
            upperFieldPrefix = fieldPrefix[0].upper() + fieldPrefix[1:]
            sectionName = '{}{}'.format(self.taskName, upperFieldPrefix)

            if field['3D']:
                fieldDepths = depths
            else:
                fieldDepths = None

            if controlConfig is None:
                remapObsSubtask = None
                galleryName = None
                refTitleLabel = 'None'

                #refFieldName = None
                #outFileLabel = 'None'
                diffTitleLabel = 'None'
                refFieldName = field['mpas']
                outFileLabel = '{}'.format(fieldPrefix)
            else:
                remapObsSubtask = None
                controlRunName = controlConfig.get('runs', 'mainRunName')
                galleryName = 'Control: {}'.format(controlRunName)
                refTitleLabel = galleryName

                refFieldName = field['mpas']
                outFileLabel = '{}SOSE'.format(fieldPrefix)
                diffTitleLabel = 'Main - Control'

            if field['3D']:
                fieldDepths = depths
            else:
                fieldDepths = [None]

            for comparisonGridName in comparisonGridNames:
                for depthIndex, depth in enumerate(fieldDepths):
                    for season in seasons:

                        subtaskName = 'plot{}_{}_{}'.format(upperFieldPrefix,
                                                            season,
                                                            comparisonGridName)
                        if depth is not None:
                            subtaskName = '{}_depth_{}'.format(subtaskName,
                                                               depth)

                        subtask = PlotClimatologyMapSubtask(
                            parentTask=self,
                            season=season,
                            comparisonGridName=comparisonGridName,
                            remapMpasClimatologySubtask=remapMpasSubtask,
                            remapObsClimatologySubtask=remapObsSubtask,
                            controlConfig=controlConfig,
                            depth=depth,
                            subtaskName=subtaskName)

                        configSectionName = 'climatologyMapTracers{}'.format(
                                upperFieldPrefix)

                        # if available, use a separate color map for shallow
                        # and deep
                        if depth is not None:
                            if shallow[depthIndex]:
                                suffix = 'Shallow'
                            else:
                                suffix = 'Deep'
                            testSectionName = '{}{}'.format(configSectionName,
                                                            suffix)
                            if config.has_section(testSectionName):
                                configSectionName = testSectionName

                        subtask.set_plot_info(
                            outFileLabel=outFileLabel,
                            fieldNameInTitle=field['titleName'],
                            mpasFieldName=field['mpas'],
                            refFieldName=refFieldName,
                            refTitleLabel=refTitleLabel,
                            diffTitleLabel=diffTitleLabel,
                            unitsLabel=field['units'],
                            imageCaption=field['titleName'],
                            galleryGroup=field['titleName'],
                            groupSubtitle=None,
                            groupLink='{}Tracers'.format(fieldPrefix),
                            galleryName=galleryName,
                            configSectionName=configSectionName)

                        self.add_subtask(subtask)


class RemapMpasVelMagClimatology(RemapDepthSlicesSubtask):
    """
    A subtask for computing climatologies of velocity magnitude from zonal
    and meridional components
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def customize_masked_climatology(self, climatology, season):
        """
        Construct velocity magnitude as part of the climatology

        Parameters
        ----------
        climatology : ``xarray.Dataset`` object
            the climatology data set

        season : str
            The name of the season to be masked

        Returns
        -------
        climatology : ``xarray.Dataset`` object
            the modified climatology data set
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the base class's version of this function so we extract
        # the desired slices.
        climatology = super(RemapMpasVelMagClimatology,
                            self).customize_masked_climatology(climatology,
                                                               season)

        if 'timeMonthly_avg_velocityZonal' in climatology and \
                'timeMonthly_avg_velocityMeridional' in climatology:
            zonalVel = climatology.timeMonthly_avg_velocityZonal
            meridVel = climatology.timeMonthly_avg_velocityMeridional
            climatology['velMag'] = numpy.sqrt(zonalVel**2 + meridVel**2)
            climatology.velMag.attrs['units'] = 'm s$^{-1}$'
            climatology.velMag.attrs['description'] = 'velocity magnitude'

        return climatology
