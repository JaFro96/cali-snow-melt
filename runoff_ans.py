""" This script consists of a snow melt model and a calibration on two parameters. """
from pcraster import *
from pcraster.framework import *
import matplotlib.pyplot as plt

class MyFirstModel(DynamicModel):
    """ This class executes the field-based model """
    def __init__(self, myParameter, cali_on_melt_rate):
        """ Initialize PCRaster class instance and pass the parameter """
        DynamicModel.__init__(self)
        setclone('dem.map')
        # Parameters to calibrate
        if(cali_on_melt_rate):
            self.meltRate = myParameter
            # original value
            self.temperatureLapseRate = 0.005
            # changed to 0.0043 after calibrating first on temperature lapse rate
            # self.temperatureLapseRate = 0.0043

        else:
            self.temperatureLapseRate = myParameter
            # original value
            self.meltRate = 0.01
            # changed to 0.16 after calibrating first on melt rate
            # self.meltRate = 0.16
    
    def initial(self):
        """ Initial definitions of inputs and parameters"""
        # the measurement location
        self.observationLocation = self.readmap('out')
        dem = self.readmap('dem')
        elevationMeteoStation = 208.1
        elevationAboveMeteoStation = dem - elevationMeteoStation
        self.temperatureCorrection = elevationAboveMeteoStation * self.temperatureLapseRate
        self.report(self.temperatureCorrection, 'tempCor')        
        self.snow = 0.0
        self.ldd = lddcreate(dem, 1e31, 1e31, 1e31, 1e31)
        self.report(self.ldd, 'ldd')
        
        # example how to calculate total precipitation
        self.totPrecip = scalar(0)
        
        # timeseries of modelled outflow at specific location
        self.modelled = TimeoutputTimeseries('modelled', self)
        
        # timeseries of observed - modelled outflow at specific location
        self.diffOut = TimeoutputTimeseries('diff_outflow', self)
        
        # initialize streamflow timeseries as numpy array for directly writing to disk
        self.simulation = numpy.zeros(self.nrTimeSteps())
        
    def dynamic(self):
        """ Transition functions and reports """
        precipitation = timeinputscalar('precip.tss', 1)
        # self.report(precipitation, 'pFromTss')
        temperatureObserved = timeinputscalar('temp.tss', 1)
        # self.report(temperatureObserved, 'tempObs')
        temperature = temperatureObserved - self.temperatureCorrection
        # self.report(temperature, 'temp')
        freezing = temperature < 0.0
        # self.report(freezing, 'freez')
        snowFall = ifthenelse(freezing, precipitation, 0.0)
        # self.report(snowFall, 'snowFall')
        rainFall = ifthenelse(pcrnot(freezing), precipitation, 0.0)
        # self.report(rainFall, 'rain')        
        self.snow = self.snow + snowFall

        potentialMelt = ifthenelse(pcrnot(freezing), temperature*self.meltRate, 0)
        # self.report(potentialMelt, 'pmelt')

        actualMelt = min(self.snow, potentialMelt)
        # self.report(actualMelt, 'amelt')
        
        self.snow = self.snow - actualMelt
        self.report(self.snow, 'snow')
        
        runoffGenerated = actualMelt + rainFall
        self.report(runoffGenerated, 'rg')
        
        discharge = accuflux(self.ldd, runoffGenerated*cellarea())
        self.report(discharge, 'q')
        
        # reading values from a map at the observation location
        runoffAtOutflowPoint = getCellValueAtBooleanLocation(self.observationLocation, discharge)
        
        self.report(scalar(runoffAtOutflowPoint), 'roatop')
        # print('modelled runoff at observation location: ', runoffAtOutflowPoint)
        # add modelled outflow to tss file
        self.modelled.sample(scalar(runoffAtOutflowPoint))
        
        # timeseries of observed outflow at specific location
        observed = timeinputscalar('observed.tss', 1)
        
        # timeseries of difference outflow
        self.outflow_diff = observed - scalar(runoffAtOutflowPoint) 
        self.diffOut.sample(self.outflow_diff)
        
        self.simulation[self.currentTimeStep() - 1] = self.outflow_diff

# Specifies if calibration is executed on melt rate, otherwise it is executed on time lapse
calibrate_melt_rate = True

# prepare a nested list to plot the values        
result = []
# prepare a numpy array with headers to get an overview on the total values
np_result = ["Melt rate parameter","Lapse rate parameter","Mean discharge difference"]

# Standard calibration is on temperature lapse rate
# coarse steps
my_range = numpy.arange(0.0041, 0.006, 0.0001)
my_label = "temperature lapse rate"
my_values = 1

# Calibrate on melt rate:
if(calibrate_melt_rate):
    # coarse steps
    my_range = numpy.arange(0.001,0.02,0.001)
    # possible finer steps:
    # my_range = numpy.arange(0.009,0.01,0.0001) 
    my_label = "melt rate"
    my_values = 0
    
# loop to run the model to calibrate the melt rate
# i specifies the current value of the parameter range
for i in my_range:
    # specify number of time steps
    NR_OF_TIME_STEPS = 181
    # run the model and pass the parameters 
    myModel = MyFirstModel(i, calibrate_melt_rate)
    dynamicModel = DynamicFramework(myModel, NR_OF_TIME_STEPS)
    dynamicModel.run()
    
    # calculate objective function: mean of absolute differences
    runOffModelled = myModel.simulation
    diff = numpy.mean(numpy.abs(runOffModelled))
    print('Melt rate parameter:', numpy.round(myModel.meltRate,4) , 
          'Temperature Lapse rate parameter:', numpy.round(myModel.temperatureLapseRate,4),
          'mean difference: ', diff, ' m^3 water')

    # add parameter values and discharge differences to the nested list
    result.append([numpy.round(myModel.meltRate,4),
                   numpy.round(myModel.temperatureLapseRate,4),
                   numpy.rint(diff)]) 
    # add parameter values and discharge differences to the table
    np_result = numpy.vstack((np_result, [numpy.round(myModel.meltRate,4),
                   numpy.round(myModel.temperatureLapseRate,4),
                   numpy.rint(diff)]))
    
# plot parameter values against outcome of goal function
plt.xlabel(my_label)
plt.ylabel("mean discharge difference in m^3 per day")
plt.title("Summary")
for i in result:
    plt.plot(i[my_values],i[2], 'bo')
plt.show()
