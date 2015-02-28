# -*- coding: utf-8 -*-
"""
Created on Thu Dec 18 08:09:12 2014

@author: tjk
"""

from pylab import *
from math import atan2
 
def fitSine(tList,yList,freq):
   '''
       freq in Hz
       tList in sec
   returns
       phase in degrees
   '''
   b = matrix(yList).T
   rows = [ [sin(freq*2*pi*t), cos(freq*2*pi*t), 1] for t in tList]
   A = matrix(rows)
   (w,residuals,rank,sing_vals) = lstsq(A,b)
   phase = atan2(w[1,0],w[0,0])*180/pi
   amplitude = norm([w[0,0],w[1,0]],2)
   bias = w[2,0]
   return (phase,amplitude,bias)
 
if __name__=='__main__':
   import random
 
   tList = arange(0.0,1.0,0.001)
   print tList
   tSamples = arange(0.0,1.0,0.01)
   random.seed(0.2)
   phase = -23
   amplitude = 1
   bias = -0.9
   frequency = 5
   yList = amplitude*sin(tList*frequency*2*pi+phase*pi/180.0)+bias
   ySamples = amplitude*sin(tSamples*frequency*2*pi+phase*pi/180.0)+bias
   yMeasured = [y+random.normalvariate(0,.01) for y in ySamples]
   #print yList
   (phaseEst,amplitudeEst,biasEst) = fitSine(tSamples,yMeasured,frequency)
   print ('Phase estimate = %f, Amplitude estimate = %f, Bias estimate = %f'
       % (phaseEst,amplitudeEst,biasEst))
        
   yEst = amplitudeEst*sin(tList*frequency*2*pi+phaseEst*pi/180.0)+biasEst
 
   figure(1)
   plot(tList,yList,'b')
   plot(tSamples,yMeasured,'+r',markersize=12,markeredgewidth=2)
   plot(tList,yEst,'-g')
   xlabel('seconds')
   legend(['True value','Measured values','Estimated value'])
   grid(True)
   show()
