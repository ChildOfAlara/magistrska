# -*- coding: utf-8 -*-
import os
import math
from lmfit import minimize, Parameters, Parameter, report_fit, fit_report
from scipy import *
from sympy import *
import scipy.optimize
import numpy as np
from numpy import array
import matplotlib.pyplot as plt
import plotly.plotly as py

path = "Magistrska"
superslovar = {}
HP_TFE_percentage = []
CCDA_TFE_percentage = []
P1_TFE_percentage = []
PAA_TFE_percentage = []
P3_TFE_percentage = []
HP_Gw0 = []
CCDA_Gw0 = []
P1_Gw0 = []
PAA_Gw0 = []
P3_Gw0 = []
HP_Hw0 = []
CCDA_Hw0 = []
P1_Hw0 = []
PAA_Hw0 = []
P3_Hw0 = []

for filename in os.listdir(path):
    f = open('{}/{}'.format(path, filename), 'r')

    file_readlines = f.readlines()

    ime = file_readlines[0]
    ime = ime[9:]

    T_line = file_readlines[20]
    T = T_line[14:16]
    T = float(T)

    lines = file_readlines[25:]

    lines = [line.strip('\n') for line in lines]

    if filename[0:2] == "HP":
        AA = 22
        M = 2487.7
        HP_TFE_percentage.append(filename[3:(len(filename)-2)])

    if filename[0:4] == "CCDA": #ccda37-72
        AA = 36
        M = 4369.82
        CCDA_TFE_percentage.append(filename[5:(len(filename)-2)])

    if filename[0:2] == "P1": #Phd
        AA = 22
        M = 2498.7
        P1_TFE_percentage.append(filename[3:(len(filename)-2)])

    if filename[0:3] == "PAA": #35-63
        AA = 29
        M = 3565.0
        PAA_TFE_percentage.append(filename[4:(len(filename)-2)])

    if filename[0:2] == "P3": #ccda37-62
        AA = 26
        M = 3191.6
        P3_TFE_percentage.append(filename[3:(len(filename)-2)])

    data_signal = {}
    data_helix = {}
    T_sequence = []
    Fhexp_sequence = []

    for line in lines:
        if "_data_end_" in line:
            break
        if ("wavelength 222.00 nm" in line) or ("_avg_time_" in line) or ("_data_" in line):
            continue
        x = line.split()[0]
        y = float(line.split()[1])
        data_signal[int(float(x))] = y
        theta_coil = 2000 + (-50)*(float(x))
        theta_helix = -40000 + (-150)*(float(x))
        elipticnost = y #elipicnost iz signala
        molarna_elipticnost = elipticnost*(M/AA)/(0.1*1)  #molarna elipticnost iz elipticnosti, c=0,1mg/ml;l=1mm
        Fhexp = (molarna_elipticnost - theta_coil)/(theta_helix - theta_coil)
        data_helix[int(float(x))] = Fhexp
        T_sequence.append(float(x))
        Fhexp_sequence.append(Fhexp)

    T_sequence=array(T_sequence)
    Fhexp_sequence=array(Fhexp_sequence)



    w = Symbol("w")
    #AA = Symbol("AA")
    Nr = AA-2
    v = 0.048
    Z = 1 + Nr*v + (Nr*(Nr-1)*v**2)/2 + v**2*w*((Nr-2)-(Nr-1)*w + w**(Nr-1))/((1-w)**2)
    #Z = 1 + (AA-2)*0.048 + ((AA-2)*((AA-2)-1)*0.048**2)/2 + 0.048**2*w*(((AA-2)-2)-((AA-2)-1)*w + w**((AA-2)-1))/((1-w)**2)
    #Fhteo = (0.002304*w*(-AA + 3 + w**(AA - 3)*(AA - 3)/w)/(-w + 1)**2 + 0.004608*w*(AA - w*(AA - 3) + w**(AA - 3) - 4)/(-w + 1)**3 + 0.002304*(AA - w*(AA - 3) + w**(AA - 3) - 4)/(-w + 1)**2)*w/(1 + (AA-2)*0.048 + ((AA-2)*((AA-2)-1)*0.048**2)/2 + 0.048**2*w*(((AA-2)-2)-((AA-2)-1)*w + w**((AA-2)-1))/((1-w)**2))
    Zprime = Z.diff(w)
    Fhteor = Zprime*w/Z
    DeltaF = Fhteor - Fhexp
    DeltaFprime = DeltaF.diff(w)
    odvod = lambdify(w, Zprime)
    Z = lambdify(w, Z)

    T0 = 298.15
    T_sequence_K = T_sequence + 273.15
    R = 8.314

    params = Parameters()
    params.add('Gw0', value= 0.0, vary=True)
    params.add('Hw0', value= 0.0, vary=True)

    print T_sequence
    print Fhexp_sequence

    def f_DeltaF(params, T_sequence_K, Fhexp_sequence):
        #Gw = Gwo*T/T0 - Hw0(1-T/T0)
        #w = np.exp(-(Gw0*T_sequence_K/T0 - Hw0*(1-T_sequence_K/T0))/(R*T_sequence_K))
        Gw0 = params['Gw0'].value
        Hw0 = params['Hw0'].value
        w = np.exp(-(Gw0*T_sequence_K/T0 - Hw0*(1-T_sequence_K/T0))/(R*T_sequence_K))
        Fhteor = odvod(w)*(np.exp(-(Gw0*T_sequence_K/T0 - Hw0*(1-T_sequence_K/T0))/(R*T_sequence_K)))/Z(w)
        return Fhteor - Fhexp_sequence

    result = minimize(f_DeltaF, params, args=(T_sequence_K, Fhexp_sequence), method="nelder")
    print result.params
    #print a.__dict__
    print result.params['Gw0']._val

    #print params.valuesdict  bound method wtf?
    #print params.valuesdict()  returna začetno vresnost: 0.0
    print result.success
    print fit_report(result)

    if filename[0:2] == "HP":
        HP_Gw0.append(result.params['Gw0']._val)
        HP_Hw0.append(result.params['Hw0']._val)

    if filename[0:4] == "CCDA":
        CCDA_Gw0.append(result.params['Gw0']._val)
        CCDA_Hw0.append(result.params['Hw0']._val)

    if filename[0:2] == "P1":
        P1_Gw0.append(result.params['Gw0']._val)
        P1_Hw0.append(result.params['Hw0']._val)

    if filename[0:3] == "PAA":
        PAA_Gw0.append(result.params['Gw0']._val)
        PAA_Hw0.append(result.params['Hw0']._val)

    if filename[0:2] == "P3":
        P3_Gw0.append(result.params['Gw0']._val)
        P3_Hw0.append(result.params['Hw0']._val)

    print ime
    print data_helix
    print "4°C: " + str(data_helix[4])
    print "26°C: " + str(data_helix[26])
    print "96°C: " + str(data_helix[96])
    print ""

#print HP_Gw0
#print HP_Hw0
#print CCDA_Gw0
#print CCDA_Hw0
#print P1_Gw0
#print P1_Hw0
#print PAA_Gw0
#print PAA_Hw0
#print P3_Gw0
#print P3_Hw0

plt.plot(HP_TFE_percentage, HP_Gw0, 'ro')
plt.set_ylim(ymin=0)
plt.show()

print max([max(HP_Gw0), max(CCDA_Gw0), max(P1_Gw0), max(PAA_Gw0), max(P3_Gw0)])
print min([min(HP_Gw0), min(CCDA_Gw0), min(P1_Gw0), min(PAA_Gw0), min(P3_Gw0)])
print max([max(HP_Hw0), max(CCDA_Hw0), max(P1_Hw0), max(PAA_Hw0), max(P3_Hw0)])
print min([min(HP_Hw0), min(CCDA_Hw0), min(P1_Hw0), min(PAA_Hw0), min(P3_Hw0)])


#sghdhf
