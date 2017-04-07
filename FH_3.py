# -*- coding: utf-8 -*-
import os
import math

from scipy import *
from sympy import *
from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
import plotly.plotly as py

path = "Magistrska"
superslovar = {}

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

    if filename[0:4] == "CCDA": #ccda37-72
        AA = 36
        M = 4369.82

    if filename[0:2] == "P1": #Phd
        AA = 22
        M = 2498.7

    if filename[0:3] == "PAA": #35-63
        AA = 29
        M = 3565.0

    if filename[0:2] == "P3": #ccda37-62
        AA = 26
        M = 3191.6

    data_signal = {}
    data_helix = {}

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


        w = Symbol("w")
        Nr= AA-2
        v=0.048
        Z = 1 + Nr*v + (Nr*(Nr-1)*v**2)/2 + v**2*w*((Nr-2)-(Nr-1)*w + w**(Nr-1))/((1-w)**2)
        Zprime = Z.diff(w)
        odvod = lambdify(w, Zprime)
        w = 0.000001
        Z = 1 + Nr*v + (Nr*(Nr-1)*v**2)/2 + v**2*w*((Nr-2)-(Nr-1)*w + w**(Nr-1))/((1-w)**2)

        Fhteor = odvod(w)*w*Z
        DeltaF=Fhteor-Fhexp
        print DeltaF

        def f(w):
            y = Zprime
            return y

    x = fsolve(f,0.5)
    print x

    w = float
    funkcija = Zprime/Z - Fhexp

    fsolve(w**2, 0.5)
    w = np.linspace(0,1)
    plt.plot(w,f(w))
    plt.plot(w,np.zeros(len(w)))

    print ime
    print Zprime
    print Fhteor
    print data_helix
    print "4°C: " + str(data_helix[4])
    print "26°C: " + str(data_helix[26])
    print "96°C: " + str(data_helix[96])
    print ""
