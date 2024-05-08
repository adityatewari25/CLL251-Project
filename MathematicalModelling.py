import math
import scipy
from pyfluids import HumidAir, InputHumidAir
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
from meteostat import Point, Daily, Hourly, Stations
import matplotlib
import matplotlib.pyplot as plt

e = math.e

L= 1
W = 1
H = 1
t = 0.1

keps = 0.036

Ta = 40 + 273.15
Ti = 5 + 273.15

def get_T_Q(L,W,H,t,keps,Ta,Ti):

    def kair(T): 
        return (HumidAir().with_state(
        InputHumidAir.altitude(300),
        InputHumidAir.temperature(T-273.15),
        InputHumidAir.relative_humidity(50)
    ).conductivity)

    def Pr(T):
        return (HumidAir().with_state(
        InputHumidAir.altitude(300),
        InputHumidAir.temperature(T-273.15),
        InputHumidAir.relative_humidity(50)
    ).prandtl)

    def gamma(T):
        return (HumidAir().with_state(
            InputHumidAir.altitude(300),
            InputHumidAir.temperature(T-273.15),
            InputHumidAir.relative_humidity(50),
        ).kinematic_viscosity)

    ri = ((L*W +W*H + H*L)/(2*e))**(0.5)
    ro = ri + t

    def Q(To,Ti):
        return (4*math.pi*ro*ri*keps*(To-Ti))/(ro-ri)


    def ha(Nu,kair):
        return Nu*kair/(2*ro)

    def Nu(Ra,Pr):
        answer = ( 0.589*(Ra**(1/4)) )
        answer = answer/ ( (1+(0.469/Pr)**(9/16))**(4/9) )
        return 2+answer

    def Ra(To,Ta,Pr,gamma):
        return -9.82*(1/(To+Ta))*(To-Ta)*((2*ro)**3)*Pr/(gamma**2)

    def Q2(To):
        return ha(Nu(Ra(To,Ta,Pr(Ta),gamma(Ta)),Pr(Ta)),kair(Ta))*4*math.pi*(ro**2)*(Ta-To)

    def balance(To):
        return Q(To,Ti)-Q2(To)

    res = scipy.optimize.fsolve(balance,100)[0]-273.15

    return res, Q(res+273.15,Ti)

L= 1
W = 1
H = 1
t = 0.1
keps = 0.036
Ta = 40
Ti = 5


kpcm=0.2
bpcm=0.125
cpcm=2000
latentpcm=200000
totalt=48*3600
dpcms=0.88
dpcml=0.77
vart=0

def mass(Ta,Ti,latentpcm,totalt):
    m=list(get_T_Q(L,W,H,t,keps,Ta+273.15,Ti+273.15))[1]*totalt/(latentpcm)
    return m
def volume(dpcml,Ta,totalt):
    vol=mass(Ta,Ti,latentpcm,totalt)/dpcml
    return vol

def get_volume():
    start = datetime(2018, 6, 18, 23, 59)
    end = datetime(2018, 6, 19, 23, 59)

    # Get hourly data
    stations = Stations()
    stations = stations.nearby(28.704060, 77.102493)
    station = stations.fetch(1)
    data = Hourly(station, start, end)
    data = data.fetch()['temp'].tolist()

    # Print DataFrame
    vol_melt = []
    tot_vol = []
    vol = 0
    for i in data:
        res = volume(dpcml,i,3600)
        vol += res
        vol_melt.append(res)
        tot_vol.append(vol)

    return vol, vol_melt, tot_vol

