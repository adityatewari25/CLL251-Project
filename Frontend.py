import streamlit as st
import matplotlib.pyplot as plt
import pandas as pd
import scipy
from datetime import datetime
from meteostat import Point, Hourly, Stations
from pyfluids import HumidAir, InputHumidAir
import math

# Constants and functions
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



# Streamlit Interface
st.title('PCM Usage Calculator')

# User Inputs
L0 = st.number_input("Length of the box (m)", value=1.0)
W0 = st.number_input("Width of the box (m)", value=1.0)
H0 = st.number_input("Height of the box (m)", value=1.0)
t0 = st.number_input("Thickness (m)", value=0.1)
keps = st.number_input("Thermal Conductivity of Insulation (W/m.K)", value=0.036)
Ta = st.number_input("Ambient Temperature (°C)", value=40) 
Ti = st.number_input("Internal Temperature (°C)", value=5) 

kpcm = st.number_input("Thermal Conductivity of PCM (W/m.K)", value=0.2)
bpcm = st.number_input("Specific Heat Capacity of PCM (J/kg.K)", value=0.125)
cpcm = st.number_input("Specific Heat of PCM (J/kg)", value=2000)
latentpcm = st.number_input("Latent Heat of PCM (J/kg)", value=200000)
totalt = st.number_input("Total Duration (seconds)", value=48 * 3600)
dpcms = st.number_input("Density of PCM Solid Phase (kg/L)", value=0.88)
dpcml = st.number_input("Density of PCM Liquid Phase (kg/L)", value=0.77)

start = datetime(2018, 6, 18, 23, 59)
end = datetime(2018, 6, 19, 23, 59)
stations = Stations()
stations = stations.nearby(28.704060, 77.102493)
station = stations.fetch(1)
data = Hourly(station, start, end).fetch()['temp'].tolist()

L,W,H,t = L0,W0,H0,t0
# Calculate PCM Usage
total_volume, vol_melt, tot_vol = get_volume()

st.subheader('Results')
st.write(f"Total Volume of PCM required: {total_volume} L")

# Plots
L,W,H = L0/2,W0/2,H0/2
vol = []
vol_pcm = []
for i in range(1, 200,5):
    L *= 1.05
    W *= 1.05
    H *= 1.05
    vol.append(L*W*H)
    vol_pcm.append(get_volume()[0])

fig1, ax1 = plt.subplots()
ax1.plot(vol,vol_pcm)
ax1.set_xlabel("Volume of Box (m3)")
ax1.set_ylabel("PCM Required (L)")
ax1.set_title("PCM Required vs. Volume of Box")
st.pyplot(fig1)


L,W,H = L0,W0,H0
fig2, ax2 = plt.subplots()
ax2.plot(range(len(tot_vol)), tot_vol)
ax2.set_xlabel("Time Duration (hr)")
ax2.set_ylabel("PCM Required (L)")
ax2.set_title("PCM Required vs. Time Duration")
st.pyplot(fig2)

L,W,H = L0,W0,H0
t = t0/2
t_l = []
vol_pcm = []
for i in range(1, 200,5):
    t *= 1.05
    t_l.append(t)
    vol_pcm.append(get_volume()[0])
fig3, ax3 = plt.subplots()
ax3.plot(t_l, vol_pcm)
ax3.set_xlabel("Thickness (m)")
ax3.set_ylabel("PCM Required (L)")
ax3.set_title("PCM Required vs. Thickness")
st.pyplot(fig3)
