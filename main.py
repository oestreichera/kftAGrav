# Python Script accopanying: (Insert paper URL/DOI)
#
# Python version: 3.8.3
# Date: 23.03.2022
# Author: Alexander Oestreicher
#         Institute for Theoretical Physics, Heidelberg University, Germany
# E-mail: oestreicher@thphys.uni-heidelberg.de

# Import libraries and seetting up plot styles.
import matplotlib.pyplot as plt
import numpy as np

plt.rc('axes', titlesize=16)          # fontsize of the axes title
plt.rc('axes', labelsize=16)          # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)         # fontsize of the tick labels
plt.rc('ytick', labelsize=14)         # fontsize of the tick labels
plt.rc('legend', fontsize=14)         # legend fontsize
plt.rc('figure', titlesize=16)        # fontsize of the figure title
plt.style.use('tableau-colorblind10') # clorblind palette for plotting


# Load GR reference spectrum and lcdm expansion function and
# gravitational constant.
k,P_GR=np.loadtxt('gr_spectrum.d', unpack = True)
a,G_GR,E_GR=np.loadtxt('gr_cosmo.d', unpack = True)
a_min=a[0]

# Load coefficients for the Taylor expansion
# x,k loaded only from the first file, since it is the same in all files.
x = np.genfromtxt('dP_dG.d', usecols = (0), unpack=True)
k = np.genfromtxt('dP_dG.d', usecols = (1), unpack=True)

c = len(x) #number of columns
b = len(k) #number of rows

dP_dG=np.loadtxt('dP_dG.d', usecols = np.arange(2,130))
dP_dE=np.loadtxt('dP_dE.d', usecols = np.arange(2,130))

A_G=np.loadtxt('A_G.d', usecols = np.arange(2,130))
A_E=np.loadtxt('A_E.d', usecols = np.arange(2,130))

# Function calculating the alternative power spectrum:
#
# Input: arrays containing the scale factor at which the functions are
# tabulated, the gravitational constant and expansion function for the desired
# alternative gravity model.
#
# Integration is done with the trapez-method.
#
# Functions are expected to be tabulated from a=0.001 to a=1.0 in 128 table
# steps logarithmically spaced in order to match the tables containing the
# taylor coefficients. If this is not the case an error is thrown.
#
# The gravitational constant is expected to be normalized to the reference
# GR value.


def P_AG (a_AG, G_AG, E_AG):

    if (np.any(abs(a_AG/a-1)>0.001)):
        raise ValueError ("Input scale factors do not match tabulated values!")

    # Calculate the change in the expansion function and grav. constant and
    # normalize them according to the convention used in KFT.

    deltaE = (E_AG-E_GR)/E_GR[0]
    deltaG = (G_AG-G_GR)

    # Integrate coefficients with changes in the expansion function
    # and gravitational coupling over the scale factor,
    # factor 1/a_min included to respect KFT scale factor normalization.

    dP_dG_int = np.zeros(b)
    j=0
    while j<c:
    	dP_dG_int[j] = np.trapz(dP_dG[:,j]*deltaG[:]/a_min,x[:])
    	j += 1

    dP_dE_int = np.zeros(b)
    j=0
    while j<c:
    	dP_dE_int[j] = np.trapz(dP_dE[:,j]*deltaE[:]/a_min,x[:])
    	j += 1

    A_G_int = np.zeros(b)
    j=0
    while j<c:
    	A_G_int[j] = np.trapz(A_G[:,j]*deltaG[:]/a_min,x[:])
    	j += 1

    A_E_int = np.zeros(b)
    j=0
    while j<c:
    	A_E_int[j] = np.trapz(A_E[:,j]*deltaE[:]/a_min,x[:])
    	j += 1

    # Return the alternative power spectrum as given by the Taylor expansion.
    # Result is an array containing the spectrum values as a function of the
    # wave number k.

    return P_GR + dP_dG_int + dP_dE_int + A_G_int + A_E_int

# Plot the results for the two different proca models shwon in the paper,
# with a range of parameters q_v.

qv=np.array(['2000','1500','1000','0500','0100','0050','0010','0001'])
qv_label=np.array([2.0,1.5,1.0,0.5,0.1,0.05,0.01,0.001])


plt.figure()
# Read the scale factor, expansion function, and gravitational constant
# for proca model 1 with different qv from the data files, pass them to
# P_AG and plot the relative change in the power spectrum.
#
# Functions are expected to be tabulated from a=0.001 to a=1.0 in 128 table
# steps logarithmically spaced in order to match the tables containing the
# taylor coefficients. If this is not the case an error is thrown.
#
# The gravitational constant is expected to be normalized to the reference
# GR value.
#
# To adapt the script to a different alternative theory of gravity import
# your file here instead.
#
i=0
while i<len(qv):
    a_AG,G_AG,E_AG=np.loadtxt('Model1/proca'+qv[i]+'.d', unpack=True)
    plt.plot(k,(P_AG(a_AG,G_AG,E_AG)/P_GR-1),label='$q_v=$'+str(qv_label[i]))
    i+=1
plt.xscale('log')
plt.xlabel('wavenumber $k$ $[h\mathrm{Mpc}^{-1}]$')
plt.ylabel('$P_{\mathrm{Proca}}(k)/P_{\mathrm{GR}}(k)-1$')
plt.title('Model 1')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('Proca1.pdf')

plt.figure()
i=0
while i<len(qv):
    a_AG,G_AG,E_AG=np.loadtxt('Model2/proca'+qv[i]+'.d', unpack=True)
    plt.plot(k,(P_AG(a_AG,G_AG,E_AG)/P_GR-1),label='$q_v=$'+str(qv_label[i]))
    i+=1
plt.xscale('log')
plt.xlabel('wavenumber $k$ $[h\mathrm{Mpc}^{-1}]$')
plt.ylabel('$P_{\mathrm{Proca}}(k)/P_{\mathrm{GR}}(k)-1$')
plt.title('Model 2')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('Proca2.pdf')
plt.show()
