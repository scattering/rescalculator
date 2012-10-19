from numpy import *

def PrefDemo(H,K,L,W,myrescal,p):
    """This is an example of a  "slowly varying prefactor and background" function for use with ConvRes.m.
    In the particular case it is a form factor squared for Ni2+.
    Two modes are assumed, polarized along the b and c axis of the crystal.
    The function also calculates the two magnetic polarization factors.
    It is meant to be used in combination with the cross section function SqwDemo.m
    All arguments are vectors, so dont forget to use ".*" instead of "*", etc."""

    I=p[5] #							% Intensity prefactor.
    bgr=p[6] #						% A flat background.

    # Here we will need to access to direct
    # and reciprocal space lattice parameters. The function GetLattice
    # can be used to get those directly from EXP.
    #[sample,rsample]=GetLattice(EXP);


    #% First, calculate the scattering vector modulus squared
    q2=myrescal.lattice_calculator.modvec(H,K,L,'latticestar')**2

    #% Now, use the Jane Brown approximate expression for Ni2+
    sd=q2/(16*pi**2)
    ff=0.0163*exp(-35.883*sd)+0.3916*exp(-13.223*sd)+0.6052*exp(-4.339*sd)-0.0133

    # Calculate the polarization factors for transverse excitations
    # For this we utilize the function angle2.m
    # Angle between (0,1,0) [direct space] and [h k l]:
    alphay=myrescal.lattice_calculator.angle2(0,1,0,H,K,L)
    # Angle between (0,0,1) [direct space] and [h k l]:
    alphaz=myrescal.lattice_calculator.angle2(0,0,1,H,K,L)
    # Angle between (1,0,0) [direct space] and [h k l]:
    alphax=myrescal.lattice_calculator.angle2(1,0,0,H,K,L)

    # Polarization factors for each of the three modes.
    polx=sin(alphax)**2
    poly=sin(alphay)**2
    polz=sin(alphaz)**2
    prefactor=zeros((3,H.shape[0]),'float64')
    prefactor[0,:]=ff**2*polx*p[5] # %Apply the overall scaling factor here
    prefactor[1,:]=ff**2*poly*p[5]
    prefactor[2,:]=ff**2*polz*p[5]

    #The background will be simply a constant
    bgr=ones(H.shape,'float64')*p[6]

    return prefactor,bgr