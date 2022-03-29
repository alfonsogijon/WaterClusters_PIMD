!*****************************************************************************
module MD_constants_module
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Created on: Wed  7 Mar 18:50:23 2018  by Eduardo R. Hernandez
!
!  Last Modified: Fri  8 Mar 12:13:20 2019
!
!*****************************************************************************
!  modules used

   use MD_precision_module

!*****************************************************************************
!  No implicits please!
   
   implicit none

!*****************************************************************************
!  variables

! numerical constants

   real (dp), parameter :: zero = 0.0_dp
   real (dp), parameter :: one = 1.0_dp
   real (dp), parameter :: two = 2.0_dp
   real (dp), parameter :: three = 3.0_dp
   real (dp), parameter :: four = 4.0_dp
   real (dp), parameter :: five = 5.0_dp
   real (dp), parameter :: six = 6.0_dp
   real (dp), parameter :: seven = 7.0_dp
   real (dp), parameter :: eight = 8.0_dp
   real (dp), parameter :: nine = 9.0_dp
   real (dp), parameter :: ten = 10.0_dp
   real (dp), parameter :: eleven = 11.0_dp
   real (dp), parameter :: twelve = 12.0_dp
   real (dp), parameter :: half = 0.5_dp
   real (dp), parameter :: sixth = 1.0_dp / 6.0_dp
   
   real (dp), parameter :: verysmall = 1.0e-7_dp
   real (dp), parameter :: tiny = 1.0e-16_dp

! pi to 50 digits
   real (dp), parameter :: pi = 3.14159265358979323846264338327950288419716939937510_dp
   real (dp), parameter :: deg = pi / 180.0_dp

!   real(dp), parameter :: zeps = 1.0d-10
!   real(dp), parameter :: eps = 1.0d-10

! unit conversion

   real (dp), parameter :: Bohr = 1.0_dp
   real (dp), parameter :: Angstrom = Bohr / 0.529177_dp
   real (dp), parameter :: Hartree = 1.0_dp
   real (dp), parameter :: Rydberg = Hartree / two
   
   real (dp), parameter :: Ang    = one / 0.529177_dp
   real (dp), parameter :: eV     = one / 13.60580_dp / two
   real (dp), parameter :: kBar   = one / 1.47108e5_dp
   real (dp), parameter :: GPa    = kBar * ten
   real (dp), parameter :: Kelvin = eV / 11604.45_dp
   real (dp), parameter :: Debye  = 0.393430_dp
   real (dp), parameter :: kcalmol = Hartree / 627.509474_dp
   real (dp), parameter :: kJoulemol = Hartree / 2625.49963865      

   real (dp), parameter :: rad_to_deg = 180.0_dp / pi
   real (dp), parameter :: au_to_joules = 4.3598e-18_dp
   real (dp), parameter :: GPa_to_au = 3.398923e-5_dp
   real (dp), parameter :: s_to_atu = 2.418884326e-17_dp

! physical constants

   real (dp), parameter :: adu = 0.52917715e-10_dp
   real (dp), parameter :: amu = 9.109535e-31_dp
   real (dp), parameter :: boltzmann_k = 3.1666590210560118e-06_dp  ! in Hartree K^-1
 
!  The easy way to make sense of units conversion:

!  real(dp), parameter :: Bohr   = 1.0_dp
!  real(dp), parameter :: Rydberg = 1.0_dp
!  real(dp), parameter :: Femtosecond = 1.0_dp
!
!  Ang = Bohr / 0.529177
!   eV = Rydberg / 13.60580
!  Joule = eV / 1.6e-19_dp
!  Meter = Ang / 1.0e-10_dp
!  Pascal = Joule/Meter**2
!   kBar  = Pascal * 1.0e4
!   .... and so on.

end module MD_constants_module
