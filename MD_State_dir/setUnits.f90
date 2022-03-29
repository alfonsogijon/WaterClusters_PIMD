!*****************************************************************************
subroutine setUnits( lengthUnit, energyUnit, pressureUnit,                   &
                     lengthName, energyName, pressureName )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: Structure (PUBLIC routine)
!
!  Created on: Thu 14 Jun 11:22:27 2018  by Eduardo R. Hernandez
!
!  Last Modified: Fri 18 Oct 12:28:27 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
! 
! This subroutine is called once at the begining; its purpose is to 
! set a number of unit conversion factors in order to manage communication
! to/from the client program. Internally, MD_State and MolecularDynamics
! work in atomic units (i.e. distance in Bohr, energy in Hartree, mass
! in units of the electron mass; all remaining units are built from these);
! the client program using MD_State and MolecularDynamics will generally work
! in arbitrary units, so the user has to 
! instruct MD_State about the units in which it is going to pass/get
! information. This is done by setting the arguments of this routine to
! the appropriate numerical conversion factor between the user's units
! and the corresponding atomic units. This must be done for length and energy;
! it is assumed that masses will be provided in Daltons (1/12 of the mass of 
! carbon). 
!
! Defining a unit conversion for input defines it also for output (if the user
! provides energies in eV, say, MD_State will also return them in eV). Note that
! if stress is provided by the client program (as required in a constant pressure
! simulation) it must come in units of energy (not pressure).
!
! Pressure is an exception. The units of pressure that are most frequently used in
! condensed matter tent to be kBar or GPa, and not the units consistent with the
! energy/length^3 in use by the user. So, optionally, the user can specify the
! output units for the internal pressure.
!
! WARNING!!!  IN CASE OF CONSTANT-PRESSURE SIMULATIONS:
!
! In the case of constant-pressure simulations, the client program will need to
! provide the stress (as well as the forces on atoms); the stress is the derivative
! of the potential energy with respect to the strain tensor components; since the 
! latter are dimensionless, the dimensions of the stress tensor components are 
! those of an energy; hence the client program must provide the stress in the same
! units as it provides energies (and not, as is sometimes done, in units of pressure!)
!
!!!!!!!!!!!!!!! MANDATORY ARGUMENTS !!!!!!!!!!!!!!!
!
! length (real): specifies client length unit; recognized standard options are
!                Bohr (atomic unit, internal unit), Angstrom 
!                   
! energy (real): specifies client energy unit; recognized standard options are
!               Hartree (atomic unit, internal unit), Rydberg, eV, kcalmol, kJmol 
!
!!!!!!!!!!!!!!! OPTIONAL ARGUMENTS !!!!!!!!!!!!!!!
!
! pressure (real): if this is a constant-pressure simulation, or if the
!                 user is going to ask for values of the instantaneous
!                 internal pressure, specify here the units pressure is to be
!                 output in kBar or GPa
!
! lengthName (char) : an optional name for the used length unit
! energyName (char) : an optional name for the used energy unit
! pressureName (char):an optional name for the used pressure unit
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
   implicit none

!*****************************************************************************
!  variables

      real (dp), intent (IN) :: lengthUnit
      real (dp), intent (IN) :: energyUnit
      real (dp), intent (IN), optional :: pressureUnit

      character (len = 20 ), intent (IN), optional :: lengthName
      character (len = 20 ), intent (IN), optional :: energyName
      character (len = 20 ), intent (IN), optional :: pressureName

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering setUnits()'

      ! set the length unit conversion factor

      lengthConversionFactor = lengthUnit

      ! see if we can identify the length units set by the user

      if ( abs( lengthUnit - Bohr ) .lt. 1.0e-4_dp ) then
         !write(*,*) 'Length units are Bohr (atomic units)'
      else if ( abs( lengthUnit - Angstrom ) .lt. 1.0e-4_dp ) then
         !write(*,*) 'Length units are Angstrom'
      else 
         write(*,*) 'Unidentified length units ', lengthName
      end if

      ! set the energy unit conversion factor

      energyConversionFactor = energyUnit

      if ( abs( energyUnit - Hartree ) .lt. 1.0e-4_dp ) then
         !write(*,*) 'Energy units are Hartree (atomic units)'
      else if ( abs( energyUnit - Rydberg ) .lt. 1.0e-4_dp ) then
         !write(*,*) 'Energy units are Rydberg'
      else if ( abs( energyUnit - eV ) .lt. 1.0e-4_dp ) then
         !write(*,*) 'Energy units are eV'
      else if ( abs( energyUnit - kcalmol ) .lt. 1.0e-4_dp ) then
         !write(*,*) 'Energy units are kcal/mol'
      else if ( abs( energyUnit - kJoulemol ) .lt. 1.0e-4_dp ) then
         !write(*,*) 'Energy units are kJ/mol'
      else 
         write(*,*) 'Unidentified energy units ', energyName
      end if

      ! with length and energy units, we can set the conversion factors
      ! for other magnitudes

      forceConversionFactor = energyConversionFactor / lengthConversionFactor
      massConversionFactor = 1822.8885302342505_dp ! convert Dalton to electron masses
      timeConversionFactor = sqrt( au_to_joules / amu ) / adu
      velConversionFactor = lengthConversionFactor / timeConversionFactor  ! is this generally correct?
      volumeConversionFactor = lengthConversionFactor * lengthConversionFactor * lengthConversionFactor

!  finally, pressure units, if required

      if ( present( pressureUnit ) ) then

         pressureConversionFactor = pressureUnit

         if ( abs( pressureUnit - kBar ) .lt. tiny ) then
            write(*,*) 'Pressure units are kBar'
         else if ( abs( pressureUnit - GPa ) .lt. tiny ) then
            write(*,*) 'Pressure units are GPa'
         else 
            write(*,*) 'Unidentified pressure units ', pressureName
         end if

      else ! use default as GPa

         pressureConversionFactor = GPa

      end if

      if ( debug ) write(*,*) 'Exiting setUnits()'

end subroutine setUnits

