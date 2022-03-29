!*****************************************************************************
subroutine inquire( state, temperature, pressure, stress,                    &
                    kineticEnergy, potentialEnergy, thermostatEnergy,        &
                    barostatEnergy, constantOfMotion, Cell,                  &
                    position, momentum, thermostat )
!*****************************************************************************
!
!  Project: MolecularDynamics 
!
!  Module: State (PUBLIC routine)
!
!  Created on: Fri 27 Apr 12:15:38 2018  by Eduardo R. Hernandez
!
!  Last Modified: Tue 16 Apr 18:49:45 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This subroutine inquires state to obtain properties of the current system
! configuration, such as the temperature, pressure, or different energy 
! components. Properties are output in the chosen units of the client program.
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  arguments

      type ( state_type ), intent (INOUT) :: state

      real (dp), intent (OUT), optional :: temperature
      real (dp), intent (OUT), optional :: pressure ! includes kinetic contribution
      real (dp), intent (OUT), dimension (3,3), optional :: stress ! includes kinetic contribution
      real (dp), intent (OUT), optional :: kineticEnergy
      real (dp), intent (OUT), optional :: potentialEnergy
      real (dp), intent (OUT), optional :: thermostatEnergy
      real (dp), intent (OUT), optional :: barostatEnergy
      real (dp), intent (OUT), optional :: constantOfMotion
      real (dp), intent (OUT), optional, dimension (:,:) :: position
      real (dp), intent (OUT), optional, dimension (:,:) :: momentum

      type (lattice_type), intent (OUT), optional :: Cell
      type (thermostat_type), dimension (:), intent (OUT), optional :: thermostat

!*****************************************************************************
!  local variables

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering inquire()'

!  make sure that the kinetic energy is consistent with current momenta

!     call calculateKineticEnergy( state )

      if ( present( temperature ) ) then
 
         temperature = two * state % energy % kinetic / ( state % nDegreesOfFreedom * boltzmann_k )

      end if

      if ( present ( pressure ) .and. state % stress % calculated ) then

         call calculateInternalPressure( state )

         pressure = state % totalInternalStress % pressure / pressureConversionFactor

         ! the total internal stress is only calculated if previously the internal pressure
         ! has been calculated too; otherwise internal stress would be undefined, as both
         ! are calculated in calculateInternalPressure

         if ( present( stress ) ) then

            stress = state % totalInternalStress % cartesian / pressureConversionFactor

         end if

      end if

      if ( present( kineticEnergy ) ) then

         kineticEnergy = state % energy % kinetic / energyConversionFactor

      end if

      if ( present( potentialEnergy ) ) then

         potentialEnergy = state % energy % potential / energyConversionFactor

      end if

      if ( present( thermostatEnergy ) ) then

         thermostatEnergy = state % energy % thermostat / energyConversionFactor

      end if

      if ( present( barostatEnergy ) ) then

         barostatEnergy = state % energy % barostat / energyConversionFactor

      end if

      if ( present( constantOfMotion ) ) then

         constantOfMotion = state % energy % constantOfMotion / energyConversionFactor

      end if   

      if ( present( position ) ) then

         position = state % position / lengthConversionFactor

      end if

      if ( present( momentum ) ) then

         momentum = state % momentum 

      end if

      if ( present( Cell ) ) then

         Cell = state % Cell

         ! convert to client length units

         Cell % a = Cell % a / lengthConversionFactor
         Cell % b = Cell % b / lengthConversionFactor
         Cell % c = Cell % c / lengthConversionFactor
         Cell % a_modulus =  Cell % a_modulus / lengthConversionFactor
         Cell % b_modulus =  Cell % b_modulus / lengthConversionFactor
         Cell % c_modulus =  Cell % c_modulus / lengthConversionFactor
         Cell % volume = Cell % volume / volumeConversionFactor
         Cell % metricTensor = Cell % metricTensor /   &
                 ( lengthConversionFactor * lengthConversionFactor )
 
      end if

      if ( present( thermostat ) ) then

         thermostat = state % thermostat

      end if

      if ( debug ) write(*,*) 'Exiting inquire()'

end subroutine inquire

