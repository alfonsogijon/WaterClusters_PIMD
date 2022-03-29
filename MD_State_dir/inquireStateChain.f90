!*****************************************************************************
subroutine inquireStateChain( state_chain, temperature, &
     kineticEnergy, potentialEnergy, thermostatEnergy,        &
     barostatEnergy, constantOfMotion, &
     qKineticEnergy, qPotentialEnergy, qKineticEnergy_v2 )
!*****************************************************************************
!
!  Project: MolecularDynamics 
!
!  Module: State (PUBLIC routine)
!
!  Created on: March 2019 by Alfonso Gijón
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
!  argument

      type ( StateChain_type ), intent (INOUT) :: state_chain

      real (dp), intent (OUT), optional :: temperature
      real (dp), intent (OUT), optional :: kineticEnergy
      real (dp), intent (OUT), optional :: potentialEnergy
      real (dp), intent (OUT), optional :: thermostatEnergy
      real (dp), intent (OUT), optional :: barostatEnergy
      real (dp), intent (OUT), optional :: constantOfMotion
      real (dp), intent (OUT), optional :: qKineticEnergy
      real (dp), intent (OUT), optional :: qPotentialEnergy
      real (dp), intent (OUT), optional :: qKineticEnergy_v2

!*****************************************************************************
!  variables

      integer :: i

      real(dp) :: const
      real(dp) :: factor

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering inquireStateChain()'

!  make sure that the kinetic energy is consistent with current momenta

      call calculateKEChain( state_chain )

      if ( present( temperature ) ) then
 
         temperature = two * state_chain % energy % kinetic / &
              ( state_chain % nDegreesOfFreedom * boltzmann_k )

      end if

      if ( present( kineticEnergy ) ) then

         kineticEnergy = state_chain % energy % kinetic / energyConversionFactor

      end if

      if ( present( potentialEnergy ) ) then

         potentialEnergy = state_chain % energy % potential / energyConversionFactor

      end if

      if ( present( thermostatEnergy ) ) then

         thermostatEnergy = state_chain % energy % thermostat / energyConversionFactor

      end if

      if ( present( barostatEnergy ) ) then

         barostatEnergy = state_chain % energy % barostat / energyConversionFactor

      end if

      if ( present( constantOfMotion ) ) then

         constantOfMotion = state_chain % energy % constantOfMotion / energyConversionFactor

      end if

      if ( present( qKineticEnergy ) ) then

         qKineticEnergy = state_chain % qKineticEnergy / energyConversionFactor

      end if

      if ( present( qPotentialEnergy ) ) then

         qPotentialEnergy = state_chain % qPotentialEnergy / energyConversionFactor

      end if

      if ( present( qKineticEnergy_v2 ) ) then

         qKineticEnergy_v2 = state_chain % qKineticEnergy_v2 / energyConversionFactor

      end if

      if ( debug ) write(*,*) 'Exiting inquireStateChain()'

end subroutine inquireStateChain
