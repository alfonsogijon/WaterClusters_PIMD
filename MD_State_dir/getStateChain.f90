!*****************************************************************************
subroutine getStateChain( state_chain, &
     iReplica, state, position, momentum, force, nReplicas, nAtoms, &
     nDegreesOfFreedom, mass, nThermostats, thermostat, energy, H0, representation, &
     qPotentialEnergy, qKineticEnergy)
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
! This subroutine extracts from an object of type StateChain_type specific values
! as requested by the caller. The first variable, which is mandatory, is the
! state chain object being interrogated. Then follows a keyword-list of desired
! values, al of which are optional, that are returned upon request.
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  variables

      type ( StateChain_type ), intent (INOUT) :: state_chain

      integer, optional, intent(IN) :: iReplica

      type ( state_type ), optional, intent(OUT) :: state

      real(dp), optional, dimension(:,:), intent (OUT) :: position

      real(dp), optional, dimension(:,:), intent (OUT) :: momentum

      real(dp), optional, dimension(:,:), intent (OUT) :: force      

      integer, optional, intent(OUT) :: nReplicas

      integer, optional, intent(OUT) :: nAtoms

      integer, optional, intent(OUT) :: nDegreesOfFreedom

      real(dp), optional, dimension(:), intent(OUT) :: mass

      integer, intent (OUT), optional :: nThermostats

      type ( thermostat_type ), dimension (:), intent (OUT), optional :: thermostat

      type ( energy_type ), intent (OUT), optional :: energy

      real (dp), intent (OUT), optional :: H0

      type ( representation_type ), intent (IN), optional :: representation

      real (dp), intent(OUT), optional :: qPotentialEnergy

      real (dp), intent(OUT), optional :: qKineticEnergy
      
!*****************************************************************************
!  variables

      integer i

      integer nat
      real(dp), dimension(3) :: a,b,c
      

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering getStateChain()'

      if (present(iReplica) .and. present(state) ) then
      
         state = state_chain % replica(iReplica)

      end if

      if (present(iReplica) .and. present(position) ) then

         !call getConfiguration( state_chain % replica(iReplica), nat, a, b, c, position )

         position = state_chain % replica(iReplica) % position

      end if

      if (present(iReplica) .and. present(momentum) ) then

         !call inquire( state_chain % replica(iReplica), momentum = momentum )

         momentum = state_chain % replica(iReplica) % momentum

      end if

      if (present(iReplica) .and. present(force) ) then

         force = state_chain % replica(iReplica) % force

      end if      

      if (present(nReplicas) ) then
      
         nReplicas = state_chain % nReplicas

      end if

      if (present(nAtoms) ) then
      
         nAtoms = state_chain % nAtoms

      end if

      if ( present(nDegreesOfFreedom) ) then
      
         nDegreesOfFreedom = state_chain % nDegreesOfFreedom

      end if

      if (present(mass) ) then
      
         mass = state_chain % mass

      end if

      if ( present( nThermostats ) ) then

         nThermostats = state_chain % nThermostats

      end if

      if ( present( thermostat ) ) then

         do i = 1, state_chain % nThermostats
 
            thermostat(i) = state_chain % thermostat(i) 

         end do

      end if

      if ( present( energy ) ) then

         energy = state_chain % energy

      end if

      if ( present( H0 ) ) then

         H0 = state_chain % H0 % E

      end if

      if ( present( qPotentialEnergy ) ) then

         qPotentialEnergy = state_chain % qPotentialEnergy

      end if

      if ( present( qKineticEnergy ) ) then

         qKineticEnergy = state_chain % qKineticEnergy

      end if

      if ( debug ) write(*,*) 'Exiting getStateChain()'
      
end subroutine getStateChain

