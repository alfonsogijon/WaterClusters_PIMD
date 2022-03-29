!*****************************************************************************
subroutine setStateChain( state_chain, iReplica, state, &
     position, momentum, force, &
     energy, thermostat, potentialEnergy, qPotentialEnergy, qKineticEnergy, &
     qKineticEnergy_v2 )
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
! This subroutine has the purpose of setting the values of variables contained
! in an instance of state_chain to new numerical values. An object of type tate_chain
! must be passed as first argument; the remaining arguments are optional and
! identified by keyword, so that a single variable, a subset or the whole set
! can be specified in one go. The updated state_chain object is then returned.
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  variables

      type ( StateChain_type ), intent (INOUT) :: state_chain

      integer, intent(IN), optional :: iReplica

      type( state_type ), intent(IN), optional :: state

      real(dp), intent(IN), dimension(:,:), optional :: position
      
      real(dp), intent(IN), dimension(:,:), optional :: momentum
      
      real(dp), intent(IN), dimension(:,:), optional :: force

      type( energy_type ), intent(IN), optional :: energy

      type( thermostat_type ), dimension(:), intent(IN), optional :: thermostat

      real (dp), intent(IN), optional :: potentialEnergy

      real (dp), intent(IN), optional :: qPotentialEnergy

      real (dp), intent(IN), optional :: qKineticEnergy

      real (dp), intent(IN), optional :: qKineticEnergy_v2
      
!*****************************************************************************
!  variables

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering setStateChain()'

      if ( present(iReplica) .and. present(state) ) then

         state_chain % replica(iReplica) = state
         
      end if

      if ( present(iReplica) .and. present(position) ) then

         state_chain % replica(iReplica) % position = position
         
      end if

      if ( present(iReplica) .and. present(momentum) ) then

         state_chain % replica(iReplica) % momentum = momentum
         
      end if      
      
      if ( present(iReplica) .and. present(force) ) then

         state_chain % replica(iReplica) % force = force

      end if

      if ( present(energy) ) then

         state_chain % energy = energy

         if ( .not. state_chain % H0 % set ) then

            state_chain % H0 % E = state_chain % energy % kinetic + &
                 state_chain % energy % potential + &
                 state_chain % energy % barostat + &
                 state_chain % energy % thermostat

            state_chain % H0 % set = .true.

         end if         
         
      end if

      if ( present(thermostat) ) then

         state_chain % thermostat = thermostat
         
      end if

      if ( present(potentialEnergy) ) then
      
         state_chain % energy % potential = potentialEnergy

      end if

      if ( present(qPotentialEnergy) ) then
      
         state_chain % qPotentialEnergy = qPotentialEnergy

      end if

      if ( present(qKineticEnergy) ) then

         state_chain % qKineticEnergy = qKineticEnergy
         
      end if

      if ( present(qKineticEnergy_v2) ) then

         state_chain % qKineticEnergy_v2 = qKineticEnergy_v2
         
      end if

      if ( debug ) write(*,*) 'Exiting setStateChain()'

end subroutine setStateChain

