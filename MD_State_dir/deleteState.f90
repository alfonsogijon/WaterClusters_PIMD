!*****************************************************************************
subroutine deleteState( state )
!*****************************************************************************
!
!  Project: MolecularDynamics 
!
!  Module: State (PUBLIC routine)
!
!  Created on: Thu 19 Jul 11:28:36 2018  by Eduardo R. Hernandez
!
!  Last Modified: Thu 17 Jan 18:57:58 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! type ( state_type ), intent (INOUT) :: state
!     this is the state that is going to be deallocated.
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  arguments

      type ( state_type ), intent (INOUT) :: state

!*****************************************************************************
!  local variables

!*****************************************************************************
!  start subroutine

      if ( debug ) write(*,*) 'Entering deleteState()'  

      deallocate( state % position )
      deallocate( state % force )
      deallocate( state % momentum )
      deallocate( state % mass )
      if ( state % nThermostats .gt. 0 ) deallocate( state % thermostat )

      if ( debug ) write(*,*) 'Exiting deleteState'

end subroutine deleteState

