!*****************************************************************************
subroutine deleteStateChain( state_chain )
!*****************************************************************************
!
!  Project: MolecularDynamics 
!
!  Module: State
!
!  Created on: Jan 2019 by Alfonso Gijón
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! type ( StateChain_type ), intent (INOUT) :: state_chain
!     chain of states that is going to be deallocated.
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  arguments

      type ( StateChain_type ), intent (INOUT) :: state_chain

!*****************************************************************************
!  local variables

      integer i

!*****************************************************************************
!  start subroutine
      
      if ( debug ) write(*,*) 'Entering deleteStateChain()'

      deallocate( state_chain % mass )

      if (state_chain % nThermostats .gt. 0) deallocate( state_chain % thermostat )

      do i = 1, state_chain % nReplicas
         call deleteState( state_chain % replica(i) )
      end do

      deallocate( state_chain % replica )

      if ( debug ) write(*,*) 'Exiting deleteStateChain'

end subroutine deleteStateChain

