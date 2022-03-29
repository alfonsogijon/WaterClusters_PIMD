!*****************************************************************************
subroutine getMass( state_chain, mass )
!*****************************************************************************
!
!  Project: MolecularDynamics
! 
!  Module: State (PUBLIC routine)
!
!  Created on: Dec 2020 by Alfonso Gijón
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  arguments

      type ( StateChain_type ), intent (INOUT) :: state_chain  ! the state to be read

      double precision, dimension(:), intent(OUT) :: mass

!*****************************************************************************
!  local variables


!*****************************************************************************
!  begin subroutine 

      if ( debug ) write(*,*) 'Entering getMass()'


      mass = state_chain % mass / massConversionFactor
      
      if ( debug ) write(*,*) 'Exiting getMass()'

end subroutine getMass
