!*****************************************************************************
subroutine PolymerDistribution( state_chain, polymer )
  
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: State (PUBLIC routine)
!
!  Created on: May 2019 by Alfonso Gijón
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This subroutine compute the spatial distribution of each atom over the polymer
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  variables

      type ( StateChain_type ), intent (IN) :: state_chain

      real(dp), dimension(:,:,:), intent (OUT) :: polymer

!*****************************************************************************
!  variables

      integer i

!*****************************************************************************
!  start of subroutine

      ! read polymer positions
     
      if ( debug ) write(*,*) 'Entering PolymerDistribution()'

      do i = 1, state_chain % nReplicas
      
         polymer(:,:,i) = state_chain % replica(i) % position
         
      end do

      if ( debug ) write(*,*) 'Exiting PolymerDistribution()'
      
end subroutine PolymerDistribution
