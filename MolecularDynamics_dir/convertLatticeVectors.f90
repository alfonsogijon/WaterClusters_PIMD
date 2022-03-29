!*****************************************************************************
subroutine convertLatticeVectors( a, b, c )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Created on: Fri  1 Jun 18:07:01 2018  by Eduardo R. Hernandez
!
!  Last Modified: Fri  1 Jun 18:17:20 2018
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This is only used in constant-pressure simulations; in this case, new 
! cell vectors are obtained by the algorithm, and these must be provided
! in the appropriate client units of length to the client program for it
! to calculate new energy, forces and stress.
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  arguments

      real(dp), intent (OUT), dimension (3) :: a
      real(dp), intent (OUT), dimension (3) :: b ! the output lattice vectors
      real(dp), intent (OUT), dimension (3) :: c

!*****************************************************************************
!  local variables

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering convertLatticeVectors()'

      a = state % cell % a / lengthConversionFactor
      b = state % cell % b / lengthConversionFactor
      c = state % cell % c / lengthConversionFactor

      if ( debug ) write(*,*) 'Exiting convertLatticeVectors()'

end subroutine convertLatticeVectors
