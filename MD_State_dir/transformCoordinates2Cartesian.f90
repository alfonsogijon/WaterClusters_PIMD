!*****************************************************************************
subroutine transformCoordinates2Cartesian( state )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: State (PRIVATE routine)
!
!  Created on: Tue 19 Jun 13:44:20 2018  by Eduardo R. Hernandez
!
!  Last Modified: Thu 17 Jan 19:05:49 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This subroutine transforms the coordinates from Lattice to Cartesian representation, 
! when the client program provides them in that representation.
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  interface variables

      type ( state_type ), intent (INOUT) :: state

!*****************************************************************************
!  internal variables

      integer :: i

      real (dp), dimension (3) :: coordinates, r, s, t

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering transformCoordinates2Cartesian()'

      if ( state % representation % position ) then  ! if coordinates are already cartesian, do nothing

         return

      else    ! coordinates are in lattice rep., so transform them

         r = (/ state % cell % a(1), state % cell % b(1), state % cell % c(1) /)
         s = (/ state % cell % a(2), state % cell % b(2), state % cell % c(2) /)
         t = (/ state % cell % a(3), state % cell % b(3), state % cell % c(3) /)

         do i=1, state % nAtoms

            coordinates = state % position(1:3,i)

            state % position(1,i) = dot_product( r, coordinates )
            state % position(2,i) = dot_product( s, coordinates )
            state % position(3,i) = dot_product( t, coordinates )

         end do

         state % representation % position = .true.

      end if

      if ( debug ) write(*,*) 'Exiting transformCoordinates2Cartesian()'

end subroutine transformCoordinates2Cartesian

