!*****************************************************************************
subroutine transformCoordinates2Lattice( state )
!*****************************************************************************
!
!  Project: MolecularDynamics 
!
!  Module: State (PRIVATE routine)
!
!  Created on: Thu  8 Mar 11:39:29 2018  by Eduardo R. Hernandez
!
!  Last Modified: Thu 17 Jan 19:05:12 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This subroutine transforms the coordinates from Cartesian to Lattice
! representation. Internally, type State contains the coordinates in Cartesian
! representation, but they may be required by MD routines or client program
! in Lattice representation. This soubroutine takes care of the transformation.
! The internal representation is not changed, only the copied array is.
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

   real (dp), dimension (3) :: coordinates

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering transformCoordinates2Lattice()'

      if ( state % representation % position ) then ! coordinates are in cartesian rep., so work to do

         do i=1, state % nAtoms
         
            coordinates = state % position(1:3,i)

            state % position(1,i) = dot_product( state % reciprocalCell % a, &
                                              coordinates )
            state % position(2,i) = dot_product( state % reciprocalCell % b, &
                                              coordinates )
            state % position(3,i) = dot_product( state % reciprocalCell % c, &
                                              coordinates )

         end do

         state % representation % position = .false.  ! now they are in lattice rep.

      else   ! coordinates are in lattice already, so there is nothing to do

         return

      end if

      if ( debug ) write(*,*) 'Exiting transformCoordinates2Lattice()'

end subroutine transformCoordinates2Lattice

