!*****************************************************************************
subroutine transformForces2Cartesian( state )
!*****************************************************************************
!
!  Project: MolecularDynamics 
!
!  Module: State (PRIVATE routine)
!
!  Created on: Thu  8 Mar 12:10:34 2018  by Eduardo R. Hernandez
!
!  Last Modified: Mon 10 Jun 12:36:35 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This is a subroutine of module configuration_module; it transforms the 
! forces from Lattice representation, if the client provides them in that way
! to Cartesian, the way in which they are stored by State type.
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

   real(dp), dimension (3) :: force, r, s, t

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering transformForces2Cartesian()'

      if ( state % representation % force ) then      ! forces already Cartesian; do nothing

         return

      else                                            ! forces in lattice rep., so transform

         r = (/ state % reciprocalCell % a(1), state % reciprocalCell % a(2), state % reciprocalCell % a(3) /)
         s = (/ state % reciprocalCell % b(1), state % reciprocalCell % b(2), state % reciprocalCell % b(3) /)
         t = (/ state % reciprocalCell % c(1), state % reciprocalCell % c(2), state % reciprocalCell % c(3) /)

         do i=1, state % nAtoms

            force = state % force(1:3,i)

            state % force(1,i) = dot_product( r, force(1:3) )
            state % force(2,i) = dot_product( s, force(1:3) )
            state % force(3,i) = dot_product( t, force(1:3) )

          end do

          state % representation % force = .true.

       end if

       if ( debug ) write(*,*) 'Exiting transformForces2Cartesian()'

end subroutine transformForces2Cartesian

