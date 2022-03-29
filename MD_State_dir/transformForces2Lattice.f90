!*****************************************************************************
subroutine transformForces2Lattice( state )
!*****************************************************************************
!
!  Project: MolecularDynamics 
!
!  Module: State (PRIVATE routine)
!
!  Created on: Thu  8 Mar 12:10:34 2018  by Eduardo R. Hernandez
!
!  Last Modified: Thu 17 Jan 19:06:17 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This subroutine transforms the forces from Cartesian representation (the 
! way they are stored internally) to Lattice, when they are requested in this
! representation by the program calling get() (usually this call will be from
! a subroutine in MolecularDynamics module).
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

   real(dp), dimension (3) :: vector

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering transformForces2Lattice()'

      if ( state % representation % force ) then    ! forces are in cartesian rep., so transform

         do i=1, state % nAtoms
         
            vector = state % force(1:3,i)

            state % force(1,i) = dot_product( state % cell % a, vector )
            state % force(2,i) = dot_product( state % cell % b, vector )
            state % force(3,i) = dot_product( state % cell % c, vector )

         end do

         state % representation % force = .false.  ! now they are in lattice rep.

      else        ! forces already in lattice rep., so do nothing

         return

      end if

      if ( debug ) write(*,*) 'Exiting transformForces2Lattice()'

end subroutine transformForces2Lattice

