!*****************************************************************************
subroutine setForce( state, potential_energy, force, stress, lattice )
!*****************************************************************************
!
!  Project: MolecularDynamics 
!
!  Module: State (PUBLIC routine)
!
!  Created on: Thu 14 Jun 09:33:11 2018  by Eduardo R. Hernandez
!
!  Last Modified: Wed 23 Jan 11:05:18 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This subroutine receives the forces, and optionally the stress tensor, as
! calculated by the client program, and stores this information in State.
! It is assumed that force and stress information are input in Cartesian
! representation. The client program will normally pass the forces in 
! Cartesian representation, but if the optional flag lattice is set tot .true.,
! then it is understood that forces are passed in lattice representation.
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  variables

      type ( state_type ), intent (INOUT) :: state

      real (dp), intent (IN) :: potential_energy 
      real (dp), dimension ( 3, state % nAtoms ), intent (IN) :: force
      real (dp), dimension ( 3, 3 ), intent (IN), optional :: stress

      logical, intent (IN), optional :: lattice

!*****************************************************************************
!  local variables

      real (dp), dimension ( 3, 3 ) :: matrix_a
      real (dp), dimension ( 3, 3 ) :: matrix_b

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering setForce()'

      state % energy % potential = energyConversionFactor * potential_energy 

      if ( present( lattice ) ) then

         if ( lattice ) then      ! forces are given in lattice rep;  
                                  ! TODO: make sure conversion factor 
                                  ! is right in this case
            state % representation % force = .false.                     

         else

            state % representation % force = .true. ! cartesian

         end if

      else

         state % representation % force = .true.

      end if

      state % force = forceConversionFactor * force

      ! if stress has been passed, convert it too

      if ( present( stress ) ) then

         state % stress % calculated = .true.

         state % stress % cartesian = energyConversionFactor * stress

         matrix_a(1:3,1) = state % reciprocalCell % a
         matrix_a(1:3,2) = state % reciprocalCell % b
         matrix_a(1:3,3) = state % reciprocalCell % c

         matrix_b = matmul( state % stress % cartesian, matrix_a )

         state % stress % lattice = half * matmul( transpose( matrix_a ), matrix_b )

      else

         state % stress % calculated = .false.
         state % stress % cartesian = zero
         state % stress % lattice = zero

      end if

      if ( debug ) write(*,*) 'Exiting setForce()'

end subroutine setForce

