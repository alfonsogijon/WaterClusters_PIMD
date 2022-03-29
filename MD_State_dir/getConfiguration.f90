!*****************************************************************************
subroutine getConfiguration( state, nAtoms, a, b, c, position, lattice )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: State (PUBLIC routine)
!
!  Created on: Wed 13 Jun 18:04:27 2018  by Eduardo R. Hernandez
!
!  Last Modified: Tue 22 Jan 11:22:43 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This subroutine is designed as a means of communicating with the client 
! program; the latter will need to be told the lattice vectors and atom positions
! in units that it understands. Lattice vectors and ionic positions are 
! privately held by state, and are thus directly not accessible to the 
! client program. This subroutine provides, given an instance of state, the
! lattice vectors a, b, c and the ionic positions, in the chosen length 
! units of the client program. Optionally, atomic positions can be returned
! in lattice coordinates, by appropriately setting the optional flag Lattice 
! as .true.
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  variables

      type ( state_type ), intent (INOUT) :: state

      integer , intent (OUT) :: nAtoms

      real (dp), dimension (3), intent (OUT) :: a
      real (dp), dimension (3), intent (OUT) :: b  ! the three lattice vectors
      real (dp), dimension (3), intent (OUT) :: c

      real (dp), dimension ( 3, state % nAtoms ), intent (OUT) :: position

      logical, optional, intent (IN) :: lattice

!*****************************************************************************
!  variables

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering getConfiguration()'

      ! convert lattice vectors to client program length units

      a = state % Cell % a / lengthConversionFactor
      b = state % Cell % b / lengthConversionFactor
      c = state % Cell % c / lengthConversionFactor

      ! the number of atoms

      nAtoms = state % nAtoms

      ! now atomic coordinates

      if ( present( lattice ) ) then  ! check if client wants coordinates in lattice rep.

         if ( lattice ) then

            call transformCoordinates2Lattice( state )  ! make sure they are in lattice rep.

            position = state % position

         else 

            call transformCoordinates2Cartesian( state ) ! make sure they are in cartesian rep.

            position = state % position / lengthConversionFactor

         end if

      else    ! the default case is that coordinates are returned in Cartesian rep.

         call transformCoordinates2Cartesian( state ) ! make sure they are in cartesian rep.

         position = state % position / lengthConversionFactor

      end if

      if ( debug ) write(*,*) 'Exiting getConfiguration()'

end subroutine getConfiguration

