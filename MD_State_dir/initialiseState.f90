!*****************************************************************************
subroutine initialiseState( state,                                           &
                        a, b, c,                                             &
                        mass, position, temperature, velocity, lattice,      &
                        periodic ) 
!*****************************************************************************
!
!  Project: MolecularDynamics 
!
!  Module: State (PUBLIC routine)
!
!  Created on: Fri 2 Nov 15:55:03 2018  by Eduardo R. Hernandez
!
!  Last Modified: Thu 17 Jan 18:59:00 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This subroutine is called to "fill" an empty but allocated instance of State
! with state information (atom coordinates, cell parameters, etc)
!
! variables:
!
! type ( state_type ), intent (INOUT) :: state
!     this is the state that is going to be allocated and filled with the data
!     passed next.
!
! integer, intent (IN) :: nAtoms
!     the number of atoms in the system; used to dimension several arrays below
!
! real (dp), dimension (3), intent (IN) :: a, b, c
!     the lattice vectors of the simulation cell
!
! real (dp), dimension (:), intent (IN) :: mass
!     an array of dimension (nAtoms) containing the masses of all atoms in 
!     Daltons (1/12 of the mass of C12)
!
! real (dp), dimension (:,:), intent (IN) :: position
!     an array of dimension (3,nAtoms) containing the coordinates of atoms;
!     if the optional flag 'lattice' (see below) is passed as .true., these
!     are assumed to be passed in lattice representation (x,y,z in [0,10])
!
! real (dp), intent(IN), optional :: temperature
!     the temperature of the simulation; the user can provide either the simulation
!     temperature (in which case random initial momenta are generated according to 
!     the Maxwel-Boltzmann distribution at that temperature), or else a set of 
!     initial velocities (see next entry)
!
! real (dp), dimension (:,:), intent (IN), optional :: velocity
!     optionally, the client program can pass initial velocities to be used
!     in the construction of the momenta that are used internally; if absent,
!     random momental will be sampled from the Maxwell-Boltzmann distribution
!     corresponding to the temperature of the simulation. If provided, velocities
!     must be in units compatible with the position; also note that if the flag
!     'lattice' is set to .true., the velocities must also be in this representation,
!     to be compatible with the atomic coordinates (see position, above)
!
! logical, intent (IN), optional :: lattice
!     if .true., this variable indicates that coordinates (position) and velocity are
!     given in lattice representation. If set to .false. (or absent), it is assumed
!     that coordinates (and velocities) are in Cartesian representation in units
!     of the client program.
!
! logical, intent (IN), optional :: periodic
!     if .true. (default) this indicates that periodic boundary conditions will be 
!     imposed; if .false. then a cluster geometry is assumed. The only relevance of this
!     here is that in the latter case the momenta are generated such that not only the
!     total linear momentum is zero, but also the total angular momentum is zero, to 
!     avoid whole-cluster rotation, which makes analysis more cumbersome. 
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  arguments

      type ( state_type ), intent (INOUT) :: state

      real (dp), dimension (3), intent (IN) :: a  ! lattice vectors in 
      real (dp), dimension (3), intent (IN) :: b  ! length units of the 
      real (dp), dimension (3), intent (IN) :: c  ! client program

      real (dp), dimension (:), intent (IN) :: mass
      real (dp), dimension (:,:), intent (IN) :: position

      real (dp), intent (IN), optional :: temperature 
      real (dp), dimension (:,:), intent (IN), optional :: velocity

      logical, intent (IN), optional :: lattice
      logical, intent (IN), optional :: periodic

!*****************************************************************************
!  local variables

      integer :: i

      real (dp), dimension (3,3) :: LatticeVectors

      real (dp), allocatable, dimension (:,:) :: momentum

!*****************************************************************************
!  start subroutine

      if ( debug ) write(*,*) 'Entering initialiseState()'  

      ! set the cell

      LatticeVectors(1:3,1) = lengthConversionFactor * a
      LatticeVectors(1:3,2) = lengthConversionFactor * b
      LatticeVectors(1:3,3) = lengthConversionFactor * c

      call initialiseCell( state, LatticeVectors )

      ! store masses internally

      state % mass = massConversionFactor * mass

      ! now the atomic positions

      if ( present( lattice ) ) then

         if ( lattice ) then  ! coordinates are in lattice representation

            state % position = position
            state % representation % position = .false.

         else

            state % position = lengthConversionFactor * position
            state % representation % position = .true.

         end if

      else

         state % position = lengthConversionFactor * position
         state % representation % position = .true.

      end if 

      if ( present( periodic ) ) then

         PBC = periodic

      else

         PBC = .true.

      end if

      if ( present( temperature ) .and. ( .not. present( velocity ) ) ) then

         call generateMomenta( state, temperature )

      else if ( present( velocity ) ) then 

         if ( present( lattice ) ) then

            if ( lattice ) then  ! velocities are in lattice units

               allocate( momentum( 3, state % nAtoms ) )

               do i = 1, state % nAtoms 
                  momentum(1:3,i) = state % mass(i) * velocity(1:3,i)
               end do

               state % momentum = momentum
               state % representation % momentum = .false.

               deallocate( momentum )

            else 

               do i = 1, state % nAtoms 
                  state % momentum(1:3,i) =                &
                     velConversionFactor * state % mass(i) * velocity(1:3,i)
               end do

               state % representation % momentum = .true.  ! cartesian rep. for momenta

            end if

         else   ! lattice not present, so Cartesian
 
            do i = 1, state % nAtoms 
               state % momentum(1:3,i) =                &
                  state % mass(i) * velocity(1:3,i)
                 !velConversionFactor * state % mass(i) * velocity(1:3,i)
            end do

            state % representation % momentum = .true.  ! cartesian rep. for momenta

         end if

         ! write(*,*) 'velConversionFactor = ', velConversionFactor
         ! write(*,*) 'vel in createState ', velocity(1:3,1)
         ! write(*,*) 'createState ', state % momentum(1:3,1)

      end if  

      if ( debug ) write(*,*) 'Exiting initialiseState'

end subroutine initialiseState

