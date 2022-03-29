!*****************************************************************************
subroutine readRestartChain( state_chain, restartFile, nSteps )
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
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  arguments

      type ( StateChain_type ), intent (INOUT) :: state_chain  ! the state to be read

      character ( len = 30 ), intent (IN) :: restartFile

      integer, intent (INOUT), optional :: nSteps

!*****************************************************************************
!  local variables

      integer :: i, irep
      integer :: n

      real (dp), dimension (3,3) :: LatticeVectors
      real (dp), allocatable, dimension (:,:) :: position
      real (dp), allocatable, dimension (:,:) :: momentum
      real (dp), allocatable, dimension (:) :: mass

      type ( representation_type ) :: representation
      logical periodic

!*****************************************************************************
!  begin subroutine 

      if ( debug ) write(*,*) 'Entering readRestartChain()'

      open( unit = 9, file = restartFile )

! read the number of steps integrated in previous runs
! this is so that the continuation run can start counting steps from a 
! the counter of a previous run, if desired

      read(9,*) nSteps

! now read the number of atoms, replicas and degrees of freedom
      
      read(9,*) state_chain % nAtoms
      read(9,*) state_chain % nReplicas
      read(9,*) state_chain % nDegreesOfFreedom

! now read the thermostat information

      read(9,*) state_chain % nThermostats

      allocate( state_chain % thermostat( state_chain % nThermostats ) )
      
      do i=1, state_chain % nThermostats

         read(9,*) state_chain % thermostat(i) % mass,                     &
                   state_chain % thermostat(i) % position,                 & 
                   state_chain % thermostat(i) % momentum,                 & 
                   state_chain % thermostat(i) % a,                        & 
                   state_chain % thermostat(i) % C
         
      end do

! initialised energies

      state_chain % energy % thermostat = zero
      state_chain % H0 % set = .false.
      state_chain % H0 % E = zero
      
! read unit cell vectors
      read(9,*) LatticeVectors(1:3,1)
      read(9,*) LatticeVectors(1:3,2)
      read(9,*) LatticeVectors(1:3,3)

! read periodic logical flag
      read(9,*) periodic
      PBC = periodic

! read representaion flag
      read(9,*) representation % position
      read(9,*) representation % momentum
      read(9,*) representation % force      
      
! allocate structures

      allocate( mass( state_chain % nAtoms ) )
      allocate( position( 3, state_chain % nAtoms ) )
      allocate( momentum( 3, state_chain % nAtoms ) )
      
      allocate( state_chain % mass( state_chain % nAtoms  ) )
      allocate( state_chain % replica( state_chain % nReplicas ) )

! read masses of atoms in internal units

      do i=1, state_chain % nAtoms
         read(9,*) mass(i)
      end do

      state_chain % mass = mass

! read positions and momenta of every replica
      
      do irep = 1, state_chain % nReplicas

         do i=1, state_chain % nAtoms
            read(9,*) position(1:3,i)
         end do

         do i=1, state_chain % nAtoms
            read(9,*) momentum(1:3,i)
         end do

         ! create and fill every replica of the chain of states

         call createState( state_chain % replica(irep), state_chain % nAtoms )

         call initialiseCell( state_chain % replica(irep), LatticeVectors ) 

         state_chain % replica(irep) % position = position
         
         state_chain % replica(irep) % momentum = momentum

         state_chain % replica(irep) % representation = representation

         state_chain % replica(irep) % mass = mass

      end do

      close( unit = 9 )

      deallocate(mass)
      deallocate(position)
      deallocate(momentum)
      
      if ( debug ) write(*,*) 'Exiting readRestartChain()'

end subroutine readRestartChain
