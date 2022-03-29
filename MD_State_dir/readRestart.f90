!*****************************************************************************
subroutine readRestart( state, restartFile, nSteps )
!*****************************************************************************
!
!  Project: MolecularDynamics
! 
!  Module: State (PUBLIC routine)
!
!  Created on: Fri  2 Nov 12:33:45 2018  by Eduardo R. Hernandez
!
!  Last Modified: Wed 10 Apr 16:09:57 2019
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

      type ( state_type ), intent (INOUT) :: state  ! the state to be read

      character ( len = 30 ), intent (IN) :: restartFile

      integer, intent (INOUT), optional :: nSteps

!*****************************************************************************
!  local variables

      integer :: i

      real (dp), dimension (3,3) :: LatticeVec

!*****************************************************************************
!  begin subroutine 

      if ( debug ) write(*,*) 'Entering readRestart()'

      open( unit = 9, file = restartFile )

! read the number of steps integrated in previous runs
! this is so that the continuation run can start counting steps from a 
! the counter of a previous run, if desired

      read(9,*) nSteps

! read the H0 constant, needed for calculating the constant of motion

      read(9,*) state % H0 % E
      state % H0 % set = .true.

! read the representation flags

      read(9,*) state % representation % position,                    &
                state % representation % momentum,                    &
                state % representation % force

! now read the number of atoms, positions, forces, momenta and mass

      read(9,*) state % nAtoms
      read(9,*) state % nDegreesOfFreedom

      do i=1, state % nAtoms
          read(9,*) state % position(1:3,i)
      end do

      do i=1, state % nAtoms
          read(9,*) state % momentum(1:3,i)
      end do

      do i=1, state % nAtoms
          read(9,*) state % mass(i)
      end do

! read the cell vectors; with this the cell
! can be easily reconstructed

      read(9,*) LatticeVec(1:3,1)
      read(9,*) LatticeVec(1:3,2)
      read(9,*) LatticeVec(1:3,3)

      call initialiseCell( state, LatticeVec )

! read associated barostat variables

      read(9,*) state % barostat % mass

      read(9,*) state % barostat % momentum(1:3,1)
      read(9,*) state % barostat % momentum(1:3,2)
      read(9,*) state % barostat % momentum(1:3,3)

! now read the thermostat information

      read(9,*) state % nThermostats

      do i=1, state % nThermostats

         read(9,*) state % thermostat(i) % mass,                     &
                   state % thermostat(i) % position,                 & 
                   state % thermostat(i) % momentum,                 & 
                   state % thermostat(i) % a,                        & 
                   state % thermostat(i) % C

      end do

      close( unit = 9 )
      
      if ( debug ) write(*,*) 'Exiting readRestart()'

end subroutine readRestart
