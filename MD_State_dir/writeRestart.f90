!*****************************************************************************
subroutine writeRestart( state, restartFileName, nSteps )
!*****************************************************************************
!
!  Project: MolecularDynamics   
!
!  Module:: State (PUBLIC routine)
!
!  Created on: Fri  2 Nov 11:05:47 2018  by Eduardo R. Hernandez
!
!  Last Modified: Thu 17 Jan 19:03:29 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
! 
! This subroutine writes to a file the current state of the system, with the
! purpose of serving as a continuation file for subsequent runs. 
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  arguments

      type ( state_type ), intent (IN) :: state  ! the state to be stored

      character ( len = 30 ), intent (IN) :: restartFileName

      integer, intent (IN), optional :: nSteps

!*****************************************************************************
!  local variables

      integer :: i
      integer :: n

!*****************************************************************************
!  begin subroutine 

      if ( debug ) write(*,*) 'Entering writeRestart()'

      open( unit = 9, file = restartFileName )

! if the number of steps has been passed, write it; else write a zero
! this is so that the continuation run can start counting steps from a 
! the counter of a previous run, if desired

      if ( present( nSteps ) ) then
         n = nSteps
      else
         n = 0
      end if

      write(9,*) n

! write the H0 constant, needed for calculating the constant of motion

      write(9,*) state % H0 % E

! write the representation flags

      write(9,*) state % representation % position,                    &
                 state % representation % momentum,                    &
                 state % representation % force

! now write the number of atoms, positions, forces, momenta and mass

      write(9,*) state % nAtoms
      write(9,*) state % nDegreesOfFreedom

      do i=1, state % nAtoms
          write(9,*) state % position(1:3,i)
      end do

      do i=1, state % nAtoms
          write(9,*) state % momentum(1:3,i)
      end do

      do i=1, state % nAtoms
          write(9,*) state % mass(i)
      end do

! write the current cell lattice vectors
! can be easily reconstructed

      write(9,*) state % Cell % a(1:3)
      write(9,*) state % Cell % b(1:3)
      write(9,*) state % Cell % c(1:3)

! write associated barostat variables

      write(9,*) state % barostat % mass

      write(9,*) state % barostat % momentum(1:3,1)
      write(9,*) state % barostat % momentum(1:3,2)
      write(9,*) state % barostat % momentum(1:3,3)

! now write the thermostat information

      write(9,*) state % nThermostats

      do i=1, state % nThermostats

          write(9,*) state % thermostat(i) % mass,                     &
                     state % thermostat(i) % position,                 & 
                     state % thermostat(i) % momentum,                 & 
                     state % thermostat(i) % a,                        & 
                     state % thermostat(i) % C

      end do

      close( unit = 9 )
      
      if ( debug ) write(*,*) 'Exiting writeRestart()'

end subroutine writeRestart
