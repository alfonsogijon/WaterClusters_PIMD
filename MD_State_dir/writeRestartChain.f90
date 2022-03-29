!*****************************************************************************
subroutine writeRestartChain( state_chain, restartFileName, nSteps )
!*****************************************************************************
!
!  Project: MolecularDynamics   
!
!  Module:: State (PUBLIC routine)
!
!  Created on: May 201 by Alfonso Gijón
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

      type ( StateChain_type ), intent (IN) :: state_chain  ! the state to be stored

      character ( len = 30 ), intent (IN) :: restartFileName

      integer, intent (IN), optional :: nSteps

!*****************************************************************************
!  local variables

      integer :: i, irep
      integer :: n

!*****************************************************************************
!  begin subroutine 

      if ( debug ) write(*,*) 'Entering writeRestartChain()'

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

! write the number of atoms, replicas and thermostats

      write(9,*) state_chain % nAtoms
      write(9,*) state_chain % nReplicas
      write(9,*) state_chain % nDegreesOfFreedom

! now write the thermostat information

      write(9,*) state_chain % nThermostats

      do i=1, state_chain % nThermostats

          write(9,*) state_chain % thermostat(i) % mass,                     &
                     state_chain % thermostat(i) % position,                 & 
                     state_chain % thermostat(i) % momentum,                 & 
                     state_chain % thermostat(i) % a,                        & 
                     state_chain % thermostat(i) % C

      end do

! write the current cell lattice vectors

      write(9,*) state_chain % replica(1) % Cell % a(1:3)
      write(9,*) state_chain % replica(1) % Cell % b(1:3)
      write(9,*) state_chain % replica(1) % Cell % c(1:3)

! write periodic flag

      write(9,*) pbc

! write the representation flags
      
      write(9,*) state_chain % replica(1) % representation % position
      write(9,*) state_chain % replica(1) % representation % momentum
      write(9,*) state_chain % replica(1) % representation % force
      
! write masses of atoms in internal units
      do i=1, state_chain % nAtoms
         write(9,*) state_chain % mass(i)
      end do

! now write the positions and momenta in internal units

      do irep = 1, state_chain % nReplicas
      
         do i=1, state_chain % nAtoms
            write(9,*) state_chain % replica(irep) % position(1:3,i)
         end do

         do i=1, state_chain % nAtoms
            write(9,*) state_chain % replica(irep) % momentum(1:3,i)
         end do

      end do

! write associated barostat variables

      !write(9,*) state % barostat % mass

      !write(9,*) state % barostat % momentum(1:3,1)
      !write(9,*) state % barostat % momentum(1:3,2)
      !write(9,*) state % barostat % momentum(1:3,3)

      close( unit = 9 )
      
      if ( debug ) write(*,*) 'Exiting writeRestartChain()'

end subroutine writeRestartChain
