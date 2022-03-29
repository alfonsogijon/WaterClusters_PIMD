!*****************************************************************************
subroutine set( state,                                                       &
                position, force, momentum, mass,                             &
                thermostat, barostat, energy,                                &
                Cell, representation )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: State (PUBLIC routine)
!
!  Created on: Fri  8 Jun 14:32:53 2018  by Eduardo R. Hernandez
!
!  Last Modified: Thu 17 Jan 19:01:27 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
! 
! This subroutine has the purpose of setting the values of variables contained
! in an instance of State to new numerical values. An object of type State
! must be passed as first argument; the remaining arguments are optional and
! identified by keyword, so that a single variable, a subset or the whole set
! can be specified in one go. The updated State object is then returned. 
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  arguments

      type ( state_type ), intent (INOUT) :: state   ! the state being modified

      real (dp), dimension ( 3, state % nAtoms ), intent (IN), optional :: position 
      real (dp), dimension ( 3, state % nAtoms ), intent (IN), optional :: force
      real (dp), dimension ( 3, state % nAtoms ), intent (IN), optional :: momentum 
      real (dp), dimension ( state % nAtoms ), intent (IN), optional :: mass

      type (thermostat_type), dimension ( state % nThermostats ), intent (IN), optional :: thermostat
      type (barostat_type), intent (IN), optional :: barostat

      type (energy_type), intent (IN), optional :: energy

      type (lattice_type), intent (IN), optional :: Cell

      type (representation_type), intent (IN), optional :: representation

!*****************************************************************************
!  local variables 

      integer :: n

!*****************************************************************************
!  start subroutine 

      if ( debug ) write(*,*) 'Entering set()'

      if ( present( position ) ) then

         if ( present( representation ) ) then

            if ( representation % position ) then
      
               state % representation % position = .true.  ! coordinates are input in cartesian rep.

            else

               state % representation % position = .false. ! coordinates are input in lattice rep. 

            end if

         else     ! default case; it is assumed coordinates provided in cartesian

            state % representation % position = .true.

         end if

         state % position = position

      end if

      if ( present( force ) ) then

         if ( present( representation ) ) then

            if ( representation % force ) then 
  
               state % representation % force = .true.  ! forces are input in cartesian rep.

            else

               state % representation % force = .false. ! they are in lattice rep.

            end if

         else 

            state % representation % force = .true.   ! default, they are in cartesian

         end if

         state % force = force

      end if

      if ( present( momentum ) ) then

         if ( present( representation ) ) then

            if ( representation % momentum ) then 

               state % representation % momentum = .true.  ! momenta are input in cartesian rep.

            else

               state % representation % momentum = .false.  ! nope, they are in lattice rep.

            end if

         else

            state % representation % momentum = .true.   ! default case, they are in cartesian

         end if

         state % momentum = momentum

      end if

      if ( present( mass ) ) state % mass = mass

      if ( present( thermostat ) ) then

         do n = 1, state % nThermostats
             state % thermostat(n) = thermostat(n)
         end do

      end if

      if ( present( barostat ) ) then
             
         state % barostat = barostat

      end if

      if ( present( energy ) ) then

         state % energy = energy

         if ( .not. state % H0 % set ) then

            state % H0 % E = state % energy % kinetic + state % energy % potential + &
                             state % energy % barostat + state % energy % thermostat

            state % H0 % set = .true.

         end if

      end if

      if ( present( Cell ) ) then

         state % Cell = Cell

         ! after the Cell has been changed/updated, we must re-calculate the reciprocal cell

         call setReciprocalCell( state )

      end if

      if ( debug ) write(*,*) 'Exiting set()'

end subroutine set

