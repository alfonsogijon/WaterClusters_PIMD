!*****************************************************************************
subroutine get( state, nAtoms, nDegreesOfFreedom, nThermostats,              &
                position, force, momentum, mass,                             &
                thermostat, barostat, energy, H0,                            &
                Cell, reciprocalCell, stress, representation  )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: State (PUBLIC routine)
!
!  Created on: Fri  8 Jun 14:53:35 2018  by Eduardo R. Hernandez
!
!  Last Modified: Wed 22 May 18:29:57 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This subroutine extracts from an object of type State_type specific values
! as requested by the caller. The first variable, which is mandatory, is the
! state object being interrogated. Then follows a keyword-list of desired 
! values (temperature, pressure, etc), all of which are optional, that are
! returned upon request.
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  arguments

      type ( state_type ), intent (INOUT) :: state   ! the state being interrogated

      integer, intent (OUT), optional :: nAtoms
      integer, intent (OUT), optional :: nDegreesOfFreedom
      integer, intent (OUT), optional :: nThermostats

      real (dp), intent (OUT), optional :: H0

      real (dp), dimension (:,:), intent (OUT), optional :: position
      real (dp), dimension (:,:), intent (OUT), optional :: force
      real (dp), dimension (:,:), intent (OUT), optional :: momentum
      real (dp), dimension (:), intent (OUT), optional :: mass

      type ( thermostat_type ), dimension (:), intent (OUT), optional :: thermostat
      type ( barostat_type ), intent (OUT), optional :: barostat

      type ( energy_type ), intent (OUT), optional :: energy

      type ( lattice_type ), intent (OUT), optional :: Cell
      type ( lattice_type ), intent (OUT), optional :: reciprocalCell

      type ( stress_type ), intent (OUT), optional :: stress

      type ( representation_type ), intent (IN), optional :: representation  

!*****************************************************************************
!  variables

      integer :: i

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering get()'

      if ( present( nAtoms ) ) then

         nAtoms = state % nAtoms

      end if

      if ( present( nDegreesOfFreedom ) ) then

         nDegreesOfFreedom = state % nDegreesOfFreedom

      end if

      if ( present( nThermostats ) ) then

         nThermostats = state % nThermostats

      end if

      if ( present( position ) ) then

         !if ( present( representation ) ) then

         !   if ( representation % position ) then    ! coordinates are needed in cartesian rep.

         !      call transformCoordinates2Cartesian( state )

         !   else                                     ! coordinates are needed in lattice rep.

         !      call transformCoordinates2Lattice( state )

         !   end if

         !else                                        ! default; make sure they are in cartesian

         !   call transformCoordinates2Cartesian( state )

         !end if

         position = state % position    

      end if

      if ( present( force ) ) then

         !if ( present( representation ) ) then

         !   if ( representation % force ) then    ! forces are needed in cartesian rep.
     
         !      call transformForces2Cartesian( state )

         !   else                                  ! forces are needed in lattice rep.

         !      call transformForces2Lattice( state )

         !   end if

         !else     ! default, we return them in cartesian, so make sure they are in that rep.

         !   call transformForces2Cartesian( state )

         !end if

         force = state % force    ! copy forces

      end if

      if ( present( momentum ) ) then

         !if ( present( representation ) ) then

         !   if ( representation % momentum ) then    ! momenta are needed in cartesian rep.

         !      call transformMomenta2Cartesian( state )

         !   else                                     ! momenta are needed in lattice rep.

         !      call transformMomenta2Lattice( state )

         !   end if

         !else                                        ! default case, we return them in cartesian

         !   call transformMomenta2Cartesian( state )

         !end if

         momentum = state % momentum    ! copy momenta

      end if

      if ( present( mass ) ) then

         mass = state % mass

      end if

      if ( present( thermostat ) ) then

         do i = 1, state % nThermostats
 
            thermostat(i) = state % thermostat(i) 

         end do

      end if

      if ( present( barostat ) ) then

         barostat = state % barostat

      end if

      if ( present( energy ) ) then

         energy = state % energy

      end if

      if ( present( H0 ) ) then

         H0 = state % H0 % E

      end if

      if ( present( Cell ) ) then

         Cell = state % Cell

      end if

      if ( present( reciprocalCell ) ) then

         reciprocalCell = state % reciprocalCell

      end if

      if ( present( stress ) ) then

         stress = state % stress

      end if

      if ( debug ) write(*,*) 'Exiting get()'

end subroutine get

