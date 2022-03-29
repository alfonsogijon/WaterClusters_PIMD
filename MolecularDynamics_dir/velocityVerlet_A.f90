!*****************************************************************************
subroutine velocityVerlet_A( state )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: MolecularDynamics (PRIVATE routine)
!
!  Created on: Wed 14 Mar 17:48:52 2018  by Eduardo R. Hernandez
!
!  Last Modified: Tue 26 Mar 14:09:21 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! Implements the update of momenta to half-step and positions to full step
! according to the classical Velocity-Verlet method for NVE (microcanonical)
! simulations. The rest of the algorithm, for the update to full-step of the
! momenta, after force evaluation at the new positions, is implemented in the
! sequel subroutine VelocityVerlet_B.
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  variables

      type ( state_type ), intent (INOUT) :: state

!*****************************************************************************
!  internal variables

      integer i

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering VelocityVerlet_A()'

! get from state the current positions, momenta and forces 
! we assume masses have not changed

      call get( state, position = position, momentum = momentum, force = force )

! now advance momenta by half a step

      momentum = momentum + half_delta_t * force

! update the positions of all atoms with current momenta

      do i=1, nAtoms

         position(1:3,i) = position(1:3,i) +            &
                 delta_t * momentum(1:3,i) / mass(i)
 
      end do

! momenta and positions have been updated, so put them in state

      call set( state, position = position, momentum = momentum )

      if ( debug ) write(*,*) 'Exiting VelocityVerlet_A()'

end subroutine velocityVerlet_A

