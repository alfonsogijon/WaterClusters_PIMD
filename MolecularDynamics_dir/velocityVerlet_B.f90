!*****************************************************************************
subroutine velocityVerlet_B( state )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: MolecularDynamics (PRIVATE routine)
!
!  Created on: Wed 14 Mar 17:48:52 2018  by Eduardo R. Hernandez
!
!  Last Modified: Thu 17 Oct 13:38:15 2019
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

      integer :: i

      real (dp) :: kinetic

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering VelocityVerlet_B()'

! get the momenta and current force from state

      call get( state, energy = energy, force = force, momentum = momentum )

! update momenta by the remaining half-step

      momentum = momentum + half_delta_t * force

! calculate the new kinetic energy

      kinetic = zero

      do i=1, nAtoms

         kinetic = kinetic + dot_product( momentum(:,i), momentum(:,i) ) / mass(i)

      end do

      energy % kinetic = half * kinetic

! put the updated energy and momentum in state

      call set( state, energy = energy, momentum = momentum )

      if ( debug ) write(*,*) 'Exiting VelocityVerlet()'

end subroutine velocityVerlet_B

