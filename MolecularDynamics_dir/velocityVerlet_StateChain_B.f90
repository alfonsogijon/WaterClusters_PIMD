!*****************************************************************************
subroutine velocityVerlet_StateChain_B( state_chain )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: MolecularDynamics (PRIVATE routine)
!
!  Created on: March 2019 by Alfonso Gijón
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

      type ( StateChain_type ), intent (INOUT) :: state_chain

!*****************************************************************************
!  internal variables

      integer i

      real(dp) :: factor

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering VelocityVerlet_StateChain_B()'

! loop over replicas

      do i = 1, nReplicas

         call getStateChain( state_chain, iReplica = i, state = state_replica(i) )
         
         ! get the momenta and current force from state
         
         call get( state_replica(i), momentum = momentum_chain(:,:,i), &
              force = force_chain(:,:,i) )

         ! update momenta by the remaining half-step

         momentum_chain(:,:,i) = momentum_chain(:,:,i) + half_delta_t * force_chain(:,:,i)

         ! put the updated momentum in state

         call set( state_replica(i), momentum = momentum_chain(:,:,i) )

         call setStateChain( state_chain, iReplica = i, state = state_replica(i) )

      end do

      if ( debug ) write(*,*) 'Exiting VelocityVerlet_StateChain_B()'

end subroutine velocityVerlet_StateChain_B

