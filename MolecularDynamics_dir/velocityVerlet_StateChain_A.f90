!*****************************************************************************
subroutine velocityVerlet_StateChain_A( state_chain )
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

      integer i, j, k

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering velocityVerlet_StateChain_A()'

! loop over replicas and atoms, applying the velocity Verlet algorithm

      do i = 1, nReplicas

         call getStateChain( state_chain, iReplica = i, state = state_replica(i) )
         
         ! get from state the current positions, momenta and forces 
         ! we assume masses have not changed

         call get( state_replica(i), position = position_chain(:,:,i), &
              momentum = momentum_chain(:,:,i), force = force_chain(:,:,i) )

         ! now advance momenta by half a step

         momentum_chain(:,:,i) = momentum_chain(:,:,i) + &
              half_delta_t * force_chain(:,:,i)

         ! update the positions of all atoms with current momenta

         do j=1, nAtoms

            position_chain(1:3,j,i) = position_chain(1:3,j,i) +            &
                 delta_t * momentum_chain(1:3,j,i) / mass(j)
 
         end do
         
         ! momenta and positions have been updated, so put them in state

         call set( state_replica(i), position = position_chain(:,:,i), &
              momentum = momentum_chain(:,:,i) )

         call setStateChain( state_chain, iReplica = i, state = state_replica(i) )
         
      end do

      if ( debug ) write(*,*) 'Exiting velocityVerlet_StateChain_A()'

end subroutine velocityVerlet_StateChain_A
