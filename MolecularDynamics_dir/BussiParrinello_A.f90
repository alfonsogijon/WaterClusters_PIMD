!*****************************************************************************
subroutine BussiParrinello_A( state_chain )
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

      real(dp) :: p_tplus

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering BussiParrinello_A()'

! loop over replicas and atoms, applying the Bussi-Parrinello algorithm 

      do i = 1, nReplicas

         !call getStateChain( state_chain, iReplica = i, state = state_replica(i) )
         
         ! get from state the current positions, momenta and forces 
         ! we assume masses have not changed

         !call get( state_replica(i), position = position_chain(:,:,i), &
         !     momentum = momentum_chain(:,:,i), force = force_chain(:,:,i) )

         call getStateChain( state_chain, iReplica = i, position = position_chain(:,:,i), &
              momentum = momentum_chain(:,:,i), force = force_chain(:,:,i) )

         !do j=1, nAtoms

            !do k=1, 3

               !p_tplus = c1 * momentum_chain(k,j,i) + c2(j) * rnorm()

               ! advance position by a time step

               !position_chain(k,j,i) = position_chain(k,j,i) + p_tplus / mass(j) * delta_t + &
               !     force_chain(k,j,i) / mass(j) * delta_t2 / two
               
               ! advance momenta by half a step

               !momentum_chain(k,j,i) = momentum_chain(k,j,i) + p_tplus + &
               !     force_chain(k,j,i) * half_delta_t
               
            !end do
            
         !end do

         !!!!! Alternative integration

         ! advance momenta by half a step
               
         momentum_chain(:,:,i) = momentum_chain(:,:,i) + &
              half_delta_t * force_chain(:,:,i)

         ! update the positions of all atoms with current momenta

         do j=1, nAtoms

            position_chain(1:3,j,i) = position_chain(1:3,j,i) +            &
                 delta_t * momentum_chain(1:3,j,i) / mass(j)

         end do

         ! momenta and positions have been updated, so put them in state

         !call set( state_replica(i), position = position_chain(:,:,i), &
         !          momentum = momentum_chain(:,:,i) )

         !call setStateChain( state_chain, iReplica = i, state = state_replica(i) )

         call setStateChain( state_chain, iReplica = i, position = position_chain(:,:,i), &
              momentum = momentum_chain(:,:,i) )

      end do

      if ( debug ) write(*,*) 'Exiting BussiParrinello_A()'

end subroutine BussiParrinello_A
