!*****************************************************************************
subroutine BussiParrinello_B( state_chain )
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

      real(dp) :: determinant, determinant_1, total_mass, r2
      
      real(dp) :: linear_momentum(3), centre_of_mass(3), angular_momentum(3), angular_velocity(3)
      real(dp) :: inertia_tensor(3,3), matrix(3,3)

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering BussiParrinello_B()'

! loop over replicas and atoms, applying the Bussi-Parrinello algorithm 

      do i = 1, nReplicas

         call getStateChain( state_chain, iReplica = i, state = state_replica(i) )
         
         ! get from state the current positions, momenta and forces (at half timestep)

         call get( state_replica(i), momentum = momentum_chain(:,:,i), force = force_chain(:,:,i) )

         !do j=1, nAtoms

            !do k=1, 3

               ! advance momenta by the remaining half a step

               !momentum_chain(k,j,i) = c1 * ( momentum_chain(k,j,i) + &
               !     force_chain(k,j,i) * half_delta_t ) + c2(j) * rnorm()
               
            !end do
            
         !end do

         !!!! Alternative integration
         
         ! update momenta by the remaining half-step

         momentum_chain(:,:,i) = momentum_chain(:,:,i) + half_delta_t * force_chain(:,:,i)

         ! Action of the thermostat

         do j = 1, nAtoms

            do k = 1, 3

               momentum_chain(k,j,i) = c1 * c1 * momentum_chain(k,j,i) + &
                    c2(j) * dsqrt( c1*c1 + one ) * gasdev() !rnorm()

            end do

         end do
         
         ! momenta and positions have been updated, so put them in state

         !call set( state_replica(i), momentum = momentum_chain(:,:,i) )

         !call setStateChain( state_chain, iReplica = i, state = state_replica(i) )

      end do
      
! Remove the total linear momentum introduced by the action of the thermostat
      
      !go to 100
      
      ! compute the total linear momentum
      
      linear_momentum = zero
      linear_momentum(1) = linear_momentum(1) + sum( momentum_chain(1,:,:) )
      linear_momentum(2) = linear_momentum(2) + sum( momentum_chain(2,:,:) )
      linear_momentum(3) = linear_momentum(3) + sum( momentum_chain(3,:,:) )      

      ! subtract part of the total momentum to each atom

      do i = 1, nReplicas
         do j = 1, nAtoms
            momentum_chain(:,j,i) = momentum_chain(:,j,i) - &
                 linear_momentum(:) / dfloat( nAtoms*nReplicas )
         end do
      end do
!100   continue

      ! just checking
      
      ! linear_momentum = zero
      ! linear_momentum(1) = linear_momentum(1) + sum( momentum_chain(1,:,:) )
      ! linear_momentum(2) = linear_momentum(2) + sum( momentum_chain(2,:,:) )
      ! linear_momentum(3) = linear_momentum(3) + sum( momentum_chain(3,:,:) )
      ! print*, "Total Linear Momentum", linear_momentum

! Remove the total angular momentum introduced by the action of the thermostat
      
      ! we will refer momentarily the positions with respect to the centre of mass

      total_mass = dfloat(nReplicas) * sum( mass(:) )
      centre_of_mass = zero

      do i = 1, nReplicas
         do j = 1,  nAtoms

            centre_of_mass(1) = centre_of_mass(1) +                              &
                 mass(j) * position_chain(1,j,i)
            centre_of_mass(2) = centre_of_mass(2) +                              &
                 mass(j) * position_chain(2,j,i)
            centre_of_mass(3) = centre_of_mass(3) +                              &
                 mass(j) * position_chain(3,j,i)
            
         end do
      end do

      centre_of_mass = centre_of_mass / total_mass

      position_chain(1,:,:) = position_chain(1,:,:) - centre_of_mass(1)
      position_chain(2,:,:) = position_chain(2,:,:) - centre_of_mass(2)
      position_chain(3,:,:) = position_chain(3,:,:) - centre_of_mass(3)

!      go to 200      

      ! compute current angular momentum and intertia tensor

      angular_momentum = zero
      inertia_tensor = zero
      
      do i = 1, nReplicas
         do j = 1, nAtoms

            angular_momentum(1) = angular_momentum(1) +                  &
                 position_chain(2,j,i) * momentum_chain(3,j,i) -         &
                 position_chain(3,j,i) * momentum_chain(2,j,i)
            angular_momentum(2) = angular_momentum(2) +                  &
                 position_chain(3,j,i) * momentum_chain(1,j,i) -         &
                 position_chain(1,j,i) * momentum_chain(3,j,i)
            angular_momentum(3) = angular_momentum(3) +                  &
                 position_chain(1,j,i) * momentum_chain(2,j,i) -         &
                 position_chain(2,j,i) * momentum_chain(1,j,i)

            r2 = dot_product( position_chain(:,j,i), position_chain(:,j,i) )

            inertia_tensor(1,1) = inertia_tensor(1,1) + mass(j) *  &
                 ( r2 - position_chain(1,j,i) * position_chain(1,j,i) )
            inertia_tensor(1,2) = inertia_tensor(1,2) - mass(j) *  &
                 position_chain(1,j,i) * position_chain(2,j,i)  
            inertia_tensor(1,3) = inertia_tensor(1,3) - mass(j) *  &
                 position_chain(1,j,i) * position_chain(3,j,i)  

            inertia_tensor(2,2) = inertia_tensor(2,2) + mass(j) *  &
                 ( r2 - position_chain(2,j,i) * position_chain(2,j,i) )
            inertia_tensor(2,3) = inertia_tensor(2,3) - mass(j) *  &
                 position_chain(2,j,i) * position_chain(3,j,i)  
        
            inertia_tensor(3,3) = inertia_tensor(3,3) + mass(j) *  &
                 ( r2 - position_chain(3,j,i) * position_chain(3,j,i) )

            inertia_tensor(2,1) = inertia_tensor(1,2)
            inertia_tensor(3,1) = inertia_tensor(1,3)
            inertia_tensor(3,2) = inertia_tensor(2,3)
            
         end do
      end do

      determinant = inertia_tensor(1,1) * (                               &
           inertia_tensor(2,2) * inertia_tensor(3,3) -                    &
           inertia_tensor(2,3) * inertia_tensor(3,2) ) -                  &
           inertia_tensor(1,2) * (                                        & 
           inertia_tensor(2,1) * inertia_tensor(3,3) -                    &
           inertia_tensor(2,3) * inertia_tensor(3,1) ) +                  &
           inertia_tensor(1,3) * (                                        &
           inertia_tensor(2,1) * inertia_tensor(3,2) -                    &
           inertia_tensor(2,2) * inertia_tensor(3,1) )

      ! calculate the components of the angular velocity

      matrix = inertia_tensor
      matrix(1:3,1) = angular_momentum

      determinant_1 = matrix(1,1) * (                                          &
           matrix(2,2) * matrix(3,3) - matrix(2,3) * matrix(3,2) ) -      &
           matrix(1,2) * (                                            &
           matrix(2,1) * matrix(3,3) - matrix(2,3) * matrix(3,1) ) +      &
           matrix(1,3) * (                                            &
           matrix(2,1) * matrix(3,2) - matrix(2,2) * matrix(3,1) )

      angular_velocity(1) = determinant_1 / determinant
      
      matrix = inertia_tensor
      matrix(1:3,2) = angular_momentum

      determinant_1 = matrix(1,1) * (                                          &
           matrix(2,2) * matrix(3,3) - matrix(2,3) * matrix(3,2) ) -      &
           matrix(1,2) * (                                            &
           matrix(2,1) * matrix(3,3) - matrix(2,3) * matrix(3,1) ) +      &
           matrix(1,3) * (                                            &
           matrix(2,1) * matrix(3,2) - matrix(2,2) * matrix(3,1) )

      angular_velocity(2) = determinant_1 / determinant

      matrix = inertia_tensor
      matrix(1:3,3) = angular_momentum

      determinant_1 = matrix(1,1) * (                                          &
           matrix(2,2) * matrix(3,3) - matrix(2,3) * matrix(3,2) ) -      &
           matrix(1,2) * (                                            &
           matrix(2,1) * matrix(3,3) - matrix(2,3) * matrix(3,1) ) +      &
           matrix(1,3) * (                                            &
           matrix(2,1) * matrix(3,2) - matrix(2,2) * matrix(3,1) )

      angular_velocity(3) = determinant_1 / determinant

      ! add its opposite angular velocity to each atom
      
      angular_velocity = -angular_velocity

      do i = 1, nReplicas
         do j = 1, nAtoms

            momentum_chain(1,j,i) = momentum_chain(1,j,i) +                   &
                 mass(j) * (                                                  &
                 position_chain(3,j,i) * angular_velocity(2) -                &
                 position_chain(2,j,i) * angular_velocity(3) )

            momentum_chain(2,j,i) = momentum_chain(2,j,i) +                   &
                 mass(j) * (                                                  &
                 position_chain(1,j,i) * angular_velocity(3) -                &
                 position_chain(3,j,i) * angular_velocity(1) )
         
            momentum_chain(3,j,i) = momentum_chain(3,j,i) +                   &
                 mass(j) * (                                                  &
                 position_chain(2,j,i) * angular_velocity(1) -                &
                 position_chain(1,j,i) * angular_velocity(2) )
         
         end do
      end do
!200   continue

      ! just checking

      !angular_momentum = zero
      
      !do i = 1, nReplicas
      !   do j = 1, nAtoms

      !      angular_momentum(1) = angular_momentum(1) +                  &
      !           position_chain(2,j,i) * momentum_chain(3,j,i) -         &
      !           position_chain(3,j,i) * momentum_chain(2,j,i)
      !      angular_momentum(2) = angular_momentum(2) +                  &
      !           position_chain(3,j,i) * momentum_chain(1,j,i) -         &
      !           position_chain(1,j,i) * momentum_chain(3,j,i)
      !      angular_momentum(3) = angular_momentum(3) +                  &
      !           position_chain(1,j,i) * momentum_chain(2,j,i) -         &
      !           position_chain(2,j,i) * momentum_chain(1,j,i)

      !   end do
      !end do      

      !print*, "Total Angular Momentum", angular_momentum

      ! go back to absolute positions
      
      position_chain(1,:,:) = position_chain(1,:,:) + centre_of_mass(1)
      position_chain(2,:,:) = position_chain(2,:,:) + centre_of_mass(2)
      position_chain(3,:,:) = position_chain(3,:,:) + centre_of_mass(3)

! momenta and positions have been updated, so put them in state

      do i = 1, nReplicas
      
         call set( state_replica(i), momentum = momentum_chain(:,:,i) )

         call setStateChain( state_chain, iReplica = i, state = state_replica(i) )

      end do
      
      if ( debug ) write(*,*) 'Exiting BussiParrinello_B()'

end subroutine BussiParrinello_B
