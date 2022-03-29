!*****************************************************************************
subroutine generateMomenta( state, target_temperature )
!*****************************************************************************
!
!  Project: MolecularDynamics 
!
!  Module: State (PRIVATE routine)
!
!  Created on: Wed 11 Apr 12:27:49 2018  by Eduardo R. Hernandez
!
!  Last Modified: Thu 17 Jan 19:04:17 2019
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

      type ( state_type ), intent (INOUT) :: state

      real (dp), intent (IN) :: target_temperature

!**************************************************************************
!  Local variables

      integer :: i, n

      real (dp) :: r, r2, rnum, total_mass, variance, v1, v2
      real (dp) :: determinant, determinant_1, E_kin, factor, temperature

      real (dp) :: angular_momentum( 3 ), angular_velocity( 3 )
      real (dp) :: centre_of_mass( 3 )
      real (dp) :: cm_velocity( 3 )
      real (dp) :: inertia_tensor(3,3)
      real (dp) :: matrix(3,3)

!**************************************************************************
!  Start of subroutine

      if ( debug ) write(*,*) 'Entering generateMomenta()'

      total_mass = zero

! loop over atoms and assign each atom a random velocity
! according to its mass and the target temperature

! let's try our luck

      do i=1, state % nAtoms

         total_mass = total_mass + state % mass(i)
         variance = dsqrt( boltzmann_k * target_temperature / state % mass(i) )

         do n= 1,3

            r = two
            do while ( r .gt. one )
               !call random_number( rnum )
               rnum = randu()
               v1 = two * rnum - one
               !call random_number( rnum )
               rnum = randu()
               v2 = two * rnum - one
               r = v1 * v1 + v2 * v2
            end do
            r = variance * v1 * sqrt( -two * log( r ) / r )
            state % momentum(n,i) = state % mass(i) * r
         end do
      end do

! now calculate the amount of velocity to take out for each atom

      do n=1, 3
         cm_velocity( n ) = zero
         do i=1, state % nAtoms
            cm_velocity(n) = cm_velocity(n) + state % momentum(n,i)
         end do
         cm_velocity(n) = cm_velocity(n) / total_mass
      end do

! now substract from each atom the velocity necessary so as
! to ensure that the total momentum is zero

      if ( state % nAtoms .gt. 1) then

         do n=1, 3
            do i=1, state % nAtoms
               state % momentum(n,i) = state % momentum(n,i) -                  &
                    state % mass(i) * cm_velocity(n)
            end do
         end do

      end if
         
! calculate the instantaneous temperature

      E_kin = zero

      do i=1, state % nAtoms

         factor = half / state % mass(i)

         E_kin = E_kin + factor * (                                          &
           dot_product( state % momentum(1:3,i), state % momentum(1:3,i) ) )

      end do

      if ( state % natoms .gt. 1) then
         temperature = two * E_kin / ( ( three * dfloat(state % nAtoms) - three ) * boltzmann_k )
      else
         temperature = two * E_kin / ( ( three * dfloat(state % nAtoms) ) * boltzmann_k )
      end if

! calculate the scaling factor and scale the velocities accordingly
     
      !write(*,*) 'temperature = ', temperature, target_temperature

      factor = dsqrt( target_temperature / temperature )

      state % momentum = factor * state % momentum

      E_kin = zero

      do i=1, state % nAtoms

         factor = half / state % mass(i)

         E_kin = E_kin + factor * (                                          &
           dot_product( state % momentum(1:3,i), state % momentum(1:3,i) ) )

      end do

      if ( state % natoms .gt. 1) then
         temperature = two * E_kin / ( ( three * dfloat(state % nAtoms) - three ) * boltzmann_k )
      else
         temperature = two * E_kin / ( ( three * dfloat(state % nAtoms) ) * boltzmann_k )
      end if

! just checking

      do n=1, 3
         do i=1, state % nAtoms
            cm_velocity(n) = cm_velocity(n) + state % momentum(n,i) 
         end do
         cm_velocity(n) = cm_velocity(n) / total_mass
      end do

      !write(*,*) 'centre of mass velocity = ', cm_velocity(1:3), temperature
      !write(*,*) 'kinetic energy = ', E_kin
      !write(*,*) 'momentum 1 ', state % momentum(1:3,1)

! if the system is not periodic in any direction (i.e. cluster or molecule)
! then we must make sure that it does have zero angular momentum (no rotation)
 
      if ( .not. PBC ) then

! we will refer momentarily the positions with respect to the centre of mass

        centre_of_mass = zero

        do i=1, state % nAtoms 

           centre_of_mass(1) = centre_of_mass(1) +                              &
               state % mass(i) * state % position(1,i)
           centre_of_mass(2) = centre_of_mass(2) +                              &
               state % mass(i) * state % position(2,i)
           centre_of_mass(3) = centre_of_mass(3) +                              &
               state % mass(i) * state % position(3,i)

        end do

        centre_of_mass = centre_of_mass / total_mass

        state % position(1,1:state % nAtoms) =                                      &
             state % position(1,1:state % nAtoms) - centre_of_mass(1)
        state % position(2,1:state % nAtoms) =                                      &
             state % position(2,1:state % nAtoms) - centre_of_mass(2)
        state % position(3,1:state % nAtoms) =                                      &
             state % position(3,1:state % nAtoms) - centre_of_mass(3)

        angular_momentum = zero
        inertia_tensor = zero

        do i=1, state % nAtoms

           angular_momentum(1) = angular_momentum(1) +                       &
               ( state % position(2,i) * state % momentum(3,i) -             &
                 state % position(3,i) * state % momentum(2,i) )
           angular_momentum(2) = angular_momentum(2) +                       &
               ( state % position(3,i) * state % momentum(1,i) -             &
                 state % position(1,i) * state % momentum(3,i) )
           angular_momentum(3) = angular_momentum(3) +                       &
               ( state % position(1,i) * state % momentum(2,i) -             &
                 state % position(2,i) * state % momentum(1,i) )

           r2 = dot_product( state % position(1:3,i), state % position(1:3,i))

           inertia_tensor(1,1) = inertia_tensor(1,1) + state % mass(i) *  &
                         ( r2 - state % position(1,i) * state % position(1,i) )
           inertia_tensor(1,2) = inertia_tensor(1,2) - state % mass(i) *  &
                           state % position(1,i) * state % position(2,i)  
           inertia_tensor(1,3) = inertia_tensor(1,3) - state % mass(i) *  &
                           state % position(1,i) * state % position(3,i)  

           inertia_tensor(2,1) = inertia_tensor(2,1) - state % mass(i) *  &
                           state % position(2,i) * state % position(1,i)  
           inertia_tensor(2,2) = inertia_tensor(2,2) + state % mass(i) *  &
                         ( r2 - state % position(2,i) * state % position(2,i) )
           inertia_tensor(2,3) = inertia_tensor(2,3) - state % mass(i) *  &
                           state % position(2,i) * state % position(3,i)  
        
           inertia_tensor(3,1) = inertia_tensor(3,1) - state % mass(i) *  &
                           state % position(3,i) * state % position(1,i)  
           inertia_tensor(3,2) = inertia_tensor(3,2) - state % mass(i) *  &
                           state % position(3,i) * state % position(2,i)  
           inertia_tensor(3,3) = inertia_tensor(3,3) + state % mass(i) *  &
                         ( r2 - state % position(3,i) * state % position(3,i) )

        end do

        determinant = inertia_tensor(1,1) * (                                    &
                  inertia_tensor(2,2) * inertia_tensor(3,3) -                    &
                  inertia_tensor(2,3) * inertia_tensor(3,2) ) -                  &
                      inertia_tensor(1,2) * (                                    &
                  inertia_tensor(2,1) * inertia_tensor(3,3) -                    &
                  inertia_tensor(2,3) * inertia_tensor(3,1) ) +                  &
                      inertia_tensor(1,3) * (                                    &
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
   
        angular_velocity = -angular_velocity

        do i=1, state % nAtoms

           state % momentum(1,i) = state % momentum(1,i) +                   &
                        state % mass(i) * (                           &
                state % position(3,i) * angular_velocity(2) -                &
                state % position(2,i) * angular_velocity(3) )

           state % momentum(2,i) = state % momentum(2,i) +                   &
                        state % mass(i) * (                           &
                state % position(1,i) * angular_velocity(3) -                &
                state % position(3,i) * angular_velocity(1) )

           state % momentum(3,i) = state % momentum(3,i) +                   &
                        state % mass(i) * (                           &
                state % position(2,i) * angular_velocity(1) -                &
                state % position(1,i) * angular_velocity(2) )

        end do

! just checking 

        angular_momentum = zero

        do i=1, state % nAtoms

           angular_momentum(1) = angular_momentum(1) +                       &
               ( state % position(2,i) * state % momentum(3,i) -             &
                 state % position(3,i) * state % momentum(2,i) )
           angular_momentum(2) = angular_momentum(2) +                       &
               ( state % position(3,i) * state % momentum(1,i) -             &
                 state % position(1,i) * state % momentum(3,i) )
           angular_momentum(3) = angular_momentum(3) +                       &
               ( state % position(1,i) * state % momentum(2,i) -             &
                 state % position(2,i) * state % momentum(1,i) )

        end do

        !write(*,*) 'angular momentum components ', angular_momentum(1:3)

        state % position(1,1:state % nAtoms) =                                      &
             state % position(1,1:state % nAtoms) + centre_of_mass(1)
        state % position(2,1:state % nAtoms) =                                      &
             state % position(2,1:state % nAtoms) + centre_of_mass(2)
        state % position(3,1:state % nAtoms) =                                      &
             state % position(3,1:state % nAtoms) + centre_of_mass(3)

      end if

! finally, calculate the kinetic energy and store it in type energy

      E_kin = zero

      do i=1, state % nAtoms

         factor = half / state % mass(i)

         E_kin = E_kin + factor * (                                          &
           dot_product( state % momentum(1:3,i), state % momentum(1:3,i) ) )

      end do

      state % energy % kinetic = E_kin

      state % representation % momentum = .true. ! momenta are in cartesian rep.

      if ( debug ) write(*,*) 'Exiting generateMomenta()'

end subroutine generateMomenta

