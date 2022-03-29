!*****************************************************************************
subroutine setupStateChain( state_chain, nReplicas, temperature, nAtoms, &
     nThermostats, thermostatMass, a, b, c, mass, position0, lattice, periodic ) 
!*****************************************************************************
!
!  Project: MolecularDynamics 
!
!  Module: State (PUBLIC routine)
!
!  Created on: March 2019  by Alfonso Gijón
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This subroutine is called to create a new chain of states; each state instance
! is allocated, and filled
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  arguments

      type ( StateChain_type ), intent (INOUT) :: state_chain

      integer, intent (IN) :: nReplicas

      real (dp), intent(IN) :: temperature
      
      integer, intent (IN) :: nAtoms

      integer, intent (IN) :: nThermostats
 
      real (dp), dimension (:), intent (IN ):: thermostatMass

      real(dp) :: linear_momentum(3)
      
      real (dp), dimension(3), intent (IN) :: a
      real (dp), dimension(3), intent (IN) :: b
      real (dp), dimension(3), intent (IN) :: c

      real (dp), dimension (:), intent (IN) :: mass

      real (dp), dimension (:,:), intent (IN) :: position0

      logical, intent(IN), optional :: lattice

      logical, intent(IN), optional :: periodic

!*****************************************************************************
!  local variables

      integer :: i, j

      real (dp) :: aux, dseed

!*****************************************************************************
!  start subroutine

      if ( debug ) write(*,*) 'Entering setupStateChain()'
      
      ! assing nReplicas, nAtoms and nThermostats and allocate variables that depend on them

      state_chain % nReplicas = nReplicas

      state_chain % nDegreesOfFreedom = 3 * nAtoms * nReplicas - 1      
      
      state_chain % nAtoms = nAtoms 

      state_chain % nThermostats = nThermostats

      allocate( state_chain % mass( nAtoms ) )

      ! store the masses internally
      
      state_chain % mass = massConversionFactor * mass

      if (nThermostats .gt. 0) then
      
         allocate( state_chain % thermostat( nThermostats ) )

         do i = 1, state_chain % nThermostats
            state_chain % thermostat(i) % mass = massConversionFactor * thermostatMass(i)
            state_chain % thermostat(i) % position = one
            state_chain % thermostat(i) % momentum = zero
            state_chain % thermostat(i) % a = one
            state_chain % thermostat(i) % C = one
            write(*,*) 'thermostat mass ', i, state_chain % thermostat(i) % mass
         end do

      end if

      ! these are initialised to zero, they will be set properly later
      
      state_chain % energy % thermostat = zero

      state_chain % H0 % set = .false.
      state_chain % H0 % E = zero

      ! initialize random number generator

      dseed = 954131
      aux = randini(dseed)
      
      ! create and fill each state of the chain

      allocate( state_chain % replica( nReplicas ) )

      if ( present(lattice) ) then

         if ( present(periodic) ) then
         
            do i = 1, nReplicas

               call createState( state_chain % replica(i), nAtoms )

               call initialiseState( state_chain % replica(i), a, b, c, mass, position0, &
                    temperature = temperature, lattice = lattice, periodic = periodic )

            end do

         else

            do i = 1, nReplicas

               call createState( state_chain % replica(i), nAtoms )

               call initialiseState( state_chain % replica(i), a, b, c, mass, position0, &
                    temperature = temperature, lattice = lattice )

            end do
            
         end if

      else

         if ( present(periodic) ) then

            do i = 1, nReplicas

               call createState( state_chain % replica(i), nAtoms )

               call initialiseState( state_chain % replica(i), a, b, c, mass, position0, &
                    temperature = temperature, periodic = periodic )

            end do
            
         else
         
            do i = 1, nReplicas

               call createState( state_chain % replica(i), nAtoms )

               call initialiseState( state_chain % replica(i), a, b, c, mass, position0, &
                    temperature = temperature )

            end do

         end if

      end if

! Remove the total linear momentum of the ring polymer

      ! compute the total linear momentum
      
      linear_momentum = zero

      do i = 1, nReplicas
         
         linear_momentum(1) = linear_momentum(1) + &
              sum( state_chain % replica(i) % momentum(1,:) )
         linear_momentum(2) = linear_momentum(2) + &
              sum( state_chain % replica(i) % momentum(2,:) )
         linear_momentum(3) = linear_momentum(3) + &
              sum( state_chain % replica(i) % momentum(3,:) )         

      end do

      ! subtract part of the total momentum to each atom

      do i = 1, nReplicas
         do j = 1, nAtoms

            state_chain % replica(i) % momentum(:,j) = &
                 state_chain % replica(i) % momentum(:,j) - linear_momentum(:) &
                 / dfloat( nAtoms*nReplicas )
            
         end do
      end do      

end subroutine setupStateChain

