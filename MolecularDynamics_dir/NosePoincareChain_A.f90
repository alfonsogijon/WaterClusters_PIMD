!*****************************************************************************
subroutine  NosePoincareChain_A( state )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: MolecularDynamics (PRIVATE routine)
!
!  Created on: Mon 23 Jul 11:32:54 2018  by Eduardo R. Hernandez
!
!  Last Modified: Wed 10 Apr 16:28:39 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This routine implements the generalisation of the Nosé-Poincaré thermostat
! to Recursive-Chain thermostats described by Leimkhuler and Sweet, SIAM J. 
! Applied Dynamical Systems Vol. 4, 187-216 (2005). When the chain is just
! one thermostat long (the physical system is coupled to a single thermostat)
! the use of routines NosePoincare_A,B is preferred; this one is used when
! the number of thermostats in the chain is larger than 1, but is otherwise
! similar in structure to NosePoincare_A,B.
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
      integer :: nThermostats

      integer, parameter :: nIterationsMax = 20

      real (dp), parameter :: tolerance = 1.0e-10_dp

      real (dp) :: kinetic
      real (dp) :: factor
      real (dp) :: thermostatProductOld, thermostatProduct

      type ( representation_type ) :: representation 

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering NosePoincareChain_A()'

! make sure that coordinates, momenta and forces are obtained in cartesian representation

      representation % position = .true.
      representation % momentum = .true.
      representation % force = .true.

! get from state all the required variables; of course we assume masses are unchanged

      call get( state, nThermostats = nThermostats,                      &
                position = position, momentum = momentum, force = force, &
                thermostat = thermostat, energy = energy, H0 = H0,       &
                representation = representation )

! get the product of thermostat positions

      thermostatProduct = product( thermostat % position )

! first advance momenta by half a step

      momentum = momentum +                            &
            half_delta_t * thermostatProduct * force

! now calculate kinetic energy with new momenta at half-step

      kinetic = 0.0_dp

      do i = 1, nAtoms

         factor = half / mass(i) 

         kinetic = kinetic + factor *   &
          dot_product( momentum(1:3,i), momentum(1:3,i) )

      end do

      energy % kinetic = kinetic / ( thermostatProduct * thermostatProduct )

! now advance the momenta of thermostats in the chain y half a step.
! due to the coupled-nature of the equations, this requires an iterative
! procedure, so we are going to do this in a subroutine

      call UpdateThermostatChainMomenta

! now advance the thermostat position to full-step; this is also going to require
! an iterative procedure, so we will pack it into the following routine

      call UpdateThermostatChainPositions

! now we can advance atomic positions

      thermostatProductOld = thermostatProduct
      thermostatProduct = product( thermostat % position )

      factor = half_delta_t * ( one / thermostatProductOld + one / thermostatProduct )

! update the positions of all atoms with current momenta

      do i=1, nAtoms

         position(1:3,i) = position(1:3,i) + factor * momentum(1:3,i) / mass(i)
 
      end do

! update the value of the kinetic energy (ionic momenta not changed, but thermostat position has)

      energy % kinetic = kinetic / ( thermostatProduct * thermostatProduct )

! finally, store the new coordinates and half-step momenta in state

      ! call set( state, position = position, momentum = momentum, thermostat = thermostat )
      call set( state, position = position, momentum = momentum, thermostat = thermostat, energy = energy )

      if ( debug ) write(*,*) 'Exiting NosePoincareChain_A()'

contains

!*****************************************************************************
subroutine UpdateThermostatChainPositions
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Created on: Mon  23 August 16:35:05 2018  by Eduardo R. Hernandez
!
!  Last Modified: Mon 23 August 16:44:51 2018
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! Because of the coupled dependence of thermostat positions in the Nosé-Poincaré
! recursive thermostat chain, thermostat positions need to be updated 
! in an iterative fashion. This subroutine implements the required procedure.
!
!*****************************************************************************
!  No implicits please!

      implicit none

!*****************************************************************************
!  shared variables

!*****************************************************************************
!  local variables

      logical :: converged

      integer :: n
      integer :: nIterations

      real (dp) :: prod
      real (dp) :: thermostatProduct, thermostatProductNew

      real (dp), allocatable, dimension(:) :: denominator
      real (dp), allocatable, dimension(:) :: denominatorNew
      real (dp), allocatable, dimension(:) :: difference
      real (dp), allocatable, dimension(:) :: thermostatPositionNew
      real (dp), allocatable, dimension(:) :: thermostatPositionTemp

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering UpdateThermostatChainPositions()'

      allocate( thermostatPositionNew( nThermostats ) )
      allocate( thermostatPositionTemp( nThermostats ) )
      allocate( denominator( nThermostats ) )
      allocate( denominatorNew( nThermostats ) )
      allocate( difference( nThermostats ) )

      thermostatPositionNew = thermostat % position  ! initialise
      thermostatPositionTemp = thermostat % position

      if ( nThermostats .gt. 1 ) then

         do n=1, nThermostats-1
         
            prod = product( thermostat(n+1:nThermostats) % position )

            denominator(n) = one / ( thermostat(n) % mass * prod * prod )

         end do

         denominator(nThermostats) = one / thermostat(nThermostats) % mass

      else 

         denominator(1) = one / thermostat(1) % mass

      end if

      denominatorNew = denominator

      thermostatProduct = product( thermostat % position )
      thermostatProductNew = thermostatProduct

      nIterations = 0

      do        ! we loop until convergence has been attained

         thermostatPositionNew = thermostat % position + half_delta_t * (     &
              thermostatProduct * denominator +                               &
              thermostatProductNew * denominatorNew ) * thermostat % momentum

         difference = thermostatPositionNew - thermostatPositionTemp

         ! check for convergence

         converged = .true.

         do n=1, nThermostats

            if ( dabs( difference(n) ) .gt. tolerance ) converged = .false.

         end do

         if ( .not. converged ) then

            thermostatPositionTemp = thermostatPositionNew

            thermostatProductNew = product( thermostatPositionTemp )

            if ( nThermostats .gt. 1 ) then

               do n=1, nThermostats-1
         
                  prod = product( thermostatPositionTemp(n+1:nThermostats) )

                  denominatorNew(n) = one / ( thermostat(n) % mass * prod * prod )

               end do 

               denominatorNew(nThermostats) = one / thermostat(nThermostats) % mass

            else

               denominatorNew(1) = one / thermostat(1) % mass

            end if

         else

            exit               ! yeah, we have converged!

         end if

         nIterations = nIterations + 1

         if ( nIterations > nIterationsMax ) then

            write(*,*) 'Maximum number of interations exceeded in '
            write(*,*) 'UpdateThermostatChainPositions: stopping '

            stop

         end if

      end do

      thermostat % position = thermostatPositionTemp

      deallocate( thermostatPositionNew )
      deallocate( thermostatPositionTemp )
      deallocate( denominator )
      deallocate( denominatorNew )
      deallocate( difference )

      if ( debug ) write(*,*) 'Exiting UpdateThermostatChainPositions()'

end subroutine UpdateThermostatChainPositions

!*****************************************************************************
subroutine UpdateThermostatChainMomenta
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Created on: Mon 23 August 18:05:05 2018  by Eduardo R. Hernandez
!
!  Last Modified: Mon 23 August 18:44:51 2018
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! Because of the coupled dependence of thermostat momenta in the Nosé-Poincaré
! recursive thermostat chain, thermostat momenta need to be updated 
! in an iterative fashion. This subroutine implements the required procedure.
!
!*****************************************************************************
!  No implicits please!

      implicit none

!*****************************************************************************
!  shared variables

!*****************************************************************************
!  local variables

      logical :: allConverged
      logical, dimension (:), allocatable :: converged

      integer :: n
      integer :: nIterations

      real (dp) :: deltaH
      real (dp) :: kbT
      real (dp) :: prod
      real (dp) :: thermostatEnergyTemp
      real (dp) :: thermostatProduct

      real (dp), dimension (:), allocatable :: denominator
      real (dp), dimension (:), allocatable :: potential
      real (dp), dimension (:), allocatable :: dpotential
      real (dp), dimension (:), allocatable :: thermostatMomentumNew
      real (dp), dimension (:), allocatable :: thermostatMomentumTemp
      real (dp), dimension (:), allocatable :: thermostatE
      real (dp), dimension (:), allocatable :: workvec1
      real (dp), dimension (:), allocatable :: workvec2

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering UpdateThermostatChainMomenta()'

      allocate( denominator( nThermostats ) )
      allocate( converged( nThermostats ) )
      allocate( potential( nThermostats ) )
      allocate( dpotential( nThermostats ) )
      allocate( thermostatMomentumNew( nThermostats ) )
      allocate( thermostatMomentumTemp( nThermostats ) )
      allocate( thermostatE( nThermostats ) )
      allocate( workvec1( nThermostats ) )
      allocate( workvec2( nThermostats ) )

      thermostatMomentumNew = thermostat % momentum
      thermostatMomentumTemp = thermostat % momentum

      kbT = boltzmann_k * simulation % temperature

      thermostatProduct = product( thermostat % position )

      if ( nThermostats .gt. 1 ) then

         do n=1, nThermostats-1

            prod = product( thermostat(n+1:nThermostats) % position )

            denominator(n) = half / ( thermostat(n) % mass * prod * prod )

            if ( n .eq. 1 ) then
               potential(n) = float( nDegreesOfFreedom ) * kbT * log( thermostat(n) % position )
               dpotential(n) =  -float( nDegreesOfFreedom ) * kbT 
            else
               potential(n) = float( nDegreesOfFreedom + n - 1 ) * kbT * log( thermostat(n) % position ) + &
                   ( thermostat(n) % a - thermostat(n) % position ) *                    &
                   ( thermostat(n) % a - thermostat(n) % position ) /                    &
                            ( two * thermostat(n) % C )
               dpotential(n) = -float( nDegreesOfFreedom + n - 1 ) * kbT +               &
                   ( thermostat(n) % a - thermostat(n) % position ) *                    &
                     thermostat(n) % position / thermostat(n) % C
            end if

         end do

         denominator(nThermostats) = half / thermostat(nThermostats) % mass
         potential(nThermostats) =                                                       &
          float( nDegreesOfFreedom + nThermostats - 1 ) * kbT * log( thermostat(nThermostats) % position ) + &
                   ( thermostat(nThermostats) % a - thermostat(nThermostats) % position ) *                  &
                   ( thermostat(nThermostats) % a - thermostat(nThermostats) % position ) /                  &
                            ( two * thermostat(nThermostats) % C )
         dpotential(nThermostats) =                                                      &
          -float( nDegreesOfFreedom + nThermostats - 1 ) * kbT +                         &
                   ( thermostat(nThermostats) % a - thermostat(nThermostats) % position ) *                  &
                     thermostat(nThermostats) % position / thermostat(nThermostats) % C

      else

         denominator(1) = half / thermostat(1) % mass
         potential(1) = float( nDegreesOfFreedom ) * kbT * log( thermostat(1) % position )
         dpotential(1) = -float( nDegreesOfFreedom ) * kbT 

      end if

      ! before updating the momenta, evaluate the thermostat energy with the current thermostat momenta

      call calculateThermostatEnergy

      ! now we can estimate the new momentum for each thermostat

      allConverged = .false.
      converged = .false.

      nIterations = 1

      do while ( .not. allConverged )

         workvec1 = thermostatMomentumTemp * thermostatMomentumTemp * denominator
         workvec2 = potential + workvec1

         thermostatEnergyTemp = sum( workvec2 )

         deltaH = energy % kinetic + energy % potential + thermostatEnergyTemp - H0

         do n=1, nThermostats
            thermostatE(n) = sum( workvec1(1:n-1) )
         end do

         thermostatMomentumNew = thermostat % momentum +                            &
                  half_delta_t * thermostatProduct * ( two * energy % kinetic +     &
                  two * thermostatE + dpotential - deltaH ) / thermostat % position

         where( dabs( thermostatMomentumNew - thermostatMomentumTemp ) .lt. tolerance ) converged = .true.

         ! now check if we have converged

         if ( all( converged ) ) allConverged = .true.

         thermostatMomentumTemp = thermostatMomentumNew

         nIterations = nIterations + 1

         if ( nIterations > nIterationsMax ) then

            write(*,*) 'Maximum number of interations exceeded in '
            write(*,*) 'UpdateThermostatChainMomenta: stopping '

            stop

         end if

      end do

      thermostat % momentum = thermostatMomentumTemp

      ! and finally, calculate the thermostat energy consistent with the new momenta

      call calculateThermostatEnergy

      deallocate( denominator )
      deallocate( converged )
      deallocate( potential )
      deallocate( dpotential )
      deallocate( thermostatMomentumNew )
      deallocate( thermostatMomentumTemp )
      deallocate( thermostatE )
      deallocate( workvec1 )
      deallocate( workvec2 )

      if ( debug ) write(*,*) 'Exiting UpdateThermostatChainMomenta()'

end subroutine UpdateThermostatChainMomenta

end subroutine NosePoincareChain_A

