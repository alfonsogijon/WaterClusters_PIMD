!*****************************************************************************
subroutine NosePoincareChain_B( state )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: MolecularDynamics (PRIVATE routine)
!
!  Created on: Tue 24 Jul 16:28:32 2018  by Eduardo R. Hernandez
!
!  Last Modified: Wed 27 Mar 10:32:36 2019
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

      integer :: n
      integer :: nThermostats

      real (dp) :: const
      real (dp) :: deltaH
      real (dp) :: factor
      real (dp) :: kbT
      real (dp) :: prod
      real (dp) :: thermostatProduct

      real (dp), dimension (:), allocatable :: denominator
      real (dp), dimension (:), allocatable :: potential
      real (dp), dimension (:), allocatable :: dpotential
      real (dp), dimension (:), allocatable :: thermostatE
      real (dp), dimension (:), allocatable :: workvec1
      real (dp), dimension (:), allocatable :: workvec2

      type ( representation_type ) :: representation 

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering NosePoincareChain_B()'

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

! first update the thermostat momenta to full step

      allocate( denominator( nThermostats ) )
      allocate( potential( nThermostats ) )
      allocate( dpotential( nThermostats ) )
      allocate( thermostatE( nThermostats ) )
      allocate( workvec1( nThermostats ) )
      allocate( workvec2( nThermostats ) )

      kbT = boltzmann_k * simulation % temperature

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

      ! we need to calculate the thermostat energy with current momenta and positions

      workvec1 = thermostat % momentum * thermostat % momentum * denominator
      ! workvec2 = potential + workvec1

      ! energy % thermostat = sum( workvec2 )

      call calculateThermostatEnergy

      deltaH = energy % kinetic + energy % potential + energy % thermostat - H0

      do n=1, nThermostats
         thermostatE(n) = sum( workvec1(1:n-1) )
      end do

      thermostat % momentum = thermostat % momentum +                                &
            half_delta_t * thermostatProduct * ( two * energy % kinetic +     &
                   two * thermostatE + dpotential - deltaH ) / thermostat % position

      ! and finally the atomic momenta 

      momentum = momentum +                            &
            half_delta_t * thermostatProduct * force

      ! now calculate kinetic energy with new momenta at full-step

      energy % kinetic = 0.0_dp

      const = half / ( thermostatProduct * thermostatProduct )

      do n = 1, nAtoms

         factor = const / mass(n) 

         energy % kinetic = energy % kinetic + factor *   &
          dot_product( momentum(1:3,n), momentum(1:3,n) )

      end do

! calculate the thermostat chain energy

      call calculateThermostatEnergy

! and the constant of motion

      energy % constantOfMotion = thermostatProduct * ( energy % kinetic +         &
                energy % potential + energy % thermostat - H0 )

! finally, store the new coordinates and half-step momenta in state ...

      call set( state, momentum = momentum, thermostat = thermostat, energy = energy )

! and deallocate temporary arrays

      deallocate( denominator )
      deallocate( potential )
      deallocate( dpotential )
      deallocate( thermostatE )
      deallocate( workvec1 )
      deallocate( workvec2 )

      if ( debug ) write(*,*) 'Exiting NosePoincareChain_B()'

end subroutine NosePoincareChain_B

