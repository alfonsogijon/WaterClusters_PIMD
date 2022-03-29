!*****************************************************************************
subroutine NPSM_B( state )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: MolecularDynamics (PRIVATE routine)
!
!  Created on: Mon 14 May 09:37:36 2018  by Eduardo R. Hernandez
!
!  Last Modified: Thu 17 Jan 19:11:17 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! Implements the updating of ionic, barostat and, if required, thermostat
! momenta to half step, and the corresponding positions to full step according
! to the Generalised Leap-Frog integration algorithm (for details see
! Hernandez, J. Chem. Phys. vol. 115, p. 10282 (2001)) as a first step to the
! numerical integration of the equations of motion in the NPH (constant pressure
! constant enthalpy) or NPT (constant pressure, constant temperature) ensembles.
! As usual, following a call to this routine, forces and stresses at the 
! new positions must be evaluated by the client program, after which NPSM_B
! can be called to complete the step by updating the momenta to full step.
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
!  local variables

      integer :: i, n

      real (dp) :: const
      real (dp) :: deltaH
      real (dp) :: E_kin_metric
      real (dp) :: factor
      real (dp) :: kbT
      real (dp) :: prod
      real (dp) :: stress_energy
      real (dp) :: thermostatProduct
      real (dp) :: trace

      real (dp), dimension (9) :: flat_metric, flat_stress
      real (dp), dimension (3,3) :: matrix_a, matrix_b, matrix_c
      real (dp), dimension (3) :: vector

      type ( representation_type ) :: representation

      real (dp), dimension (:), allocatable :: denominator
      real (dp), dimension (:), allocatable :: potential
      real (dp), dimension (:), allocatable :: dpotential
      real (dp), dimension (:), allocatable :: thermostatE
      real (dp), dimension (:), allocatable :: workvec1
      real (dp), dimension (:), allocatable :: workvec2

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering NPSM_B'

      representation % position = .false.     ! in this routine we need a lattice rep.
      representation % force = .false.
      representation % momentum = .false.

!  make sure coordinate, momenta and forces are in lattice representation

      call get( state, nThermostats = nThermostats,                              &
                position = position, force = force, momentum = momentum,         &
                thermostat = thermostat, barostat = barostat, energy = energy,   &
                H0 = H0, Cell = cell, reciprocalCell = reciprocalCell,           &
                stress = stress, representation = representation )

!  update the thermostat momentum

      thermostatProduct = product( thermostat % position )

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

      call calculateThermostatEnergy

      deltaH = energy % kinetic + energy % potential + energy % barostat + energy % thermostat - H0

      do n=1, nThermostats
         thermostatE(n) = sum( workvec1(1:n-1) )
      end do

      thermostat % momentum = thermostat % momentum +                                &
            half_delta_t * thermostatProduct * ( two * energy % kinetic +     &
                   two * thermostatE + dpotential - deltaH ) / thermostat % position

! now update the metric tensor momenta

      factor = one /                                                         & 
         ( barostat % mass * cell % volume * cell % volume ) 

      matrix_a = factor * matmul( barostat % momentum, cell % metricTensor )

      matrix_b = matmul( matrix_a, barostat % momentum )

      matrix_c = transpose( matrix_a )

      factor = half / factor

      E_kin_metric = factor * dot_product( pack( matrix_a, .true. ), pack( matrix_c, .true. ))

      barostat % momentum =  barostat % momentum -           &
         half_delta_t * thermostatProduct * ( stress % lattice -  &
                    kinetic_stress + matrix_b +                              & 
          ( half * simulation % pressure * cell % volume -          & 
            E_kin_metric ) * reciprocalCell % metricTensor +       &
                    half * simulation % stress % lattice ) 

! again, if this is an isoshape calculation, we must impose restrictions on 
! the metric tensor momenta

      if ( simulation % isoshape ) then

         barostat % momentum(1,2) = zero
         barostat % momentum(1,3) = zero
         barostat % momentum(2,3) = zero
         barostat % momentum(2,1) = zero
         barostat % momentum(3,1) = zero
         barostat % momentum(3,2) = zero

         trace = barostat % momentum(1,1) +             &
                 barostat % momentum(2,2) +             &
                 barostat % momentum(3,3) 

         barostat % momentum(1,1) = trace / three
         barostat % momentum(2,2) = trace / three
         barostat % momentum(3,3) = trace / three

      end if

! now update the metric tensor kinetic energy

      factor = half / ( barostat % mass * cell % volume * cell % volume ) 

      matrix_a = matmul( barostat % momentum, cell % metricTensor )

      matrix_c = transpose( matrix_a )

      E_kin_metric = factor * dot_product( pack( matrix_a, .true. ), pack( matrix_c, .true. ))

      flat_stress = pack( simulation % stress % lattice, .true. )
      flat_metric = pack( cell % metricTensor, .true. )

      stress_energy = half * dot_product( flat_stress, flat_metric )

! and update the atomic momenta

      momentum = momentum + half_delta_t * thermostatProduct * force

      energy % kinetic = 0.0_dp

      const = half / ( thermostatProduct * thermostatProduct )

!  in the above, th_position = 1 if no thermostat is used

      do i = 1, nAtoms

         factor = const / mass(i) 

         vector = matmul( reciprocalCell % metricTensor, momentum(1:3,i) )

         energy % kinetic = energy % kinetic + factor *                    &
          dot_product( vector, momentum(1:3,i) )

      end do

      call calculateThermostatEnergy

      energy % barostat = E_kin_metric + simulation % pressure * cell % volume + &
                stress_energy

      energy % constantOfMotion = thermostatProduct * (                       &
         energy % kinetic + energy % potential +                              &
         energy % thermostat + energy % barostat - H0 )

! before leaving, store the updated variables in state

      call set( state, momentum = momentum, thermostat = thermostat,  &
                barostat = barostat, energy = energy, representation = representation )

! and deallocate temporary arrays 

      deallocate( denominator )
      deallocate( potential )
      deallocate( dpotential )
      deallocate( thermostatE )
      deallocate( workvec1 )
      deallocate( workvec2 )

      if ( debug ) write(*,*) 'Exiting NPSM_B'

end subroutine NPSM_B

