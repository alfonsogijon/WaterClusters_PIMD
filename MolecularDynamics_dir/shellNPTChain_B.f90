!*****************************************************************************
subroutine shellNPTChain_B( state )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: MolecularDynamics (PRIVATE routine)
!
!  Created on: Fri 22 Mar 16:05:36 2018  by Eduardo R. Hernandez
!
!s  Last Modified: Thu  4 Jul 17:23:44 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! Implements the updating of ionic and, if required, thermostat
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

      integer :: i, j, n

      real (dp) :: const
      real (dp) :: deltaH
      real (dp) :: factor
      real (dp) :: kbT
      real (dp) :: prod
      real (dp) :: shellKinetic
      real (dp) :: stress_energy
      real (dp) :: thermostatProduct

      real (dp), dimension (3) :: vector
      real (dp), dimension (3) :: p, r, s, t

      type ( representation_type ) :: representation

      real (dp), dimension (:), allocatable :: denominator
      real (dp), dimension (:), allocatable :: potential
      real (dp), dimension (:), allocatable :: dpotential
      real (dp), dimension (:), allocatable :: thermostatE
      real (dp), dimension (:), allocatable :: workvec1
      real (dp), dimension (:), allocatable :: workvec2
      real (dp), dimension (3,3) :: CartesianForce
      real (dp), dimension (3,3) :: CartesianMomentum
      real (dp), dimension (3,3) :: cellForce
      real (dp), dimension (9) :: flatStress, flatMetric
      real (dp), dimension (9) :: flatCellForce

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering shellNPTChain_B'

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

!  we need to calculate the thermostat energy with current momenta and positions

      workvec1 = thermostat % momentum * thermostat % momentum * denominator

      call calculateThermostatEnergy

      deltaH = energy % kinetic + energy % potential + energy % barostat + energy % thermostat - H0

      do n=1, nThermostats
         thermostatE(n) = sum( workvec1(1:n-1) )
      end do

      thermostat % momentum = thermostat % momentum +                                &
            half_delta_t * thermostatProduct * ( two * energy % kinetic +     &
                   two * thermostatE + dpotential - deltaH ) / thermostat % position

!  now advance the shell-atoms momenta

! convert the shell-atom momenta to Cartesian; likewise with the forces

      ! this transforms m_i G_ab \dot{s_ib} first to m_i \dot{s_ia}
      ! and then to m_i H_ab \dot_{s_ib} = m_i v_ia

      r = (/ cell % a(1), cell % b(1), cell % c(1) /)
      s = (/ cell % a(2), cell % b(2), cell % c(2) /)
      t = (/ cell % a(3), cell % b(3), cell % c(3) /)

      j = 1

      do i=nAtoms-2, nAtoms

         p = matmul( reciprocalCell % metricTensor, momentum(1:3,i) )

         CartesianMomentum(1,j) = dot_product( r, p )
         CartesianMomentum(2,j) = dot_product( s, p )
         CartesianMomentum(3,j) = dot_product( t, p )

         j = j + 1

      end do

! now convert the forces

      r = (/ reciprocalCell % a(1), reciprocalCell % a(2), reciprocalCell % a(3) /)
      s = (/ reciprocalCell % b(1), reciprocalCell % b(2), reciprocalCell % b(3) /)
      t = (/ reciprocalCell % c(1), reciprocalCell % c(2), reciprocalCell % c(3) /)

      j = 1

      do i=nAtoms-2, nAtoms

         CartesianForce(1,j) = dot_product( r, force(1:3,i) )
         CartesianForce(2,j) = dot_product( s, force(1:3,i) )
         CartesianForce(3,j) = dot_product( t, force(1:3,i) )

         j = j + 1

      end do

! now calculate the contributions to the shell-atom forces coming from the cell

      cellForce = -( stress % lattice + kinetic_stress +                    &
                       half * simulation % pressure * cell % volume *       &
                       reciprocalCell % metricTensor +                      &
                       half * simulation % stress % lattice )

      flatCellForce = pack( cellForce, .true. )

! with this we have everything we need to advance the shell-atom momenta to half step

      do i=1, 3

         CartesianMomentum(1,i) = CartesianMomentum(1,i) +                  &
                half_delta_t * thermostatProduct * ( CartesianForce(1,i) +  &
                      dot_product( flatCellForce, dGdr(:,1,i) ) )

         CartesianMomentum(2,i) = CartesianMomentum(2,i) +                  &
                half_delta_t * thermostatProduct * ( CartesianForce(2,i) +  &
                      dot_product( flatCellForce, dGdr(:,2,i) ) )

         CartesianMomentum(3,i) = CartesianMomentum(3,i) +                  &
                half_delta_t * thermostatProduct * ( CartesianForce(3,i) +  &
                      dot_product( flatCellForce, dGdr(:,3,i) ) )

      end do

! finally, transform the shell-atom momenta back to lattice reppresentation, 
! and in passing, calculate the shell-atom contribution to the kinetic energy

      shellKinetic = zero

      j = 1

      do i=nAtoms-2, nAtoms

         shellKinetic = shellKinetic +                                      &
           dot_product( CartesianMomentum(:,j), CartesianMomentum(:,j) ) /  &
                          mass(i)

         p(1) = dot_product( reciprocalCell % a,         &
                                       CartesianMomentum(1:3,j) )
         p(2) = dot_product( reciprocalCell % b,         &
                                       CartesianMomentum(1:3,j) )
         p(3) = dot_product( reciprocalCell % c,         &
                                       CartesianMomentum(1:3,j) )

         momentum(:,i) = matmul( cell % metricTensor, p )

         j = j + 1

      end do

      shellKinetic = half * shellKinetic / thermostatProduct / thermostatProduct

!  and update the atomic momenta

      do i=1, nAtoms-3

         momentum(:,i) = momentum(:,i) + half_delta_t * thermostatProduct * force(:,i)

      end do

!  finally, calculate the total kinetic energy

      energy % kinetic = zero

      const = half / ( thermostatProduct * thermostatProduct )

!  in the above, th_position = 1 if no thermostat is used

      do i = 1, nAtoms - 3

         factor = const / mass(i) 

         vector = matmul( reciprocalCell % metricTensor, momentum(1:3,i) )

         energy % kinetic = energy % kinetic + factor *                    &
              dot_product( vector, momentum(1:3,i) )

      end do

      energy % kinetic = energy % kinetic + shellKinetic

      call calculateThermostatEnergy

      flatStress = pack( simulation % stress % lattice, .true. )
      flatMetric = pack( cell % metricTensor, .true. )

      stress_energy = half * dot_product( flatStress, flatMetric )

      energy % barostat = simulation % pressure * cell % volume + &
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

      if ( debug ) write(*,*) 'Exiting shellNPTChain_B'

end subroutine shellNPTChain_B

