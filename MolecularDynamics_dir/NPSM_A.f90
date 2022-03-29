!*****************************************************************************
subroutine NPSM_A( state )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: MolecularDynamics (PRIVATE routine)
!
!  Created on: Wed  9 May 16:35:05 2018  by Eduardo R. Hernandez
!
!  Last Modified: Fri 10 May 11:01:07 2019
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

      integer :: i

      integer, parameter :: nIterationsMax = 20

      real (dp) :: E_kin_metric
      real (dp) :: factor
      real (dp) :: kinetic
      real (dp) :: stress_energy
      real (dp) :: thermostatProductOld, thermostatProduct
      real (dp) :: trace

      real (dp), parameter :: tolerance = 1.0e-7_dp

      real (dp), allocatable, dimension (:,:) :: lattice_velocity
      real (dp), dimension (3) :: vector
      real (dp), dimension (9) :: flat_stress, flat_metric

      type ( representation_type ) :: representation

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering NPSM_A'

      representation % position = .false.     ! in this routine we need a lattice rep.
      representation % force = .false.
      representation % momentum = .false.

!  the following array will be useful for updating ionic positions

      allocate( lattice_velocity( 3, nAtoms ) )

!  make sure coordinate, momenta and forces are in lattice representation

      call get( state, nThermostats = nThermostats,                              &
                position = position, force = force, momentum = momentum,         &
                thermostat = thermostat, barostat = barostat, energy = energy,   &
                H0 = H0, Cell = cell, reciprocalCell = reciprocalCell,           &
                stress = stress, representation = representation )

!  get the product of thermostat positions 

      thermostatProduct = product( thermostat % position )

!  first of all, advance momenta of atoms by half step

      momentum = momentum +                            &
          half_delta_t * thermostatProduct * force

!  now construct the lattice_velocity array

      do i = 1, nAtoms

          lattice_velocity(1,i) =                                                &
                    dot_product( reciprocalCell % metricTensor(1,1:3),         &
                                 momentum(1:3,i) ) / mass(i)

          lattice_velocity(2,i) =                                                &
                    dot_product( reciprocalCell % metricTensor(2,1:3),         &
                                 momentum(1:3,i) ) / mass(i)

          lattice_velocity(3,i) =                                                &
                    dot_product( reciprocalCell % metricTensor(3,1:3),         &
                                 momentum(1:3,i) ) / mass(i)

      end do

!  and calculate the kinetic energy and kinetic stress

      kinetic = zero
      kinetic_stress = zero

      do i=1, nAtoms

         vector = matmul( reciprocalCell % metricTensor, momentum(1:3,i) ) / mass(i)

         kinetic = kinetic +  dot_product( vector, momentum(1:3,i) )

         kinetic_stress(1,1) = kinetic_stress(1,1) + mass(i) * vector(1) * vector(1)
         kinetic_stress(1,2) = kinetic_stress(1,2) + mass(i) * vector(1) * vector(2)
         kinetic_stress(1,3) = kinetic_stress(1,3) + mass(i) * vector(1) * vector(3)

         kinetic_stress(2,1) = kinetic_stress(2,1) + mass(i) * vector(2) * vector(1)
         kinetic_stress(2,2) = kinetic_stress(2,2) + mass(i) * vector(2) * vector(2)
         kinetic_stress(2,3) = kinetic_stress(2,3) + mass(i) * vector(2) * vector(3)

         kinetic_stress(3,1) = kinetic_stress(3,1) + mass(i) * vector(3) * vector(1)
         kinetic_stress(3,2) = kinetic_stress(3,2) + mass(i) * vector(3) * vector(2)
         kinetic_stress(3,3) = kinetic_stress(3,3) + mass(i) * vector(3) * vector(3)

      end do

      factor = half / ( thermostatProduct * thermostatProduct )

      energy % kinetic = factor * kinetic

      kinetic_stress = factor * kinetic_stress

!  now update the metric tensor momenta to half-step; this needs to be 
!  done iteratively, so it is done in the following call

      call updateMetricTensorMomenta

!  if the shape of the cell is fixed (this may be needed e.g. for simulating
!  liquids at constant pressure), impose the necessary restrictions on the
!  metric tensor momenta

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

! make sure the barostat energy contains the updated metric tensor kinetic energy

      flat_stress = pack( simulation % stress % lattice, .true. )
      flat_metric = pack( cell % metricTensor, .true. )

      stress_energy = half * dot_product( flat_stress, flat_metric )

!     write(*,*) 'E_kin_metric = ', E_kin_metric
!     write(*,*) 'simulation % pressure = ', simulation % pressure
!     write(*,*) 'cell % volume = ', cell % volume
!     write(*,*) 'stress_energy = ', stress_energy

      energy % barostat = E_kin_metric + simulation % pressure * cell % volume + &
                stress_energy

! now advance the momenta of thermostats in the chain y half a step.
! due to the coupled-nature of the equations, this requires an iterative
! procedure, so we are going to do this in a subroutine

      call updateThermostatChainMomenta

! now advance the thermostat position to full-step; this is also going to require
! an iterative procedure, so we will pack it into the following routine

      call updateThermostatChainPositions

      thermostatProductOld = thermostatProduct
      thermostatProduct = product( thermostat % position )

!  now we can advance the components of the metric tensor to full-step; this, like the 
!  momenta, requires an iterative procedure that is implemented in the following call

      call updateMetricTensorComponents

!  update the metric kinetic energy; this changes because the volume changes when we update cell

      E_kin_metric = E_kin_metric * cell % volume * cell % volume 
 
!  update the cell according to the new metric tensor components

      call updateCell( cell, cell % metricTensor )

      E_kin_metric = E_kin_metric / ( cell % volume * cell % volume ) 

      flat_stress = pack( simulation % stress % lattice, .true. )
      flat_metric = pack( cell % metricTensor, .true. )

      stress_energy = half * dot_product( flat_stress, flat_metric )

      energy % barostat = E_kin_metric + simulation % pressure * cell % volume + &
                stress_energy

!  finally advance the lattice components of the atomic positions

      factor = half_delta_t * ( one / thermostatProductOld + one / thermostatProduct )

      position = position + factor * lattice_velocity

!  update kinetic enegy and kinetic stress (they need to be re-scaled with present thermostatPosition)

      energy % kinetic = kinetic / ( two * thermostatProduct * thermostatProduct )

      factor = thermostatProductOld * thermostatProductOld / thermostatProduct / thermostatProduct 

      kinetic_stress = factor * kinetic_stress

!  before returning, store updated variables in state

      call set( state, position = position, momentum = momentum, thermostat = thermostat,  &
                barostat = barostat, energy = energy, Cell = cell,                         &
                representation = representation )

!  before leaving, deallocate lattice_velocity storage

      deallocate( lattice_velocity )

      if ( debug ) write(*,*) 'Exiting NPSM_A'

contains     ! the following routines are specific to NPSM_A, so we keep them here rather than as module routines

!*****************************************************************************
subroutine updateMetricTensorComponents
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Created on: Wed  9 May 16:35:05 2018  by Eduardo R. Hernandez
!
!  Last Modified: Thu 10 May 09:44:51 2018
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
!     Purpose:          This subroutine updates the metric 
!                       tensor components to full-step. The equations that 
!                       appear in the Souza-Martins algorithm for constant
!                       pressure molecular dynamics can be integrated using
!                       the Generalised Leap-Frog (GLF) integration method.
!                       When this is done, the equations that appear for 
!                       the numerical integration of the components of the 
!                       metric tensor are implicit, and therefore
!                       we must use a Newton-Raphson procedure to solve them.
!                       This is what this subroutine does.
!
!*****************************************************************************
!  No implicits please!
   
   implicit none

!*****************************************************************************
!  shared variables

!*****************************************************************************
!  local variables

   logical :: converged

   integer :: i, n_newton_raphson_steps

   integer, parameter :: n_newton_raphson_steps_max = 100

   real (dp) :: det_metric, det_metric_tmp, factor_1, factor_2
   real (dp) :: cofactor_11, cofactor_12, cofactor_13
   real (dp) :: cofactor_22, cofactor_23, cofactor_33

   real (dp), parameter :: tolerance = 1.0e-7_dp

   real (dp) :: difference(9)
   real (dp) :: gpi(3,3)
   real (dp) :: pig(3,3)
   real (dp) :: gpig(3,3)
   real (dp) :: gpig_old(3,3)
   real (dp) :: metricTensor_new(3,3)
   real (dp) :: metricTensor_tmp(3,3)
   real (dp) :: recipr_metricTensor_tmp(3,3)

!*****************************************************************************
!  start of subroutine

   if ( debug ) write(*,*) 'Entering updateMetricTensorComponents()'

! extract from state 

! the first thing is to generate the half-step velocities of the metric 
! in real space, as up to now we have been using them in reciprocal space

   gpi = matmul( cell % metricTensor, barostat % momentum )
   gpig = matmul( gpi, cell % metricTensor )

! initialise the tmp momenta making them equal to the old

   metricTensor_tmp = cell % metricTensor
   gpig_old = gpig

   det_metric = cell % volume * cell % volume

   n_newton_raphson_steps = 0

   do 

! first calculate the reciprocal metric tensor, the determinant and the
! gpig, gpi and pig matrices based on the trial metric tensor

      cofactor_11 = metricTensor_tmp(2,2) * metricTensor_tmp(3,3) -        &
                    metricTensor_tmp(2,3) * metricTensor_tmp(2,3)
      cofactor_12 = metricTensor_tmp(1,3) * metricTensor_tmp(2,3) -        &
                    metricTensor_tmp(1,2) * metricTensor_tmp(3,3)
      cofactor_13 = metricTensor_tmp(1,2) * metricTensor_tmp(2,3) -        &
                    metricTensor_tmp(1,3) * metricTensor_tmp(2,2)

      cofactor_22 = metricTensor_tmp(1,1) * metricTensor_tmp(3,3) -        &
                    metricTensor_tmp(1,3) * metricTensor_tmp(1,3)
      cofactor_23 = metricTensor_tmp(1,2) * metricTensor_tmp(1,3) -        &
                    metricTensor_tmp(1,1) * metricTensor_tmp(2,3)

      cofactor_33 = metricTensor_tmp(1,1) * metricTensor_tmp(2,2) -        &
                    metricTensor_tmp(1,2) * metricTensor_tmp(1,2)

      det_metric_tmp = metricTensor_tmp(1,1) * cofactor_11 +                &
                       metricTensor_tmp(1,2) * cofactor_12 +                &
                       metricTensor_tmp(1,3) * cofactor_13

      recipr_metricTensor_tmp(1,1) = cofactor_11 / det_metric_tmp
      recipr_metricTensor_tmp(1,2) = cofactor_12 / det_metric_tmp
      recipr_metricTensor_tmp(1,3) = cofactor_13 / det_metric_tmp

      recipr_metricTensor_tmp(2,1) = cofactor_12 / det_metric_tmp
      recipr_metricTensor_tmp(2,2) = cofactor_22 / det_metric_tmp
      recipr_metricTensor_tmp(2,3) = cofactor_23 / det_metric_tmp
 
      recipr_metricTensor_tmp(3,1) = cofactor_13 / det_metric_tmp
      recipr_metricTensor_tmp(3,2) = cofactor_23 / det_metric_tmp
      recipr_metricTensor_tmp(3,3) = cofactor_33 / det_metric_tmp

      gpi = matmul( metricTensor_tmp, barostat % momentum ) 
      gpig = matmul( gpi, metricTensor_tmp )
      pig = transpose( gpi )

! we then predict the tmp metricTensor

      factor_1 = half_delta_t * thermostatProductOld / ( det_metric * barostat % mass )
      factor_2 = half_delta_t * thermostatProduct /                           &
                                           ( det_metric_tmp * barostat % mass )

      metricTensor_new = cell% metricTensor + factor_1 * gpig_old + factor_2 * gpig

! lets calculate the differences between predicted and actual momenta
! remember: metric tensor and its momentum componets are symmetric

      difference(1) = metricTensor_new(1,1) -                               &
                      metricTensor_tmp(1,1)
      difference(2) = metricTensor_new(1,2) -                               &
                      metricTensor_tmp(1,2)
      difference(3) = metricTensor_new(1,3) -                               &
                      metricTensor_tmp(1,3)
      difference(4) = metricTensor_new(2,1) -                               &
                      metricTensor_tmp(2,1)
      difference(5) = metricTensor_new(2,2) -                               &
                      metricTensor_tmp(2,2)
      difference(6) = metricTensor_new(2,3) -                               &
                      metricTensor_tmp(2,3)
      difference(7) = metricTensor_new(3,1) -                               &
                      metricTensor_tmp(3,1)
      difference(8) = metricTensor_new(3,2) -                               &
                      metricTensor_tmp(3,2)
      difference(9) = metricTensor_new(3,3) -                               &
                      metricTensor_tmp(3,3)

      converged = .true.

      do i=1, 9

         if ( dabs( difference(i) ) .gt. tolerance ) converged = .false.

      end do
 
      if ( .not. converged ) then

         metricTensor_tmp = metricTensor_new

         metricTensor_tmp = half * ( metricTensor_tmp +                  &
                           transpose( metricTensor_tmp ) )
         
      else

         exit 

      end if

      n_newton_raphson_steps = n_newton_raphson_steps + 1

      if ( n_newton_raphson_steps > n_newton_raphson_steps_max ) then

         write(*,*) 'n_newton_raphson_max exceeded A: stopping'
         stop

      end if

   end do

   cell % metricTensor = metricTensor_tmp

   if ( debug ) write(*,*) 'Exiting updateMetricTensorComponents()'

end subroutine updateMetricTensorComponents

!*****************************************************************************
subroutine updateMetricTensorMomenta
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Created on: Wed  9 May 16:35:05 2018  by Eduardo R. Hernandez
!
!  Last Modified: Thu 10 May 09:44:51 2018
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
!                       This subroutine updates the momenta of the metric 
!                       tensor components to half-step. The equations that 
!                       appear in the Souza-Martins algorithm for constant
!                       pressure molecular dynamics can be integrated using
!                       the Generalised Leap-Frog (GLF) integration method.
!                       When this is done, the equations that appear for 
!                       the numerical integration of the momenta of the 
!                       metric tensor components are implicit, and therefore
!                       we must use a Newton-Raphson procedure to solve them.
!                       This is what this subroutine does.
!
!                       IMPORTANT: this subroutine updates the components
!                                  of the momenta for the reciprocal metric
!                                  tensor
!*****************************************************************************
!  No implicits please!
   
   implicit none

!*****************************************************************************
!  shared variables

!*****************************************************************************
!  local variables

   logical :: converged

   integer :: i, n_newton_raphson_steps

   integer, parameter :: n_newton_raphson_steps_max = 100

   real (dp) :: det_metric, factor

   real (dp), parameter :: tolerance = 1.0e-7_dp

   real (dp) :: difference(9)
   real (dp) :: matrix_1(3,3)
   real (dp) :: matrix_2(3,3)
   real (dp) :: matrix_3(3,3)
   real (dp) :: momentum_new(3,3)
   real (dp) :: momentum_tmp(3,3)

!*****************************************************************************
!  start of subroutine

   if ( debug ) write(*,*) 'Entering updateMetricTensorMomenta()'

   det_metric = cell % volume * cell % volume

! initialise the tmp momenta making them equal to the old

   momentum_tmp = barostat % momentum

!  write(*,'(3f15.7)') momentum_tmp(1:3,1)
!  write(*,'(3f15.7)') momentum_tmp(1:3,2)
!  write(*,'(3f15.7)') momentum_tmp(1:3,3)

   n_newton_raphson_steps = 0

   do 

! we first predict the tmp momenta

      factor = one / ( barostat % mass * det_metric )
 
!     write(*,*) 'factor = ', factor, barostat % mass, det_metric

      matrix_1 = factor * matmul( momentum_tmp, cell % metricTensor )

      matrix_2 = matmul( matrix_1, momentum_tmp )

      matrix_3 = transpose( matrix_1 )

      factor = half * barostat % mass * det_metric

      E_kin_metric = factor * dot_product( pack( matrix_1, .true. ), pack( matrix_3, .true. ) )

!     write(*,*) n_newton_raphson_steps, E_kin_metric

      momentum_new = barostat % momentum -     &
                 half_delta_t * thermostatProduct * (               &
                stress % lattice - kinetic_stress + matrix_2               & 
       + ( half * simulation % pressure * cell % volume - E_kin_metric ) * &
         reciprocalCell % metricTensor +                           & 
                              half * simulation % stress % lattice )

      momentum_new = half * ( momentum_new + transpose( momentum_new ) )

! lets calculate the differences between predicted and actual momenta
! remember: metric tensor and its momentum componets are symmetric

      difference(1) = momentum_new(1,1) -                       &
                      momentum_tmp(1,1)
      difference(2) = momentum_new(1,2) -                       &
                      momentum_tmp(1,2)
      difference(3) = momentum_new(1,3) -                       &
                      momentum_tmp(1,3)
      difference(4) = momentum_new(2,1) -                       &
                      momentum_tmp(2,1)
      difference(5) = momentum_new(2,2) -                       &
                      momentum_tmp(2,2)
      difference(6) = momentum_new(2,3) -                       &
                      momentum_tmp(2,3)
      difference(7) = momentum_new(3,1) -                       &
                      momentum_tmp(3,1)
      difference(8) = momentum_new(3,2) -                       &
                      momentum_tmp(3,2)
      difference(9) = momentum_new(3,3) -                       &
                      momentum_tmp(3,3)

! the correspondence of indices is: 1,1-> 1, 1,2-> 2, 1,3-> 3, 
!                                   2,1-> 4, 2,2-> 5, 2,3-> 6,
!                                   3,1-> 7, 3,2-> 8, 3,3-> 6

      converged = .true.

      do i=1, 9

         if ( dabs( difference(i) ) .gt. tolerance ) converged = .false.

      end do

      if ( .not. converged ) then

         momentum_tmp = momentum_new

         momentum_tmp = half * ( momentum_tmp + transpose( momentum_tmp ) )
         
      else

         exit 

      end if

      n_newton_raphson_steps = n_newton_raphson_steps + 1

      if ( n_newton_raphson_steps > n_newton_raphson_steps_max ) then

         write(6,1) momentum_new(1,1:3)
         write(6,1) momentum_new(2,1:3)
         write(6,1) momentum_new(3,1:3)

         write(*,*) 'n_newton_raphson_max exceeded B: stopping'
         stop

      end if

   end do

   factor = one / ( barostat % mass * det_metric )

   matrix_1 = factor * matmul( momentum_tmp, cell % metricTensor )

   matrix_2 = matmul( matrix_1, momentum_tmp )

   matrix_3 = transpose( matrix_1 )

   factor = half * barostat % mass * det_metric

   E_kin_metric = factor * dot_product( pack( matrix_1, .true. ), pack( matrix_3, .true. ) )

!  write(*,*) 'in update E_kin_metric = ', E_kin_metric

   barostat % momentum = momentum_tmp

 1   format(3f15.7)

   if ( debug ) write(*,*) 'Exiting updateMetricTensorMomenta()'

end subroutine updateMetricTensorMomenta

!*****************************************************************************
subroutine updateThermostatChainPositions
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

      if ( debug ) write(*,*) 'Entering updateThermostatChainPositions()'

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
            write(*,*) 'updateThermostatChainPositions: stopping '

            stop

         end if

      end do

      thermostat % position = thermostatPositionTemp

      deallocate( thermostatPositionNew )
      deallocate( thermostatPositionTemp )
      deallocate( denominator )
      deallocate( denominatorNew )
      deallocate( difference )

      if ( debug ) write(*,*) 'Exiting updateThermostatChainPositions()'

end subroutine updateThermostatChainPositions

!*****************************************************************************
subroutine updateThermostatChainMomenta
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

      if ( debug ) write(*,*) 'Entering updateThermostatChainMomenta()'

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

!     write(*,*) 'kinetic ', energy % kinetic
!     write(*,*) 'potential ', energy % potential
!     write(*,*) 'barostat ', energy % barostat

      do while ( .not. allConverged )

         workvec1 = thermostatMomentumTemp * thermostatMomentumTemp * denominator
         workvec2 = potential + workvec1

         thermostatEnergyTemp = sum( workvec2 )

!        write(*,*) nIterations, thermostatEnergyTemp

         deltaH = energy % kinetic + energy % potential + energy % barostat + thermostatEnergyTemp - H0

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
            write(*,*) 'updateThermostatChainMomenta: stopping '

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

      if ( debug ) write(*,*) 'Exiting updateThermostatChainMomenta()'

end subroutine updateThermostatChainMomenta

end subroutine NPSM_A

