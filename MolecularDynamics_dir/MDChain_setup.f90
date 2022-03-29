!*****************************************************************************
subroutine MDChain_setup( state_chain, dyntype, timestep, temperature, &
     pressure, stress, fixCellShape, friction )                    
!*****************************************************************************
!
!  Project: MolecularDynamics 
!
!  Module: MolecularDynamics (PUBLIC routine)
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
! This subroutine must be called once, by the MD-driver routine, previous to 
! entering the MD loop (loop over time-steps). The purpose of this routine is
! to specify a number of system parameters (number of atoms, atomic masses,
! cell vectors of the simulation, cell, etc), and so on. 
!
! INPUT VARIABLES
! 
!!!!!!!!!!!!!!! MANDATORY ARGUMENTS !!!!!!!!!!!!!!!
!
! state (StateChain_type) :: the structure containing the system information
! 
! start (bool) :: indicates if the run is a start (initial) run or a 
!                 restart (continuation) run; this info is needed to properly
!                 set up the run
!
! dyntype (int) :: the type of dynamics, NVE (0), NVT (1), NPH (2), NPT (3)
!
! timestep (real): the time step, in fs (no choice of units here)
!
! temperature (real): the target temperature of the simulation in Kelvin
!
!!!!!!!!!!!!!!! OPTIONAL ARGUMENTS !!!!!!!!!!!!!!!
!
! pressure (real): the external pressure imposed on the system in the chosen 
!                  input units of pressure
!
! stress (real 3x3): the external imposed stress (non-hydrostatic part)
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  arguments

      type ( StateChain_type ), intent (INOUT) :: state_chain

      integer, intent (IN) :: dyntype

      real (dp), intent (IN) :: timestep

      real (dp), intent (IN) :: temperature
      real (dp), optional :: pressure
      real (dp), dimension (3,3), optional :: stress

      logical, intent (IN), optional :: fixCellShape

      real (dp), intent(IN), optional :: friction

!*****************************************************************************
!  local variables

      integer :: i, j
      integer :: n

      integer :: iclock, isize
      integer, dimension(:), allocatable :: iseed

      real(dp) :: kinetic_iReplica
      
      real (dp) :: stress_energy

      real (dp), dimension (9) :: flat_stress, flat_metric

      real (dp) :: beta

      real (dp) :: aux, dseed

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering MDChain_setup()'

      simulation % dyntype = dyntype

      ! define the time-step in internal units

      delta_t = timeConversionFactor * timestep * 1.0e-15_dp
      half_delta_t = half * delta_t
      delta_t2 = delta_t * delta_t

      ! set in the external thermodynamic parameters of the simulation

      simulation % temperature = temperature

      if ( present( pressure ) ) then

         write(*,*) 'Pressure conversion factor = ', pressureConversionFactor

         simulation % pressure = pressureConversionFactor * pressure

      else 

         simulation % pressure = zero

      end if

      if ( present( stress ) ) then

         simulation % stress % cartesian = pressureConversionFactor * stress

      else

         simulation % stress % cartesian = zero 
         simulation % stress % lattice = zero

      end if

      ! check if the cell shape is to be kept fixed (only applicable in NPH/NPT simulations)

      if ( present( fixCellShape ) ) then

         simulation % isoshape = fixCellShape 

      else

         simulation % isoshape = .false.

      end if

      ! now set local copies of number of atoms, and allocate appropriate arrays
      ! note that it is assumed that nAtoms stays constant during the course of
      ! a simulation; likewise happens with masses of atoms and extended variables

      call getStateChain( state_chain, nReplicas = nReplicas, nAtoms = nAtoms, &
           nDegreesOfFreedom = nDegreesOfFreedom ) 

      ! now we can allocate local arrays
    
      allocate( mass( nAtoms ) )

      allocate( state_replica(nReplicas) )
      allocate( position_chain( 3, nAtoms, nReplicas ) )
      allocate( force_chain( 3, nAtoms, nReplicas ) )
      allocate( momentum_chain( 3, nAtoms, nReplicas ) )

      allocate( c2( nAtoms) )

      ! and calculate the kinetic energy

      call getStateChain( state_chain, energy = energy )
      
      energy % kinetic = zero

      do i = 1, nReplicas

         kinetic_iReplica = zero
         
         ! store the current values

         call getStateChain( state_chain, iReplica = i, state = state_replica(i) )

         call get( state_replica(i), position = position_chain(:,:,i), &
              momentum = momentum_chain(:,:,i), mass = mass, &
              Cell = cell, reciprocalCell = reciprocalCell )

         do j=1, nAtoms

            kinetic_iReplica = kinetic_iReplica + half * &
                 dot_product( momentum_chain(1:3,j,i), momentum_chain(1:3,j,i) ) / mass(j)

         end do

         energy % kinetic = energy % kinetic + kinetic_iReplica
         
      end do

      ! now, if the kinetic energy is small (a signal that momenta have not been 
      ! read from elsewhere) then generate atomic momenta here; the criterion is
      ! that if the temperature deduced from the current kinetic energy is less
      ! than 10% of the target temperature, and this is a start calculation, we
      ! generate momenta here

      ! if dyntype > 0 we might have one or more thermostats

      if ( simulation % dyntype > 0 ) then

         call getStateChain( state_chain, nThermostats = nThermostats )

         if ( nThermostats  > 0 ) then

            allocate( thermostat( nThermostats ) )

            call getStateChain( state_chain, thermostat = thermostat )

            ! calculate the initial thermostat energy

            call calculateThermostatEnergy

         end if

      end if

      if ( ( simulation % dyntype == 2 ) .or. ( simulation % dyntype == 3 ) ) then   ! we have a barostat

         flat_stress = pack( simulation % stress % lattice, .true. )
         flat_metric = pack( cell % metricTensor, .true. )

         stress_energy = half * dot_product( flat_stress, flat_metric )

         energy % barostat = simulation % pressure * cell % volume + stress_energy

      end if

      if ( simulation % dyntype == 4 ) then ! Bussi-Parrinello thermostat

         beta = one / (boltzmann_k*temperature)

         gamma = friction
         
         c1 = dexp( -gamma * half_delta_t )
         c2(:) = dsqrt( (one-c1*c1) * mass(:) / beta )
         
      end if

      call setStateChain( state_chain, energy = energy )

      ! initialize random number generator

      call random_seed(SIZE=isize)
      allocate(iseed(isize))
      call system_clock(COUNT=iclock)
      iseed = iclock + 37 * (/ (i - 1, i = 1, isize) /)
      call random_seed(PUT=iseed)
      
      call random_number(dseed)
      !dseed = 954132
      aux = randini(dseed)

      if ( debug ) write(*,*) 'Exiting MDChain_setup()'

end subroutine MDChain_setup
