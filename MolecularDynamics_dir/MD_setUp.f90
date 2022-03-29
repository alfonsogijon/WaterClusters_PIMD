!*****************************************************************************
subroutine MD_setUp( state, dyntype, timestep, temperature,                  &
                      pressure, stress, fixCellShape )                                         
!*****************************************************************************
!
!  Project: MolecularDynamics 
!
!  Module: MolecularDynamics (PUBLIC routine)
!
!  Created on: Fri 16 Mar 14:35:15 2018  by Eduardo R. Hernandez
!
!  Last Modified: Mon 15 Apr 17:00:31 2019
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
! to specify a number of system 
! parameters (number of atoms, atomic masses, cell vectors of the simulation
! cell, etc), and so on. A detailed list of input variables with appropriate
! explanations is given below. Unless otherwise specified, a given input
! variable is mandatory.
!
! INPUT VARIABLES
! 
!!!!!!!!!!!!!!! MANDATORY ARGUMENTS !!!!!!!!!!!!!!!
!
! state (state_type) :: the structure containing the system information
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
! fixCellShape (bool): logical flag indicating if the cell shape is to be held
!                      fixed during constant pressure simulations 
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  arguments

      type ( state_type ), intent (INOUT) :: state

      integer, intent (IN) :: dyntype

      real (dp), intent (IN) :: timestep

      real (dp) :: temperature
      real (dp), optional :: pressure
      real (dp), dimension (3,3), optional :: stress

      logical, intent (IN), optional :: fixCellShape

!*****************************************************************************
!  local variables

      integer :: i

      real (dp) :: stress_energy

      real (dp), dimension (9) :: flat_stress, flat_metric

      type ( representation_type ) :: representation

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering MD_setUp()'

      simulation % dyntype = dyntype

      ! define the time-step in internal units

      delta_t = timeConversionFactor * timestep * 1.0e-15_dp
      half_delta_t = half * delta_t

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
      ! a simulation; likewise happes with masses of atoms and extended variables

      call get( state, nAtoms = nAtoms, nDegreesOfFreedom = nDegreesOfFreedom ) 

      ! now we can allocate local arrays
    
      allocate( position( 3, nAtoms ) )
      allocate( force( 3, nAtoms ) )
      allocate( momentum( 3, nAtoms ) )
      allocate( mass( nAtoms ) )

      ! set the representation type for positions

      if ( ( simulation % dyntype .eq. 2 ) .or. &
           ( simulation % dyntype .eq. 3 ) .or. &
           ( simulation % dyntype .eq. 4 ) ) then

           representation % position = .false.
           representation % force = .true.
           representation % momentum = .true. 

      else

           representation % position = .true.
           representation % force = .true.
           representation % momentum = .true. 

      end if

      ! and store the current values

      call get( state, position = position, momentum = momentum, mass = mass,  &
                energy = energy, Cell = cell, reciprocalCell = reciprocalCell, &
                representation = representation )
            ! no need for forces just yet

      ! calculate the kinetic energy

      energy % kinetic = zero

      if ( simulation % dyntype == 4 ) then 

         do i=1, nAtoms - 2

            energy % kinetic = energy % kinetic +       &
              half * dot_product( momentum(1:3,i), momentum(1:3,i) ) / mass(i)

         end do

         ! in this case the shell atoms have no kinetic energy

         momentum(:,nAtoms-1) = zero
         momentum(:,nAtoms) = zero

         shell_stress = zero

         do i=nAtoms-1, nAtoms

            shell_stress(1,1) = shell_stress(1,1) + mass(i) * position(1,i) * position(1,i)
            shell_stress(1,2) = shell_stress(1,2) + mass(i) * position(1,i) * position(2,i)
            shell_stress(1,3) = shell_stress(1,3) + mass(i) * position(1,i) * position(3,i)

            shell_stress(2,1) = shell_stress(2,1) + mass(i) * position(2,i) * position(1,i)
            shell_stress(2,2) = shell_stress(2,2) + mass(i) * position(2,i) * position(2,i)
            shell_stress(2,3) = shell_stress(2,3) + mass(i) * position(2,i) * position(3,i)

            shell_stress(3,1) = shell_stress(3,1) + mass(i) * position(3,i) * position(1,i)
            shell_stress(3,2) = shell_stress(3,2) + mass(i) * position(3,i) * position(2,i)
            shell_stress(3,3) = shell_stress(3,3) + mass(i) * position(3,i) * position(3,i)

         end do

      else

         do i=1, nAtoms 

            energy % kinetic = energy % kinetic +       &
              half * dot_product( momentum(1:3,i), momentum(1:3,i) ) / mass(i)

         end do

      end if

      ! now, if the kinetic energy is small (a signal that momenta have not been 
      ! read from elsewhere) then generate atomic momenta here; the criterion is
      ! that if the temperature deduced from the current kinetic energy is less
      ! than 10% of the target temperature, and this is a start calculation, we
      ! generate momenta here

      ! if dyntype > 0 we might have one or more thermostats

      if ( simulation % dyntype > 0 ) then

         call get( state, nThermostats = nThermostats )

         if ( nThermostats  > 0 ) then

            allocate( thermostat( nThermostats ) )

            call get( state, thermostat = thermostat )

            ! calculate the initial thermostat energy

            call calculateThermostatEnergy

         end if

      end if

      if ( ( simulation % dyntype == 2 ) .or.                          &
           ( simulation % dyntype == 3 ) .or.                          &
           ( simulation % dyntype == 4 ) ) then   ! we have a barostat

         flat_stress = pack( simulation % stress % lattice, .true. )
         flat_metric = pack( cell % metricTensor, .true. )

         stress_energy = half * dot_product( flat_stress, flat_metric )

         energy % barostat = simulation % pressure * cell % volume + stress_energy

      end if

      call set( state, energy = energy )

      if ( debug ) write(*,*) 'Exiting MD_setUp()'

end subroutine MD_setUp

