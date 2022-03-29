!*****************************************************************************
subroutine createState( state, nAtoms,                                       &
     nThermostats, thermostatMass, barostatMass) 
!*****************************************************************************
!
!  Project: MolecularDynamics 
!
!  Module: State (PUBLIC routine)
!
!  Created on: Tue 12 Jun 11:42:03 2018  by Eduardo R. Hernandez
!
!  Last Modified: Thu 10 Oct 09:32:47 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This subroutine is called to create a new instance of State; the state instance
! is allocated, but not "filled"; this must be done by either initialiseState (in
! case of a brand new simulation) or by readRestart (in case of a continuation
! simulation).
!
! variables:
!
! type ( state_type ), intent (INOUT) :: state
!     this is the state that is going to be allocated and filled with the data
!     passed next.
!
! integer, intent (IN) :: nAtoms
!     the number of atoms in the system; used to dimension several arrays below
!
! integer, intent (IN), optional :: nThermostats
!     the number of thermostats in the thermostat chain. For NVE or NPH simulations
!     that have no thermostat, this variable is either set to zero or absent (default = 0);
!     for simulations with a thermostat or a thermostat chain (NVT, NPT), nThermostats
!     is a positive integer indicating the number of Nosé-Poincaré thermostats in the 
!     chain. 
!
! real (dp), dimension (:), intent (IN), optional :: thermostatMass
!     if nThermostats > 0, thermostatMass is an array of length nThermostats providing the
!     masses of each thermostat in the chain, in Daltons (1/12 of the mass of C12)
!
! real (dp), intent (IN), optional :: barostatMass
!     if the simulation is a constant-pressure one, this specifies the mass of the barostat
!     in Daltons (1/12 of the mass of C12)
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  arguments

      type ( state_type ), intent (INOUT) :: state

      integer, intent (IN) :: nAtoms

      integer, intent (IN), optional :: nThermostats
 
      real (dp), dimension (:), intent (IN), optional :: thermostatMass
      real (dp), intent (IN), optional :: barostatMass

!*****************************************************************************
!  local variables

      integer :: i

!*****************************************************************************
!  start subroutine

      if ( debug ) write(*,*) 'Entering createState()'  

      ! assign nAtoms and allocate variables that depend on it
       
      state % nAtoms = nAtoms
      if (nAtoms .gt. 1) then
         state % nDegreesOfFreedom = 3 * state % nAtoms - 3
      else
         state % nDegreesOfFreedom = 3 * state % nAtoms
      end if

      allocate( state % position( 3, state % nAtoms ) )
      allocate( state % force( 3, state % nAtoms ) )
      allocate( state % momentum( 3, state % nAtoms ) )
      allocate( state % mass( state % nAtoms ) )

      if ( present( nThermostats ) ) then

         state % nThermostats = nThermostats

         allocate( state % thermostat( nThermostats ) )

         do i = 1, state % nThermostats
            state % thermostat(i) % mass = massConversionFactor * thermostatMass(i)
            state % thermostat(i) % position = one
            state % thermostat(i) % momentum = zero
            state % thermostat(i) % a = one
            state % thermostat(i) % C = one
            write(*,*) 'thermostat mass ', i, state % thermostat(i) % mass
         end do 

      else

         state % nThermostats = 0  ! no thermostat is defined

      end if

      if ( present( barostatMass ) ) then

         state % barostat % mass = massConversionFactor * barostatMass
         write(*,*) 'barostat mass ', state % barostat % mass

         state % barostat % momentum = zero

      end if

      ! these are initialised to zero; they will be set prorperly if needed later

      state % energy % thermostat = zero
      state % energy % barostat = zero

      state % H0 % set = .false.
      state % H0 % E = zero

      if ( debug ) write(*,*) 'Exiting createState'

end subroutine createState

