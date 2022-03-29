!*****************************************************************************
subroutine calculateThermostatEnergy
!*****************************************************************************
!
!  Project: MolecularDynamics 
!
!  Module: MolecularDynamics (PRIVATE routine)
!
!  Created on: Thu 26 Jul 11:03:55 2018  by Eduardo R. Hernandez
!
!  Last Modified: Thu 17 Jan 19:10:06 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! A simple routine for calculating the energy of the Nosé-Poincaré thermostat
! or of a chain of such thermostats. This needs to be done at several places
! so it is convenient to bundle the calculation into a separate private 
! routine; this is it.
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  variables

      integer :: n

      real (dp) :: kbT
      real (dp) :: prod

      real (dp), dimension (:), allocatable :: denominator
      real (dp), dimension (:), allocatable :: potential
      real (dp), dimension (:), allocatable :: workvec

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering calculateThermostatEnergy()'

      allocate( denominator( nThermostats ) )
      allocate( potential( nThermostats ) )
      allocate( workvec( nThermostats ) )

      kbT = boltzmann_k * simulation % temperature

      if ( nThermostats .gt. 1 ) then

         do n=1, nThermostats-1

            prod = product( thermostat(n+1:nThermostats) % position )

            denominator(n) = half / ( thermostat(n) % mass * prod * prod )

            if ( n .eq. 1 ) then
               potential(n) = float( nDegreesOfFreedom ) * kbT * log( thermostat(n) % position )
            else
               potential(n) = float( nDegreesOfFreedom + n - 1 ) * kbT * log( thermostat(n) % position ) + &
                   ( thermostat(n) % a - thermostat(n) % position ) *                    &
                   ( thermostat(n) % a - thermostat(n) % position ) /                    &
                            ( two * thermostat(n) % C )
            end if

         end do

         denominator(nThermostats) = half / thermostat(nThermostats) % mass
         potential(nThermostats) =                                                       &
          float( nDegreesOfFreedom + nThermostats - 1 ) * kbT * log( thermostat(nThermostats) % position ) + &
                   ( thermostat(nThermostats) % a - thermostat(nThermostats) % position ) *                  &
                   ( thermostat(nThermostats) % a - thermostat(nThermostats) % position ) /                  &
                            ( two * thermostat(nThermostats) % C )

      else

         denominator(1) = half / thermostat(1) % mass
         potential(1) = float( nDegreesOfFreedom ) * kbT * log( thermostat(1) % position )

      end if

      workvec = thermostat % momentum * thermostat % momentum * denominator + potential

      energy % thermostat = sum( workvec )

      deallocate( denominator )
      deallocate( potential )
      deallocate( workvec )

      if ( debug ) write(*,*) 'Exiting calculateThermostatEnergy()'

end subroutine calculateThermostatEnergy

