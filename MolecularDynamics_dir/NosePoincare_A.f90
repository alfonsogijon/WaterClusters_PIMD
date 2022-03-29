!*****************************************************************************
subroutine NosePoincare_A( state )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Created on: Fri  4 May 12:19:12 2018  by Eduardo R. Hernandez
!
!  Last Modified: Fri 24 Aug 10:25:04 2018
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! Implements the update of momenta to half-step and positions to full step
! according to the classical Nose-Poincare method for NVT (canonical)
! simulations. The rest of the algorithm, for the update to full-step of the
! momenta, after force evaluation at the new positions, is implemented in the
! sequel subroutine NosePoincare_B.
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

      real (dp), parameter :: tolerance = 1.0e-7_dp

      real (dp) :: const
      real (dp) :: factor
      real (dp) :: thermostat_temp, thermostat_p_new

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering NosePoincare_A()'

! get from state all the required variables; of course we assume masses are unchanged

      call get( state, position = position, momentum = momentum, force = force, &
                energy = energy, thermostat = thermostat )

! first advance momenta by half a step

      momentum = momentum +                            &
            half_delta_t * thermostat(1) % position * force

! now calculate kinetic energy with new momenta at half-step

      energy % kinetic = 0.0_dp

      const = half / ( thermostat(1) % position * thermostat(1) % position )

      do i = 1, nAtoms

         factor = const / mass(i) 

         energy % kinetic = energy % kinetic + factor *   &
          dot_product( momentum(1:3,i), momentum(1:3,i) )

      end do

! now advance the thermostat momentum by half a step

      factor = half_delta_t * ( nDegreesOfFreedom * boltzmann_k                              &
                        * simulation % temperature *                                         &
                 ( one + log( thermostat(1) % position ) ) -                                 &
                 energy % kinetic + energy % potential - energy % H0 ) -                     &
                 thermostat(1) % momentum 

      thermostat(1) % momentum =                                                                   &
        -two * factor / ( one + sqrt( one - factor * delta_t / thermostat(1) % mass ) ) 

! now advance the thermostat position to full-step

      thermostat_temp = thermostat(1) % position
 
      do 

         thermostat_p_new = thermostat(1) % position + half_delta_t * ( thermostat_temp +          &
                thermostat(1) % position ) * thermostat(1) % momentum / thermostat(1) % mass

         if ( abs( thermostat_p_new - thermostat_temp ) < tolerance ) exit

         thermostat_temp = thermostat_p_new

      end do

! now we can advance atomic positions

      factor = half_delta_t * ( one / thermostat(1) % position + one / thermostat_p_new )

! update the positions of all atoms with current momenta

      do i=1, nAtoms

         position(1:3,i) = position(1:3,i) + factor * momentum(1:3,i) / mass(i)
 
      end do

! update the value of the kinetic energy (ionic momenta not changed, but thermostat position has)

      factor = thermostat(1) % position * thermostat(1) % position / ( thermostat_p_new * thermostat_p_new )

      energy % kinetic = factor * energy % kinetic

      thermostat(1) % position = thermostat_p_new

! finally, store the new coordinates and half-step momenta in state

      call set( state, position = position, momentum = momentum, thermostat = thermostat, energy = energy )

      if ( debug ) write(*,*) 'Exiting NosePoincare_A()'

end subroutine NosePoincare_A

