!*****************************************************************************
subroutine NosePoincare_B( state )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Created on: Fri  4 May 12:19:12 2018  by Eduardo R. Hernandez
!
!  Last Modified: Fri 24 Aug 10:26:10 2018
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
!  shared variables

      type ( state_type ), intent (INOUT) :: state

!*****************************************************************************
!  internal variables

      integer :: i

      real (dp) :: const
      real (dp) :: factor

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering NosePoincare_B()'

! get from state all the required variables; of course we assume masses are unchanged

      call get( state, force = force, momentum = momentum, thermostat = thermostat, &
                energy = energy )

! update the thermostat momentum

      thermostat(1) % momentum = thermostat(1) % momentum +           &
        half_delta_t * ( energy % H0 + energy % kinetic -             &
                            energy % potential -                      &
        nDegreesOfFreedom * boltzmann_k * simulation % temperature * &
         ( one + log( thermostat(1) % position ) ) -                         &
           half * thermostat(1) % momentum * thermostat(1) % momentum /      &
                  thermostat(1) % mass )

! update the ion momenta to full step

      momentum = momentum + half_delta_t * thermostat(1) % position * force

! update the energy components

! first calculate kinetic energy with momenta at full-step

      energy % kinetic = 0.0_dp

      const = half / ( thermostat(1) % position * thermostat(1) % position )

      do i = 1, nAtoms

         factor = const / mass(i) 

         energy % kinetic = energy % kinetic + factor *   &
          dot_product( momentum(1:3,i), momentum(1:3,i) )

      end do

      energy % thermostat =                                                   &
          half * thermostat(1) % momentum * thermostat(1) % momentum /        &
            thermostat(1) % mass +    &
        boltzmann_k * simulation % temperature * nDegreesOfFreedom * log( thermostat(1) % position ) 

      energy % constantOfMotion =                                             &
          thermostat(1) % position * ( energy % kinetic +                     &
                                   energy % potential +                       &
                                   energy % thermostat - energy % H0 ) 

! finally, store all updated variables in state

      call set( state, momentum = momentum, thermostat = thermostat, energy = energy )

      if ( debug ) write(*,*) 'Exiting NosePoincare_B()'

end subroutine NosePoincare_B

