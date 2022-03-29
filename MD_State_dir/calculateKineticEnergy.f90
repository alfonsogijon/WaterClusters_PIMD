!*****************************************************************************
subroutine calculateKineticEnergy( state )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: State (PRIVATE routine) 
!
!  Created on: Wed 27 Jun 17:37:12 2018  by Eduardo R. Hernandez
!
!  Last Modified: Thu 17 Jan 19:03:50 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! Given the momentum contained in state (and if appropriate, any thermostats)
! this routine calculates the kinetic energy of the system.
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
!  variables

      integer :: i, n

      real (dp) E_kin, factor, prefactor

!*****************************************************************************
!  Start of subroutine

      if ( debug ) write(*,*) 'Entering calculateKineticEnergy()'

      ! consider the different possible cases regarding thermostats

      E_kin = zero

      ! make sure momenta are in Cartesian coordinates

      call transformMomenta2Cartesian( state )

      do i=1, state % nAtoms

         factor = half / state % mass(i)

         E_kin = E_kin + factor *                                         &
         dot_product( state % momentum(1:3,i), state % momentum(1:3,i) )

      end do

      prefactor = one

      if ( state % nThermostats > 0 ) then    ! this is the standard NVE case

         do n=1, state % nThermostats

            prefactor = prefactor / state % thermostat(n) % position

         end do

         prefactor = prefactor * prefactor

      end if

      state % energy % kinetic = prefactor * E_kin

      if ( debug ) write(*,*) 'Exiting calculateKineticEnergy()'

end subroutine calculateKineticEnergy

