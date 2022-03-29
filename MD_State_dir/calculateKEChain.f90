!*****************************************************************************
subroutine calculateKEChain( state_chain )
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
! Given the momentum contained in each state of the chain, this routine
! calculates the kinetic energy of the system.
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  variables

      type ( StateChain_type ), intent (INOUT) :: state_chain

!*****************************************************************************
!  variables

      integer :: i, n, j

      real (dp) E_kin, factor, prefactor

      real (dp) :: kinetic_energy

!*****************************************************************************
!  Start of subroutine

      if ( debug ) write(*,*) 'Entering calculateKEChain()'

      ! consider the different possible cases regarding thermostats

      E_kin = zero
      
      do j = 1, state_chain % nReplicas

         ! make sure momenta are in Cartesian coordinates

         !call transformMomenta2Cartesian( state_chain % replica(j) )

         do i=1, state_chain % nAtoms

            factor = half / state_chain % mass(i)

            E_kin = E_kin + factor *                                         &
                 dot_product( state_chain % replica(j) % momentum(1:3,i), &
                 state_chain % replica(j) % momentum(1:3,i) )

         end do

      end do

      ! compute prefactor due to the thermostat chain

      prefactor = one

      if ( state_chain % nThermostats > 0 ) then    ! this is the standard NVE case

         do n=1, state_chain % nThermostats

            prefactor = prefactor / state_chain % thermostat(n) % position

         end do

         prefactor = prefactor * prefactor

      end if

      kinetic_energy = prefactor * E_kin

      state_chain % energy % kinetic = kinetic_energy

      if ( debug ) write(*,*) 'Exiting calculateKEChain()'

end subroutine calculateKEChain
