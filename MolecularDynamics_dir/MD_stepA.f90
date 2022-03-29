!*****************************************************************************
subroutine MD_stepA( state )
!*****************************************************************************
!
!  Project: MolecularDynamics    
! 
!  Module: MolecularDynamics (PUBLIC routine)
!
!  Created on: Fri  9 Mar 12:11:01 2018  by Eduardo R. Hernandez
!
!  Last Modified: Wed  9 Oct 11:16:15 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This subroutine performs an MD step in the numerical
! soluction of the equations of motion, namely, using the current positions, 
! momenta and forces, advances the momenta by half a time step; with the
! new velocities positions can be projected forward by a full 
! step; subsequently velocities must be updated by the remaining half-time step.
! In this scheme the whole MD step is carried out by a single routine call, but
! momenta and positions are out of phase by half a time step.
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  arguments

      type ( state_type ), intent (INOUT) :: state ! the current state of the system

!*****************************************************************************
!  local variables

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering MD_stepA()'

!  now perform appropriate dynamics step (first part)

      select case( simulation % dyntype )

      case( 0 )  ! this is the standard velocity-verlet NVE (microcanonical )dynamics

          call velocityVerlet_A( state )

      case( 1 )  ! this case is NVT dynamics implemented via the Nosé-Poincaré thermostat as
                 ! described in Bond, Leimkuhler and Laird (J. Comput. Phys. vol 151, 114 (1999))

          call NosePoincareChain_A( state )

      case( 2 )  ! this case is NPH dynamics using the Souza-Martins barostat 
                 ! Souza and Martins (Phys. Rev. B 55, 8733 (1997)),

          call NPSM_A( state )
          
      case( 3 )  ! this case is NPT dynamics, combining a Nosé-Poincaré thermostat and
                 ! a Souza-Martins barostat, see Hernández, J. Chem. Phys. 115, 10282 (2001)

          call NPSM_A( state )

      end select 

end subroutine MD_stepA

