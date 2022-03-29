!*****************************************************************************
subroutine MD_cleanUp
!*****************************************************************************
!
!  Project: MolecularDynamics 
!
!  Module: MolecularDynamics (PUBLIC routine)
!
!  Created on: Fri 16 Mar 14:35:15 2018  by Eduardo R. Hernandez
!
!  Last Modified: Wed 10 Apr 11:55:27 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This subroutine must be called at the end of an MD simulation to de-allocate
! temporary arrays that were allocated at the start of the simulation.
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
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  arguments

!*****************************************************************************
!  local variables

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering MD_cleanUp()'

      deallocate( position )
      deallocate( force )
      deallocate( momentum )
      deallocate( mass )

      if ( debug ) write(*,*) 'Exiting MD_cleanUp()'

end subroutine MD_cleanUp

