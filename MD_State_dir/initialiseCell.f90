!*****************************************************************************
subroutine initialiseCell( state, LatticeVectors )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: State (PRIVATE routine)
!
!  Created on: Mon 19 Mar 10:53:13 2018  by Eduardo R. Hernandez
!
!  Last Modified: Wed 10 Apr 16:00:42 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! From the given lattice vectors this subroutine constructs an instance 
! of the Cell structure (cell vectors, angles, moduli, etc) as well as 
! its reciprocal cell
! 
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  arguments

      type ( state_type ), intent (INOUT) :: state

      real (dp), dimension (3,3), intent (IN) :: LatticeVectors 
            ! a lattice vector per column, a(1:3) = LatticeVectors(1:3,1), etc

!*****************************************************************************
!  local variables

      real (dp) :: cos_angle

      real (dp) :: vector(3)

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering initialiseCell()'

      state % Cell % a(1:3) = LatticeVectors(1:3,1)
      state % Cell % b(1:3) = LatticeVectors(1:3,2)
      state % Cell % c(1:3) = LatticeVectors(1:3,3)

      state % Cell % metricTensor(1,1) = dot_product( state % Cell % a, state % Cell % a )
      state % Cell % metricTensor(1,2) = dot_product( state % Cell % a, state % Cell % b )
      state % Cell % metricTensor(1,3) = dot_product( state % Cell % a, state % Cell % c )

      state % Cell % metricTensor(2,1) = dot_product( state % Cell % b, state % Cell % a )
      state % Cell % metricTensor(2,2) = dot_product( state % Cell % b, state % Cell % b )
      state % Cell % metricTensor(2,3) = dot_product( state % Cell % b, state % Cell % c )

      state % Cell % metricTensor(3,1) = dot_product( state % Cell % c, state % Cell % a )
      state % Cell % metricTensor(3,2) = dot_product( state % Cell % c, state % Cell % b )
      state % Cell % metricTensor(3,3) = dot_product( state % Cell % c, state % Cell % c )

      state % Cell % a_modulus = sqrt( dot_product( state % Cell % a, state % Cell % a ) )
      state % Cell % b_modulus = sqrt( dot_product( state % Cell % b, state % Cell % b ) )
      state % Cell % c_modulus = sqrt( dot_product( state % Cell % c, state % Cell % c ) )

      cos_angle = dot_product( state % Cell % b, state % Cell % c ) /       &
                ( state % Cell % b_modulus * state % Cell % c_modulus )

      state % Cell % alpha = acos ( cos_angle ) * rad_to_deg

      cos_angle = dot_product( state % Cell % a, state % Cell % c ) /       &
                ( state % Cell % a_modulus * state % Cell % c_modulus )

      state % Cell % beta = acos ( cos_angle ) * rad_to_deg

      cos_angle = dot_product( state % Cell % a, state % Cell % b ) /       &
                ( state % Cell % a_modulus * state % Cell % b_modulus )

      state % Cell % gamma = acos ( cos_angle ) * rad_to_deg

      vector(1) = state % Cell % b(2) * state % Cell % c(3) -               &
                  state % Cell % c(2) * state % Cell % b(3)
      vector(2) = state % Cell % b(3) * state % Cell % c(1) -               &
                  state % Cell % c(3) * state % Cell % b(1)
      vector(3) = state % Cell % b(1) * state % Cell % c(2) -               &
                  state % Cell % c(1) * state % Cell % b(2)

      state % Cell % volume = dot_product( vector, state % Cell % a )

! now that the Cell has been constructed, obtain the reciprocal cell

      call setReciprocalCell( state ) 

      if ( debug ) write(*,*) 'Exiting initialiseCell()'

end subroutine initialiseCell

