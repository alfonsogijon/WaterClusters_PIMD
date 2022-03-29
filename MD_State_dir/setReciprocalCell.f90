!*****************************************************************************
subroutine setReciprocalCell( state )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: State (PRIVATE routine)
!
!  Created on: Mon 19 Mar 10:53:13 2018  by Eduardo R. Hernandez
!
!  Last Modified: Wed 10 Apr 15:59:29 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! Once a Cell has been defined for the current instance of State (input), this
! subroutine constructs the corresponding reciprocal cell, contained in 
! state % reciprocalCell. The reciprocal cell vectors are defined in the usual
! way except for a factor of 2pi, which is not used in our case (this makes the
! matrix of reciprocal lattice vectors the inverse of the matrix of direct 
! lattice vectors). 
! 
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  arguments

      type ( state_type ), intent (INOUT) :: state

!*****************************************************************************
!  local variables

      real (dp) :: cos_angle

      real (dp) :: vector(3)

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering setReciprocalCell()'

! calculate the reciprocal lattice vectors
! these are the reciprocal lattice vectors but for a factor of 2pi

      vector(1) = state % Cell % b(2) * state % Cell % c(3) -               &
                  state % Cell % c(2) * state % Cell % b(3)
      vector(2) = state % Cell % b(3) * state % Cell % c(1) -               &
                  state % Cell % c(3) * state % Cell % b(1)
      vector(3) = state % Cell % b(1) * state % Cell % c(2) -               &
                  state % Cell % c(1) * state % Cell % b(2)

      state % reciprocalCell % a = vector / state % Cell % volume

      vector(1) = state % Cell % c(2) * state % Cell % a(3) - state % Cell % a(2) * state % Cell % c(3)
      vector(2) = state % Cell % c(3) * state % Cell % a(1) - state % Cell % a(3) * state % Cell % c(1)
      vector(3) = state % Cell % c(1) * state % Cell % a(2) - state % Cell % a(1) * state % Cell % c(2)

      state % reciprocalCell % b = vector / state % Cell % volume

      vector(1) = state % Cell % a(2) * state % Cell % b(3) - state % Cell % b(2) * state % Cell % a(3)
      vector(2) = state % Cell % a(3) * state % Cell % b(1) - state % Cell % b(3) * state % Cell % a(1)
      vector(3) = state % Cell % a(1) * state % Cell % b(2) - state % Cell % b(1) * state % Cell % a(2)

      state % reciprocalCell % c = vector / state % Cell % volume

      state % reciprocalCell % a_modulus =                                    &
          sqrt( dot_product( state % reciprocalCell % a, state % reciprocalCell % a ) )
      state % reciprocalCell % b_modulus =                                    &
          sqrt( dot_product( state % reciprocalCell % b, state % reciprocalCell % b ) )
      state % reciprocalCell % c_modulus =                                    &
          sqrt( dot_product( state % reciprocalCell % c, state % reciprocalCell % c ) )

      cos_angle = dot_product( state % reciprocalCell % b, state % reciprocalCell % c ) /   &
                ( state % reciprocalCell % b_modulus * state % reciprocalCell % c_modulus )

      state % reciprocalCell % alpha = acos ( cos_angle ) * rad_to_deg

      cos_angle = dot_product( state % reciprocalCell % a, state % reciprocalCell % c ) /   &
                ( state % reciprocalCell % a_modulus * state % reciprocalCell % c_modulus )

      state % reciprocalCell % beta = acos ( cos_angle ) * rad_to_deg

      cos_angle = dot_product( state % reciprocalCell % a, state % reciprocalCell % b ) /   &
                ( state % reciprocalCell % a_modulus * state % reciprocalCell % b_modulus )
   
      state % reciprocalCell % gamma = acos ( cos_angle ) * rad_to_deg

      state % reciprocalCell % volume = one / state % Cell % volume

! calculate the reciprocal metric tensor

      state % reciprocalCell % metricTensor(1,1) =                                 &
           dot_product( state % reciprocalCell % a, state % reciprocalCell % a )
      state % reciprocalCell % metricTensor(1,2) =                                 &
           dot_product( state % reciprocalCell % a, state % reciprocalCell % b )
      state % reciprocalCell % metricTensor(1,3) =                                 &
           dot_product( state % reciprocalCell % a, state % reciprocalCell % c )

      state % reciprocalCell % metricTensor(2,1) =                                 &
           dot_product( state % reciprocalCell % b, state % reciprocalCell % a )
      state % reciprocalCell % metricTensor(2,2) =                                 &
           dot_product( state % reciprocalCell % b, state % reciprocalCell % b )
      state % reciprocalCell % metricTensor(2,3) =                                 &
           dot_product( state % reciprocalCell % b, state % reciprocalCell % c )

      state % reciprocalCell % metricTensor(3,1) =                                 &
           dot_product( state % reciprocalCell % c, state % reciprocalCell % a )
      state % reciprocalCell % metricTensor(3,2) =                                 &
           dot_product( state % reciprocalCell % c, state % reciprocalCell % b )
      state % reciprocalCell % metricTensor(3,3) =                                 &
           dot_product( state % reciprocalCell % c, state % reciprocalCell % c )

      if ( debug ) write(*,*) 'Exiting setReciprocalCell()'

end subroutine setReciprocalCell

