!*****************************************************************************
subroutine updateCell( Cell, MetricTensor, updateReciprocal )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: MolecularDynamics (PRIVATE routine)
!
!  Created on: Mon 19 Mar 10:53:13 2018  by Eduardo R. Hernandez
!
!  Last Modified: Fri 21 Jun 12:09:28 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! From the current metric tensor, this subroutine updates (recalculates) the
! current cell (cell vectors, angles, moduli, etc) as well as its reciprocal
! cell; this is only needed in case of flexible-cell dynamics.
! 
! The new cell is obtained from the current metric tensor and the previous
! cell by calculating a symmetric strain tensor (hence containing only 
! deformation and no rotatations) that gives new lattice vectors compatible
! with the current metric tensor. The current metric tensor is:
! 
!     G = H^t H 
!
! where G is the metric tensor, and H is the sought cell vector matrix (each vector
! a column of H). The idea is to obtain H from H_0, the previous cell, and 
! a symmetric (hence no rotation, only distortion) strain tensor 
!
!     H = (I + epsilon) H_0
!
! where I is the 3x3 unit matrix and the strain tensor is expected to have
! components abs(epsilon_ij) << 1, since from one step to the next we expect
! the lattice vectors to have changed little. Then, we have that
!     
!     G = H^t H = G_0 + 2 H_0^t epsilon H_0 + O(epsilon^2)
!
! where we disregard terms higher than linear in epsilon (assumed small). Then
! we can solve for epsilon as
!
!     epsilon ~ 1/2 (H_0^t)^-1[G - G_0] H_0^-1
! 
! Since H_0 is expected to be close to H, and consequently epsilon small, 
! the above expression should be a very good guess, but just to be sure that 
! we approach the correct H to within a specified tolerance, the above 
! expression can be used iteratively until the tolerance is met. 
! from which the new cell is obtained as H = e H_o. 
!
! All above matrices are simple 3x3 matrices. 
! 
! As of 2 Apr 2019, this routine optionally calculates the reciprocal cell,
! which is sometimes needed. 
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  arguments

      type ( lattice_type ), intent (INOUT) :: Cell

      real (dp), dimension (3,3), intent (IN) :: MetricTensor 

      logical, optional, intent (IN) :: updateReciprocal

!*****************************************************************************
!  local variables

      logical :: converged

      integer :: nIterations

      integer, parameter :: nMaxIterations = 20

      real (dp) :: cos_angle, deltaG2, volume

      real (dp), parameter :: tolerance = 1.0e-18_dp

      real (dp), dimension (3) :: a, b, c, vector
      real (dp), dimension (3) :: a0, b0, c0

      real (dp), dimension (3,3) :: epsilon, epsilon2, H0inv, MetricTensor0
      real (dp), dimension (3,3) :: deltaG, identity, matrix, HGH

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering updateCell()'

! enter an iterative loop until we home-in on lattice vectors compatible
! with the current MetricTensor      

      nIterations = 1
      converged = .false.

      identity = zero
      identity(1,1) = one
      identity(2,2) = one
      identity(3,3) = one

      epsilon2 = zero

! we initialise the process with the current cell vectors

      a0 = Cell % a
      b0 = Cell % b
      c0 = Cell % c

      a = a0
      b = b0
      c = c0

! first we must calculate the reciprocal cell vectors

      vector(1) = b(2) * c(3) - c(2) * b(3)
      vector(2) = b(3) * c(1) - c(3) * b(1)
      vector(3) = b(1) * c(2) - c(1) * b(2)

      volume = dot_product( a, vector )

      H0inv(:,1) = vector / volume

      vector(1) = c(2) * a(3) - a(2) * c(3)
      vector(2) = c(3) * a(1) - a(3) * c(1)
      vector(3) = c(1) * a(2) - a(1) * c(2)

      H0inv(:,2) = vector / volume

      vector(1) = a(2) * b(3) - b(2) * a(3)
      vector(2) = a(3) * b(1) - b(3) * a(1)
      vector(3) = a(1) * b(2) - b(1) * a(2)

      H0inv(:,3) = vector / volume

      matrix = matmul( H0inv, MetricTensor )
      
      HGH = matmul( matrix, transpose(H0inv) )

      do while ( .not. converged )

! check for convergence 

         MetricTensor0(1,1) = dot_product( a, a )
         MetricTensor0(1,2) = dot_product( a, b )
         MetricTensor0(1,3) = dot_product( a, c )

         MetricTensor0(2,1) = dot_product( b, a )
         MetricTensor0(2,2) = dot_product( b, b )
         MetricTensor0(2,3) = dot_product( b, c )

         MetricTensor0(3,1) = dot_product( c, a )
         MetricTensor0(3,2) = dot_product( c, b )
         MetricTensor0(3,3) = dot_product( c, c )

         deltaG = MetricTensor - MetricTensor0

         deltaG2 = deltaG(1,1) * deltaG(1,1) +         &
                   deltaG(2,2) * deltaG(2,2) +         &
                   deltaG(3,3) * deltaG(3,3) +         &
                   deltaG(1,2) * deltaG(1,2) +         &
                   deltaG(1,3) * deltaG(1,3) +         &
                   deltaG(2,3) * deltaG(2,3)

         if ( deltaG2 .le. tolerance ) converged = .true.

! now we can obtain the linear approximation to the strain

         epsilon = half * ( HGH - identity - epsilon2 )

! epsilon must be a symmetric tensor, so lets impose it

         matrix = epsilon

         epsilon = half * ( transpose( matrix ) + matrix )

         epsilon2 = matmul( epsilon, epsilon )

! now obtain the new lattice vectors 

         a = a0 + matmul( epsilon, a0 )         
         b = b0 + matmul( epsilon, b0 )         
         c = c0 + matmul( epsilon, c0 )         

         nIterations = nIterations + 1

         if ( nIterations .gt. nMaxIterations ) then 

            write(*,*) 'Maximum number of iterations exceeded '
            write(*,*) 'exceeded in updateCell: stopping '

            stop

         end if

      end do

      ! now we can construct the new cell parameters

      Cell % a = a
      Cell % b = b
      Cell % c = c

      ! and then proceed to calculate the rest of the contents of Cell

      Cell % metricTensor(1,1) = dot_product( Cell % a, Cell % a )
      Cell % metricTensor(1,2) = dot_product( Cell % a, Cell % b )
      Cell % metricTensor(1,3) = dot_product( Cell % a, Cell % c )

      Cell % metricTensor(2,1) = dot_product( Cell % b, Cell % a )
      Cell % metricTensor(2,2) = dot_product( Cell % b, Cell % b )
      Cell % metricTensor(2,3) = dot_product( Cell % b, Cell % c )

      Cell % metricTensor(3,1) = dot_product( Cell % c, Cell % a )
      Cell % metricTensor(3,2) = dot_product( Cell % c, Cell % b )
      Cell % metricTensor(3,3) = dot_product( Cell % c, Cell % c )

      Cell % a_modulus = sqrt( dot_product( Cell % a, Cell % a ) )
      Cell % b_modulus = sqrt( dot_product( Cell % b, Cell % b ) )
      Cell % c_modulus = sqrt( dot_product( Cell % c, Cell % c ) )

      cos_angle = dot_product( Cell % b, Cell % c ) /                       &
                ( Cell % b_modulus * Cell % c_modulus )

      Cell % alpha = acos ( cos_angle ) * rad_to_deg

      cos_angle = dot_product( Cell % a, Cell % c ) /                       &
                ( Cell % a_modulus * Cell % c_modulus )

      Cell % beta = acos ( cos_angle ) * rad_to_deg

      cos_angle = dot_product( Cell % a, Cell % b ) /                       &
                ( Cell % a_modulus * Cell % b_modulus )

      Cell % gamma = acos ( cos_angle ) * rad_to_deg

      vector(1) = Cell % b(2) * Cell % c(3) - Cell % c(2) * Cell % b(3)
      vector(2) = Cell % b(3) * Cell % c(1) - Cell % c(3) * Cell % b(1)
      vector(3) = Cell % b(1) * Cell % c(2) - Cell % c(1) * Cell % b(2)

      Cell % volume = dot_product( vector, Cell % a )

! the reciprocal cell is now calculated in setReciprocalCell in the MD_State module, 
! but, if required, calculate it here

      if ( present( updateReciprocal  ) ) then 

         if ( updateReciprocal ) then

            vector(1) = Cell % b(2) * Cell % c(3) - Cell % c(2) * Cell % b(3)
            vector(2) = Cell % b(3) * Cell % c(1) - Cell % c(3) * Cell % b(1)
            vector(3) = Cell % b(1) * Cell % c(2) - Cell % c(1) * Cell % b(2)

            a = vector / Cell % volume
            reciprocalCell % a = a

            vector(1) = Cell % c(2) * Cell % a(3) - Cell % a(2) * Cell % c(3)
            vector(2) = Cell % c(3) * Cell % a(1) - Cell % a(3) * Cell % c(1)
            vector(3) = Cell % c(1) * Cell % a(2) - Cell % a(1) * Cell % c(2)

            b = vector / Cell % volume
            reciprocalCell % b = b

            vector(1) = Cell % a(2) * Cell % b(3) - Cell % b(2) * Cell % a(3)
            vector(2) = Cell % a(3) * Cell % b(1) - Cell % b(3) * Cell % a(1)
            vector(3) = Cell % a(1) * Cell % b(2) - Cell % b(1) * Cell % a(2)

            c = vector / Cell % volume
            reciprocalCell % c = c

            reciprocalCell % MetricTensor(1,1) = dot_product( a, a )
            reciprocalCell % MetricTensor(1,2) = dot_product( a, b )
            reciprocalCell % MetricTensor(1,3) = dot_product( a, c )

            reciprocalCell % MetricTensor(2,1) = dot_product( b, a )
            reciprocalCell % MetricTensor(2,2) = dot_product( b, b )
            reciprocalCell % MetricTensor(2,3) = dot_product( b, c )

            reciprocalCell % MetricTensor(3,1) = dot_product( c, a )
            reciprocalCell % MetricTensor(3,2) = dot_product( c, b )
            reciprocalCell % MetricTensor(3,3) = dot_product( c, c )

            reciprocalCell % a_modulus = sqrt( dot_product( a, a ) )
            reciprocalCell % b_modulus = sqrt( dot_product( b, b ) )
            reciprocalCell % c_modulus = sqrt( dot_product( c, c ) )

            cos_angle = dot_product( reciprocalCell % b, reciprocalCell % c ) /     &
                   ( reciprocalCell % b_modulus * reciprocalCell % c_modulus )

            reciprocalCell % alpha = acos ( cos_angle ) * rad_to_deg

            cos_angle = dot_product( reciprocalCell % a, reciprocalCell % c ) /  &
                   ( reciprocalCell % a_modulus * reciprocalCell % c_modulus )

            reciprocalCell % beta = acos ( cos_angle ) * rad_to_deg

            cos_angle = dot_product( reciprocalCell % a, reciprocalCell % b ) /  &
                   ( reciprocalCell % a_modulus * reciprocalCell % b_modulus )

            reciprocalCell % gamma = acos ( cos_angle ) * rad_to_deg

            vector(1) = b(2) * c(3) - c(2) * b(3)
            vector(2) = b(3) * c(1) - c(3) * b(1)
            vector(3) = b(1) * c(2) - c(1) * b(2)

            reciprocalCell % volume = dot_product( vector, reciprocalCell % a )

         end if
      end if

      if ( debug ) write(*,*) 'Exiting updateCell()'

end subroutine updateCell
