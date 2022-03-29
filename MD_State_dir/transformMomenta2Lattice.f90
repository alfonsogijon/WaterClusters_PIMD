!*****************************************************************************
subroutine transformMomenta2Lattice( state )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: State (PRIVATE routine)
!
!  Created on: Thu  8 Mar 12:23:32 2018  by Eduardo R. Hernandez
!
!  Last Modified: Fri 22 Feb 16:10:04 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This subroutine transforms the atomic momenta from Cartesian to Lattice
! representation. Internally, State stores all atomic variables in Cartesian
! representation, but they may be requested by the calling program in 
! Lattice representation; this routine performs the appropriate transformation
! for momenta.
! 
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
   implicit none

!*****************************************************************************
!  interface variables

   type ( state_type ), intent (INOUT) :: state

!*****************************************************************************
!  local variables

   integer :: i

   real (dp), dimension (3) :: vector

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering transformMomenta2Lattice()'

      if ( state % representation % momentum ) then    ! momenta are in cartesian, so transform

         ! what we are referring to as lattice representation is 
         ! m_i G_ab \dot{s_i}^b 
         ! where the dot is a time derivative and summation over repeated indices is implied

         do i=1, state % nAtoms
         
            vector(1) = dot_product( state % reciprocalCell % a,         &
                                       state % momentum(1:3,i) )
            vector(2) = dot_product( state % reciprocalCell % b,         &
                                       state % momentum(1:3,i) )
            vector(3) = dot_product( state % reciprocalCell % c,         &
                                       state % momentum(1:3,i) )

            state % momentum(1:3,i) = matmul( state % cell % metricTensor, vector )

         end do

         state % representation % momentum = .false.   ! now they are in lattice

      else       ! momenta already in lattice rep, so nothing to do

         return

      end if

      if ( debug ) write(*,*) 'Exiting transformMomenta2Lattice()'

end subroutine transformMomenta2Lattice

