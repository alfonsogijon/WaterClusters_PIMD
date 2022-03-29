!*****************************************************************************
subroutine transformMomenta2Cartesian( state )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: State (PRIVATE routine)
!
!  Created on: Thu  8 Mar 12:23:32 2018  by Eduardo R. Hernandez
!
!  Last Modified: Fri 22 Feb 16:48:07 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
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

   real(dp), dimension (3) :: p, r, s, t

!*****************************************************************************
!  begin subroutine

      if ( debug ) write(*,*) 'Entering transformMomenta2Cartesian()'

      if ( state % representation % momentum ) then    ! momenta already cartesian; do nothing

         return

      else           ! no, we need to transform to cartesian, so do it here

         ! this transforms m_i G_ab \dot{s_ib} first to m_i \dot{s_ia}
         ! and then to m_i H_ab \dot_{s_ib} = m_i v_ia

         r = (/ state % cell % a(1), state % cell % b(1), state % cell % c(1) /)
         s = (/ state % cell % a(2), state % cell % b(2), state % cell % c(2) /)
         t = (/ state % cell % a(3), state % cell % b(3), state % cell % c(3) /)

         do i=1, state % nAtoms

            p = matmul( state % reciprocalCell % metricTensor, state % momentum(1:3,i) )

            state % momentum(1,i) = dot_product( r, p )
            state % momentum(2,i) = dot_product( s, p )
            state % momentum(3,i) = dot_product( t, p )

         end do

         state % representation % momentum = .true. 

      end if

      if ( debug ) write(*,*) 'Exiting transformMomenta2Cartesian()'

end subroutine transformMomenta2Cartesian

