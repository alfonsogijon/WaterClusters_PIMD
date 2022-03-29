!*****************************************************************************
subroutine calculateInternalPressure( state )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: State (PRIVATE routine)
!
!  Created on: Fri 18 Jan 11:29:50 2019  by Eduardo R. Hernandez
!
!  Last Modified: Wed 10 Apr 16:09:11 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! As its name indicates, this subroutine calculates the internal pressure on
! the system, resulting from the ionic kinetic contribution and the cell
! stress tensor. Obviously, it only makes sense to call this routine (from 
! the inquire routine) if the stress tensor has been computed by the client
! program and passed on to the state instance; otherwise only the kinetic 
! contribution will be correctly calculated, and the calculated pressure will
! be smaller than the real value.
!
! Modified to compute both the pressure derived from the potential part of the 
! stress as well as that derived from the total internal stress (including both
! the kinetic part and the potential part).
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
!  local variables

      integer :: i, n

      real (dp) :: factor

      real (dp), dimension (3) :: vector

      real (dp), dimension (3,3) :: kinetic_stress

!*****************************************************************************
!  Start of subroutine

      if ( debug ) write(*,*) 'Entering calculateInternalPressure()'

      ! we need the atomic momenta in Lattice reprsentation; make sure they are

      call transformMomenta2Lattice( state )

      ! calculate the kinetic stress in lattice rep

      kinetic_stress = zero

      do i=1, state % nAtoms 

         factor = one / state % mass(i)

         vector = matmul( state % reciprocalCell % metricTensor,       &
                          state % momentum(1:3,i) )

         kinetic_stress(1,1) = kinetic_stress(1,1) +                   &
            factor * vector(1) * vector(1) 
         kinetic_stress(1,2) = kinetic_stress(1,2) +                   &
            factor * vector(1) * vector(2) 
         kinetic_stress(1,3) = kinetic_stress(1,3) +                   &
            factor * vector(1) * vector(3) 

         kinetic_stress(2,1) = kinetic_stress(2,1) +                   &
            factor * vector(2) * vector(1) 
         kinetic_stress(2,2) = kinetic_stress(2,2) +                   &
            factor * vector(2) * vector(2) 
         kinetic_stress(2,3) = kinetic_stress(2,3) +                   &
            factor * vector(2) * vector(3) 

         kinetic_stress(3,1) = kinetic_stress(3,1) +                   &
            factor * vector(3) * vector(1) 
         kinetic_stress(3,2) = kinetic_stress(3,2) +                   &
            factor * vector(3) * vector(2) 
         kinetic_stress(3,3) = kinetic_stress(3,3) +                   &
            factor * vector(3) * vector(3) 

      end do

      factor = one

      if ( state % nThermostats > 0 ) then

         do n=1, state % nThermostats

            factor = factor * state % thermostat(n) % position

         end do

      end if

      kinetic_stress = kinetic_stress / ( factor * factor )

      state % totalInternalStress % calculated = .true.

      state % totalInternalStress % lattice = ( kinetic_stress -       &
            two * state % stress % lattice ) / state % Cell % volume

      state % totalInternalStress % cartesian =                        &
         matmul( state % totalInternalStress % lattice, state % Cell % metricTensor )

      state % totalInternalStress % pressure = ( one / three ) * (     &
              state % totalInternalStress % cartesian(1,1) +           &
              state % totalInternalStress % cartesian(2,2) +           &
              state % totalInternalStress % cartesian(3,3) )

      ! the following is the same, but without the kinetic contribution

      state % stress % pressure = ( two / three ) *                    &
         dot_product( pack( state % stress % lattice, .true. ),        &
                      pack( state % Cell % metricTensor, .true. ) )

      if ( debug ) write(*,*) 'Exiting calculateInternalPressure()'

end subroutine calculateInternalPressure

