!*****************************************************************************
subroutine MSD( nAtoms, mass, centroid0, centroid, msd_water )
  
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: State (PUBLIC routine)
!
!  Created on: May 2019 by Alfonso Gijón
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This subroutine calculates the MSD of the water molecules
!
!*****************************************************************************
!  modules used

      use omp_lib  

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  variables

      integer, intent(IN) :: nAtoms

      real (dp), dimension (:), intent(IN) :: mass

      real (dp), dimension (:,:), intent(IN) :: centroid0

      real (dp), dimension (:,:), intent(IN) :: centroid      

      real(dp), intent (OUT) :: msd_water

!*****************************************************************************
!  local variables

      integer imol
      integer iat
      integer j

      integer nmols

      real(dp), parameter :: mh2o = 18.02d0      !!!!!!!!!! Definir correctamente

      real (dp), dimension (3) :: rcm, rcm0, rdif

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering MSD()'

      nmols = nAtoms/3

! Compute MSD of water molecules

      msd_water = zero
      
      do imol = 1, nmols

         ! Compute center of mass
         
         rcm = zero
         rcm0 = zero
         do j = 1, 3
            iat = j + 3*(imol-1)
            rcm(:) = rcm(:) + mass(iat)*centroid(:,iat)
            rcm0(:) = rcm0(:) + mass(iat)*centroid0(:,iat)
         end do
         rcm(:) = rcm(:) / mh2o
         rcm0(:) = rcm0(:) / mh2o         

         rdif = rcm - rcm0
         msd_water = msd_water + dot_product(rdif,rdif)
            
      end do
      
      msd_water = msd_water / nmols

      if ( debug ) write(*,*) 'Exiting MSD()'
      
end subroutine MSD
