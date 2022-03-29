!*****************************************************************************
subroutine PolymerProperties( nAtoms, nRep, polymer, &
     ! rmax_histo, density_o, density_h, rhor_o, rhor_h, goo, ghh, goh )
     rmax_histo, goo, ghh, goh )

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
! This subroutine calculates the  intermolecular distances, the
! spatial distribution of the atoms and the radial and pair distribution
! functions.
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

      integer, intent(IN) :: nRep      

      real (dp), dimension (:,:,:), intent(IN) :: polymer

      real(dp), optional, intent (IN) :: rmax_histo

      ! real(dp), optional, dimension(:,:,:), intent (OUT) :: density_o

      ! real(dp), optional, dimension(:,:,:), intent (OUT) :: density_h

      ! real(dp), optional, dimension(:), intent (OUT) :: rhor_o

      ! real(dp), optional, dimension(:), intent (OUT) :: rhor_h

      real(dp), optional, dimension(:), intent (OUT) :: goo

      real(dp), optional, dimension(:), intent (OUT) :: ghh

      real(dp), optional, dimension(:), intent (OUT) :: goh

!*****************************************************************************
!  local variables

      integer i, j, iRep
      integer ix, iy, iz, ir
      integer imol, jmol
      integer iat, jat

      integer iox, ih1, ih2
      integer jox, jh1, jh2      
      integer itype, jtype

      integer nbinx, nbinr

      integer nmols

      double precision, parameter :: ang_to_au = one / 0.529177d0      

      real (dp), dimension (3) :: rcm_cluster, deltar
      real (dp), dimension (3) :: voh1, voh2
      real (dp), dimension (3) :: vi, vj
      real (dp), dimension (3) :: vox, vh1, vh2, vh2o
      real (dp), dimension (3) :: voo, vhh, voh

      real(dp) :: dx, xmin, xmax, rx, ry, rz
      real(dp) :: rmax_rho, dr_rho, rmax_g, dr_g, rl, rh
      real(dp) :: voh1_sq, voh2_sq
      real(dp) :: dij
      real(dp) :: dox, dh1, dh2, dh2o
      real(dp) :: doo, doh, dhh        

      real(dp), parameter :: qh = 0.41d0
      real(dp), parameter :: qo = -0.82d0

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering PolymerProperties()'

      nmols = nAtoms/3

! If the keyword rmax_histo is present, compute densities and pair distributions

      if ( present(rmax_histo) ) then

         ! nbinx = size( density_o(1,1,:) )
         ! xmax =  rmax_histo / dsqrt(three) * ang_to_au
         ! xmin = -xmax
         ! dx = (xmax-xmin)/dfloat(nbinx)
         ! density_o = zero
         ! density_h = zero

         ! nbinr = size( rhor_o(:) )
         ! rmax_rho = rmax_histo * ang_to_au
         ! dr_rho = rmax_rho/dfloat(nbinr)         
         ! rhor_o = zero
         ! rhor_h = zero

         nbinr = size( goo )
         rmax_g = rmax_histo * ang_to_au
         dr_g = rmax_g / dfloat(nbinr)
         goo = zero
         ghh = zero
         goh = zero         

         do iRep = 1, nRep
         
            do imol = 1, nmols

               ! ! Compute density projections of the polymer on cartesian planes,
               ! ! because the 3D density would need too much memory

               ! iox = 1 + 3*(imol-1)
               ! ih1 = iox + 1
               ! ih2 = iox + 2

               ! ! O projections

               ! rx = polymer(1,iox,iRep) 
               ! ry = polymer(2,iox,iRep) 
               ! rz = polymer(3,iox,iRep)
               ! ix = floor( (rx - xmin) / dx )
               ! iy = floor( (ry - xmin) / dx )
               ! iz = floor( (rz - xmin) / dx )

               ! ! xy projection
               ! if( 0.le.ix .and. ix.le.(nbinx-1) .and. 0.le.iy .and. iy.le.(nbinx-1) ) then
               !    density_o(1,ix+1,iy+1) = density_o(1,ix+1,iy+1) + one
               ! end if
               ! ! yz projection
               ! if( 0.le.iy .and. iy.le.(nbinx-1) .and. 0.le.iz .and. iz.le.(nbinx-1) ) then
               !    density_o(2,iy+1,iz+1) = density_o(2,iy+1,iz+1) + one
               ! end if
               ! ! zx projection
               ! if( 0.le.iz .and. iz.le.(nbinx-1) .and. 0.le.ix .and. ix.le.(nbinx-1) ) then
               !    density_o(3,iz+1,ix+1) = density_o(3,iz+1,ix+1) + one
               ! end if

               ! ! H1 projections

               ! rx = polymer(1,ih1,iRep) 
               ! ry = polymer(2,ih1,iRep) 
               ! rz = polymer(3,ih1,iRep)
               ! ix = floor( (rx - xmin) / dx )
               ! iy = floor( (ry - xmin) / dx )
               ! iz = floor( (rz - xmin) / dx )

               ! ! xy projection
               ! if( 0.le.ix .and. ix.le.(nbinx-1) .and. 0.le.iy .and. iy.le.(nbinx-1) ) then
               !    density_h(1,ix+1,iy+1) = density_h(1,ix+1,iy+1) + one
               ! end if
               ! ! yz projection
               ! if( 0.le.iy .and. iy.le.(nbinx-1) .and. 0.le.iz .and. iz.le.(nbinx-1) ) then
               !    density_h(2,iy+1,iz+1) = density_h(2,iy+1,iz+1) + one
               ! end if
               ! ! zx projection
               ! if( 0.le.iz .and. iz.le.(nbinx-1) .and. 0.le.ix .and. ix.le.(nbinx-1) ) then
               !    density_h(3,iz+1,ix+1) = density_h(3,iz+1,ix+1) + one
               ! end if

               ! ! H2 projections

               ! rx = polymer(1,ih2,iRep) 
               ! ry = polymer(2,ih2,iRep) 
               ! rz = polymer(3,ih2,iRep)
               ! ix = floor( (rx - xmin) / dx )
               ! iy = floor( (ry - xmin) / dx )
               ! iz = floor( (rz - xmin) / dx )

               ! ! xy projection
               ! if( 0.le.ix .and. ix.le.(nbinx-1) .and. 0.le.iy .and. iy.le.(nbinx-1) ) then
               !    density_h(1,ix+1,iy+1) = density_h(1,ix+1,iy+1) + one
               ! end if
               ! ! yz projection
               ! if( 0.le.iy .and. iy.le.(nbinx-1) .and. 0.le.iz .and. iz.le.(nbinx-1) ) then
               !    density_h(2,iy+1,iz+1) = density_h(2,iy+1,iz+1) + one
               ! end if
               ! ! zx projection
               ! if( 0.le.iz .and. iz.le.(nbinx-1) .and. 0.le.ix .and. ix.le.(nbinx-1) ) then
               !    density_h(3,iz+1,ix+1) = density_h(3,iz+1,ix+1) + one
               ! end if
            
               ! ! Compute radial distributions of O and H
            
               ! iox = 1 + 3*(imol-1)
               ! ih1 = iox + 1
               ! ih2 = iox + 2

               ! vox(:) = polymer(:,iox,iRep)
               ! vh1(:) = polymer(:,ih1,iRep)
               ! vh2(:) = polymer(:,ih2,iRep)
            
               ! dox = dsqrt( dot_product( vox, vox ) )
               ! dh1 = dsqrt( dot_product( vh1, vh1 ) )
               ! dh2 = dsqrt( dot_product( vh2, vh2 ) )

               ! ir = floor( dox / dr_rho )
               ! if ( ir.le.(nbinr-1) ) rhor_o(ir+1) = rhor_o(ir+1) + one

               ! ir = floor( dh1 / dr_rho )
               ! if ( ir.le.(nbinr-1) ) rhor_h(ir+1) = rhor_h(ir+1) + one

               ! ir = floor( dh2 / dr_rho )
               ! if ( ir.le.(nbinr-1) ) rhor_h(ir+1) = rhor_h(ir+1) + one

               ! !do ir = 0, nbinr-1 ! This way is faster in this case
               
               ! !   rl = dfloat(ir) * dr_rho
               ! !   rh = rl + dr_rho

               ! !   if( rl.lt.dox .and. dox.lt.rh ) rhor_o(ir+1) = rhor_o(ir+1) + one
               
               ! !   if( rl.lt.dh1 .and. dh1.lt.rh ) rhor_h(ir+1) = rhor_h(ir+1) + one
               
               ! !   if( rl.lt.dh2 .and. dh2.lt.rh ) rhor_h(ir+1) = rhor_h(ir+1) + one

               ! !end do

               ! Compute the PDFs (pair distribution functions)

               iox = 1 + 3*(imol-1)
               ih1 = iox + 1
               ih2 = iox + 2

               ! intramolecular O-H1
            
               voh(:) = polymer(:,ih1,iRep) - polymer(:,iox,iRep)
               doh = dsqrt( dot_product(voh,voh) )
               j = floor( doh / dr_g )
               if ( j.le.(nbinr-1) ) goh(j+1) = goh(j+1) + one

               ! intramolecular O-H2

               voh(:) = polymer(:,ih2,iRep) - polymer(:,iox,iRep)
               doh = dsqrt( dot_product(voh,voh) )
               j = floor( doh / dr_g )
               if ( j.le.(nbinr-1) ) goh(j+1) = goh(j+1) + one                        

               ! intramolecular H-H
               
               vhh(:) = polymer(:,ih2,iRep) - polymer(:,ih1,iRep)
               dhh = dsqrt( dot_product(vhh,vhh) )
               j = floor( dhh / dr_g )
               if ( j.le.(nbinr-1) ) ghh(j+1) = ghh(j+1) + one

               do jmol = imol+1, nmols
            
                  jox = 1 + 3*(jmol-1)
                  jh1 = jox + 1
                  jh2 = jox + 2

                  ! intermolecular O(imol)-O(jmol)
                  voo(:) = polymer(:,jox,iRep) - polymer(:,iox,iRep)
                  doo = dsqrt( dot_product(voo,voo) )
                  j = floor( doo / dr_g )
                  if ( j.le.(nbinr-1) ) goo(j+1) = goo(j+1) + one

                  ! intermolecular O(imol)-H1(jmol)

                  voh(:) = polymer(:,jh1,iRep) - polymer(:,iox,iRep)
                  doh = dsqrt( dot_product(voh,voh) )
                  j = floor( doh / dr_g )
                  if ( j.le.(nbinr-1) ) goh(j+1) = goh(j+1) + one               

                  ! intermolecular O(imol)-H2(jmol)

                  voh(:) = polymer(:,jh2,iRep) - polymer(:,iox,iRep)
                  doh = dsqrt( dot_product(voh,voh) )
                  j = floor( doh / dr_g )
                  if ( j.le.(nbinr-1) ) goh(j+1) = goh(j+1) + one               

                  ! intermolecular H1(imol)-H1(jmol)
                  
                  vhh(:) = polymer(:,jh1,iRep) - polymer(:,ih1,iRep)
                  dhh = dsqrt( dot_product(vhh,vhh) )
                  j = floor( dhh / dr_g )
                  if ( j.le.(nbinr-1) ) ghh(j+1) = ghh(j+1) + one               

                  ! intermolecular H1(imol)-H2(jmol)

                  vhh(:) = polymer(:,jh2,iRep) - polymer(:,ih1,iRep)
                  dhh = dsqrt( dot_product(vhh,vhh) )
                  j = floor( dhh / dr_g )
                  if ( j.le.(nbinr-1) ) ghh(j+1) = ghh(j+1) + one               

                  ! intermolecular H2(imol)-H1(jmol)

                  vhh(:) = polymer(:,jh1,iRep) - polymer(:,ih2,iRep)
                  dhh = dsqrt( dot_product(vhh,vhh) )
                  j = floor( dhh / dr_g )
                  if ( j.le.(nbinr-1) ) ghh(j+1) = ghh(j+1) + one               

                  ! intermolecular H2(imol)-H2(jmol)
                  
                  vhh(:) = polymer(:,jh2,iRep) - polymer(:,ih2,iRep)
                  dhh = dsqrt( dot_product(vhh,vhh) )
                  j = floor( dhh / dr_g )
                  if ( j.le.(nbinr-1) ) ghh(j+1) = ghh(j+1) + one               
               
               end do ! end jmol            

            end do ! end imol

         end do ! end iRep

      end if ! end present(rmax_histo)

      if ( debug ) write(*,*) 'Exiting PolymerProperties()'
      
end subroutine PolymerProperties
