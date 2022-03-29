!*****************************************************************************
subroutine CentroidProperties_dipole( nAtoms, mass, centroid, radius2, radmax2, roh2, theta, &
     rmax_histo, density_o, density_h, rhor_o, rhor_h, rhor_h2o, goo, ghh, goh, &
     mutot, mumod)

! subroutine CentroidProperties( nAtoms, mass, centroid, radius2, roh2, theta, &
!       rmax_histo, rhor_o, rhor_h, rhor_h2o, goo, ghh, goh )  
  
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
! This subroutine calculates the radius of the centroid of the cluster, with
! respect to the center of masses. Also the intermolecular distances, the
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

      real (dp), dimension (:), intent(IN) :: mass

      real (dp), dimension (:,:), intent(IN) :: centroid

      real(dp), intent (OUT) :: radius2

      real(dp), intent (OUT) :: radmax2      

      real(dp), intent (OUT) :: roh2

      real(dp), intent (OUT) :: theta

      real(dp), optional, intent (IN) :: rmax_histo

      real(dp), optional, dimension(:,:,:), intent (OUT) :: density_o

      real(dp), optional, dimension(:,:,:), intent (OUT) :: density_h

      real(dp), optional, dimension(:), intent (OUT) :: rhor_o

      real(dp), optional, dimension(:), intent (OUT) :: rhor_h

      real(dp), optional, dimension(:), intent (OUT) :: rhor_h2o      

      real(dp), optional, dimension(:), intent (OUT) :: goo

      real(dp), optional, dimension(:), intent (OUT) :: ghh

      real(dp), optional, dimension(:), intent (OUT) :: goh

      real(dp), optional, dimension(3), intent (OUT) :: mutot

      real(dp), optional, intent (OUT) :: mumod

!*****************************************************************************
!  local variables

      integer i, j
      integer ix, iy, iz, ir
      integer imol, jmol
      integer iat, jat

      integer iox, ih1, ih2
      integer jox, jh1, jh2      
      integer itype, jtype

      integer nbinx, nbinr
      
      integer nmols

      double precision, parameter :: ang_to_au = one / 0.529177d0      

      real(dp) :: mtotal, mh2o
      real (dp), dimension (3) :: rcm_cluster, deltar
      real (dp), dimension (3) :: voh1, voh2
      real (dp), dimension (3) :: vi, vj
      real (dp), dimension (3) :: vox, vh1, vh2, vh2o
      real (dp), dimension (3) :: voo, vhh, voh

      real (dp), allocatable, dimension(:,:) :: rcm_mol

      real(dp) :: dx, xmin, xmax, rx, ry, rz
      real(dp) :: rmax_rho, dr_rho, rmax_g, dr_g, rl, rh
      real(dp) :: voh1_sq, voh2_sq
      real(dp) :: dij
      real(dp) :: dox, dh1, dh2, dh2o
      real(dp) :: doo, doh, dhh

      double precision :: r2_aux, r2_max

      real(dp), parameter :: qh = 0.41d0
      real(dp), parameter :: qo = -0.82d0

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering CentroidProperties()'

      nmols = nAtoms/3
      
      allocate( rcm_mol( 3, nmols ) )

! radius2
      
      ! center of masses of the centroid cluster
         
      rcm_cluster = zero
      mtotal = zero
      do j = 1, nAtoms
         rcm_cluster = rcm_cluster + mass(j)*centroid(:,j)
         mtotal = mtotal + mass(j)
      end do
      rcm_cluster = rcm_cluster / mtotal

      !print*, ""
      !print*, "CM", rcm_cluster
      !print*, ""

      ! mean distance of water molecules to the center of masses of the cluster

      radius2 = zero
      mh2o = mass(1) + two * mass(2)

      rcm_mol = zero

      r2_max = zero
      
      do imol = 1, nmols
         
         do j = 1, 3
            iat = j + 3*(imol-1)
            rcm_mol(:,imol) = rcm_mol(:,imol) + mass(iat)*centroid(:,iat)
         end do
         rcm_mol(:,imol) = rcm_mol(:,imol) / mh2o

         deltar = rcm_mol(:,imol) - rcm_cluster
         r2_aux = dot_product( deltar, deltar )
         radius2 = radius2 + r2_aux

         if (r2_aux.gt.r2_max) r2_max = r2_aux
            
      end do
      radius2 = radius2 / dfloat( nmols )
      radmax2 = r2_max

! intramolecular properties: OH distance and HOH angle 

      roh2 = zero
      theta = zero
      
      do imol = 1, nmols
         
         iox = 1 + 3*(imol-1)
         ih1 = iox + 1
         ih2 = iox + 2

         ! roh and theta
         
         voh1 = centroid(:,ih1) - centroid(:,iox)
         voh1_sq = dot_product(voh1,voh1)
         roh2 = roh2 + voh1_sq
         voh2 = centroid(:,ih2) - centroid(:,iox)
         voh2_sq = dot_product(voh2,voh2)
         roh2 = roh2 + voh2_sq

         theta = theta + dacos( dot_product(voh1,voh2) / ( dsqrt(voh1_sq) * dsqrt(voh2_sq) ) )
         
      end do

      roh2 = roh2 * half / dfloat( nmols )
      theta = theta / dfloat( nmols )            

! If the keyword rmax_histo is present, compute densities and pair distributions

      if ( present(rmax_histo) ) then

         nbinx = size( density_o(1,1,:) )
         xmax =  rmax_histo / dsqrt(three) * ang_to_au
         xmin = -xmax
         dx = (xmax-xmin)/dfloat(nbinx)
         density_o = zero
         density_h = zero

         nbinr = size( rhor_o(:) )
         rmax_rho = rmax_histo * ang_to_au
         dr_rho = rmax_rho/dfloat(nbinr)         
         rhor_o = zero
         rhor_h = zero
         rhor_h2o = zero

         rmax_g = rmax_histo * ang_to_au
         dr_g = rmax_g / dfloat(nbinr)
         goo = zero
         ghh = zero
         goh = zero         

         mutot = zero
         
         !!$OMP PARALLEL DO &
         !!$OMP PRIVATE( imol, iox, ih1, ih2, rx, ry, rz, ix, iy, iz ) &
         !!$OMP PRIVATE( vox, vh1, vh2, vh2o, dox, dh1, dh2, dh2o, ir ) &
         !!$OMP PRIVATE( voh, doh, j, vhh, dhh, jmol, jox, jh1, jh2, voo, doo ) &
         !!$OMP REDUCTION( + : density_o, density_h, rhor_o, rhor_h, rhor_h2o, goo, goh, ghh )
         do imol = 1, nmols

         ! Compute density projections of the centroid on cartesian planes, because the 3D density
         ! would need too much memory

            iox = 1 + 3*(imol-1)
            ih1 = iox + 1
            ih2 = iox + 2

            ! O projections

            rx = centroid(1,iox) 
            ry = centroid(2,iox) 
            rz = centroid(3,iox)
            ix = floor( (rx - xmin) / dx )
            iy = floor( (ry - xmin) / dx )
            iz = floor( (rz - xmin) / dx )

            ! xy projection
            if( 0.le.ix .and. ix.le.(nbinx-1) .and. 0.le.iy .and. iy.le.(nbinx-1) ) then
               density_o(1,ix+1,iy+1) = density_o(1,ix+1,iy+1) + one
            end if
            ! yz projection
            if( 0.le.iy .and. iy.le.(nbinx-1) .and. 0.le.iz .and. iz.le.(nbinx-1) ) then
               density_o(2,iy+1,iz+1) = density_o(2,iy+1,iz+1) + one
            end if
            ! zx projection
            if( 0.le.iz .and. iz.le.(nbinx-1) .and. 0.le.ix .and. ix.le.(nbinx-1) ) then
               density_o(3,iz+1,ix+1) = density_o(3,iz+1,ix+1) + one
            end if                        

            ! H1 projections

            rx = centroid(1,ih1) 
            ry = centroid(2,ih1) 
            rz = centroid(3,ih1)
            ix = floor( (rx - xmin) / dx )
            iy = floor( (ry - xmin) / dx )
            iz = floor( (rz - xmin) / dx )

            ! xy projection
            if( 0.le.ix .and. ix.le.(nbinx-1) .and. 0.le.iy .and. iy.le.(nbinx-1) ) then
               density_h(1,ix+1,iy+1) = density_h(1,ix+1,iy+1) + one
            end if
            ! yz projection
            if( 0.le.iy .and. iy.le.(nbinx-1) .and. 0.le.iz .and. iz.le.(nbinx-1) ) then
               density_h(2,iy+1,iz+1) = density_h(2,iy+1,iz+1) + one
            end if
            ! zx projection
            if( 0.le.iz .and. iz.le.(nbinx-1) .and. 0.le.ix .and. ix.le.(nbinx-1) ) then
               density_h(3,iz+1,ix+1) = density_h(3,iz+1,ix+1) + one
            end if                                    

            ! H2 projections

            rx = centroid(1,ih2) 
            ry = centroid(2,ih2) 
            rz = centroid(3,ih2)
            ix = floor( (rx - xmin) / dx )
            iy = floor( (ry - xmin) / dx )
            iz = floor( (rz - xmin) / dx )

            ! xy projection
            if( 0.le.ix .and. ix.le.(nbinx-1) .and. 0.le.iy .and. iy.le.(nbinx-1) ) then
               density_h(1,ix+1,iy+1) = density_h(1,ix+1,iy+1) + one
            end if
            ! yz projection
            if( 0.le.iy .and. iy.le.(nbinx-1) .and. 0.le.iz .and. iz.le.(nbinx-1) ) then
               density_h(2,iy+1,iz+1) = density_h(2,iy+1,iz+1) + one
            end if
            ! zx projection
            if( 0.le.iz .and. iz.le.(nbinx-1) .and. 0.le.ix .and. ix.le.(nbinx-1) ) then
               density_h(3,iz+1,ix+1) = density_h(3,iz+1,ix+1) + one
            end if                                                
            
         ! Compute radial distributions of O and H
            
            iox = 1 + 3*(imol-1)
            ih1 = iox + 1
            ih2 = iox + 2

            vox(:) = centroid(:,iox) - rcm_cluster(:)
            vh1(:) = centroid(:,ih1) - rcm_cluster(:)
            vh2(:) = centroid(:,ih2) - rcm_cluster(:)
            vh2o(:) = rcm_mol(:,imol) - rcm_cluster(:)

            ! Compute total dipole of the cluster

            mutot = mutot + qo * vox + qh * vh1 + qh * vh2
            
            dox = dsqrt( dot_product( vox, vox ) )
            dh1 = dsqrt( dot_product( vh1, vh1 ) )
            dh2 = dsqrt( dot_product( vh2, vh2 ) )
            dh2o = dsqrt( dot_product( vh2o, vh2o ) )

            ir = floor( dox / dr_rho )
            if ( ir.le.(nbinr-1) ) rhor_o(ir+1) = rhor_o(ir+1) + one

            ir = floor( dh1 / dr_rho )
            if ( ir.le.(nbinr-1) ) rhor_h(ir+1) = rhor_h(ir+1) + one

            ir = floor( dh2 / dr_rho )
            if ( ir.le.(nbinr-1) ) rhor_h(ir+1) = rhor_h(ir+1) + one

            ir = floor( dh2o / dr_rho )
            if ( ir.le.(nbinr-1) ) rhor_h2o(ir+1) = rhor_h2o(ir+1) + one                        
            
            !do ir = 0, nbinr-1 ! This way is faster in this case
               
            !   rl = dfloat(ir) * dr_rho
            !   rh = rl + dr_rho

            !   if( rl.lt.dox .and. dox.lt.rh ) rhor_o(ir+1) = rhor_o(ir+1) + one

            !   if( rl.lt.dh1 .and. dh1.lt.rh ) rhor_h(ir+1) = rhor_h(ir+1) + one

            !   if( rl.lt.dh2 .and. dh2.lt.rh ) rhor_h(ir+1) = rhor_h(ir+1) + one

            !   if( rl.lt.dh2o .and. dh2o.lt.rh ) rhor_h2o(ir+1) = rhor_h2o(ir+1) + one
            
            !end do

         ! Compute the PDFs (pair distribution functions)

            iox = 1 + 3*(imol-1)
            ih1 = iox + 1
            ih2 = iox + 2

            ! intramolecular O-H1
            
            voh(:) = centroid(:,ih1) - centroid(:,iox)
            doh = dsqrt( dot_product(voh,voh) )
            j = floor( doh / dr_g )
            if ( j.le.(nbinr-1) ) goh(j+1) = goh(j+1) + one

            ! intramolecular O-H2

            voh(:) = centroid(:,ih2) - centroid(:,iox)
            doh = dsqrt( dot_product(voh,voh) )
            j = floor( doh / dr_g )
            if ( j.le.(nbinr-1) ) goh(j+1) = goh(j+1) + one                        

            ! intramolecular H-H

            vhh(:) = centroid(:,ih2) - centroid(:,ih1)
            dhh = dsqrt( dot_product(vhh,vhh) )
            j = floor( dhh / dr_g )
            if ( j.le.(nbinr-1) ) ghh(j+1) = ghh(j+1) + one

            do jmol = imol+1, nmols
            
               jox = 1 + 3*(jmol-1)
               jh1 = jox + 1
               jh2 = jox + 2

               ! intermolecular O(imol)-O(jmol)
               voo(:) = centroid(:,jox) - centroid(:,iox)
               doo = dsqrt( dot_product(voo,voo) )
               j = floor( doo / dr_g )
               if ( j.le.(nbinr-1) ) goo(j+1) = goo(j+1) + one

               ! intermolecular O(imol)-H1(jmol)

               voh(:) = centroid(:,jh1) - centroid(:,iox)
               doh = dsqrt( dot_product(voh,voh) )
               j = floor( doh / dr_g )
               if ( j.le.(nbinr-1) ) goh(j+1) = goh(j+1) + one               

               ! intermolecular O(imol)-H2(jmol)

               voh(:) = centroid(:,jh2) - centroid(:,iox)
               doh = dsqrt( dot_product(voh,voh) )
               j = floor( doh / dr_g )
               if ( j.le.(nbinr-1) ) goh(j+1) = goh(j+1) + one               

               ! intermolecular H1(imol)-H1(jmol)

               vhh(:) = centroid(:,jh1) - centroid(:,ih1)
               dhh = dsqrt( dot_product(vhh,vhh) )
               j = floor( dhh / dr_g )
               if ( j.le.(nbinr-1) ) ghh(j+1) = ghh(j+1) + one               

               ! intermolecular H1(imol)-H2(jmol)

               vhh(:) = centroid(:,jh2) - centroid(:,ih1)
               dhh = dsqrt( dot_product(vhh,vhh) )
               j = floor( dhh / dr_g )
               if ( j.le.(nbinr-1) ) ghh(j+1) = ghh(j+1) + one               

               ! intermolecular H2(imol)-H1(jmol)

               vhh(:) = centroid(:,jh1) - centroid(:,ih2)
               dhh = dsqrt( dot_product(vhh,vhh) )
               j = floor( dhh / dr_g )
               if ( j.le.(nbinr-1) ) ghh(j+1) = ghh(j+1) + one               

               ! intermolecular H2(imol)-H2(jmol)

               vhh(:) = centroid(:,jh2) - centroid(:,ih2)
               dhh = dsqrt( dot_product(vhh,vhh) )
               j = floor( dhh / dr_g )
               if ( j.le.(nbinr-1) ) ghh(j+1) = ghh(j+1) + one               
               
            end do ! end jmol            

         end do ! end imol
         !!$OMP END PARALLEL DO
         mumod = dsqrt( dot_product( mutot, mutot ) )
         
      end if ! end present(rmax_histo)

      deallocate( rcm_mol )      
      
      if ( debug ) write(*,*) 'Exiting CentroidProperties()'
      
end subroutine CentroidProperties_dipole
