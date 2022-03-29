!*****************************************************************************
module SPCF
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Created on: Feb 2020 by Alfonso Gijón
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This module provides a public subroutine, evaluate_energy_SPCF, which
! evaluates the energy and force of a water cluster system subject to a
! flexible SPCF potential, given the positions of the atoms
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  private module variables

      ! constants to convert parameters of SPCF potential

      double precision, private, parameter :: zero = 0.0d0
      double precision, private, parameter :: one = 1.0d0
      double precision, private, parameter :: two = 2.0d0
      double precision, private, parameter :: four = 4.0d0
      double precision, private, parameter :: six = 6.0d0
      double precision, private, parameter :: twelve = 12.0d0
      double precision, private, parameter :: pi = 3.14159265358979323846264338327950288419716939937510d0
      double precision, private, parameter :: kelvin_to_au = 3.16683E-6
      double precision, private, parameter :: mdyn_to_au = one / 8.2388584d0
      double precision, private, parameter :: mdyntimesang_to_au = one / 4.3598149d0
      double precision, private, parameter :: mdynoverang_to_au = one / 15.5691905d0
      double precision, private, parameter :: kcalmol_to_au = one / 627.510d0
      double precision, private, parameter :: ang_to_au = one / 0.529177d0
      double precision, private, parameter :: deg_to_rad = pi / 180.0d0
      double precision, private, parameter :: eV = one / 13.60580d0 / two
      double precision, private, parameter :: Kelvin = eV / 11604.45d0

      ! Parameters of the potential
      
      double precision, private :: epsilon, sigma, qh, qo
      double precision, private :: rhow, dw, roh_eq, rhh_eq, theta_eq, bc, cc, dc
      double precision, private :: alj, blj
      double precision, private :: kw, kb, kc, kd
      !double precision, dimension(3) :: charge


!*****************************************************************************
!  private module subroutines

contains

!*****************************************************************************
subroutine SPCF_setup()
!*****************************************************************************
!
!  Project: MolecularDynamics    
!
!  Created on: Sep 2019 by Alfonso Gijón
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
! 
! Define the parameters of the SPCF potential as in J. Chem. Phys 106, 2400 (1997)  
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  variables

!*****************************************************************************
!  local variables

!*****************************************************************************
!  start of subroutine

      ! potential parameters
      
      epsilon = 78.22d0       ! kelvins
      sigma = 3.165d0         ! Angstrom
      qh = 0.41d0             ! electron charge
      qo = -two * qh           ! electron charge

      ! SPC/F
      rhow = 2.566d0          ! Angstroms^-1. Bug in J. Chem. Phys., Vol 106, No. 6, 1997 !!!
      dw = 0.708d0            ! mdyn * Angstrom
      roh_eq = 1.0d0          ! Angstroms
      theta_eq = 109.47d0     ! º
      bc = 2.283d0            ! mdyn / Angstrom
      cc = -1.469d0           ! mdyn / Angstrom
      dc = 0.776d0            ! mdyn / Angstrom      

      ! SPC/F2
      !rhow = 2.361d0          ! Angstroms^-1. Bug in J. Chem. Phys., Vol 106, No. 6, 1997 !!!
      !dw = 0.708d0            ! mdyn * Angstrom
      !roh_eq = 1.0d0          ! Angstroms
      !theta_eq = 108.0d0      ! º
      !bc = 1.803d0            ! mdyn / Angstrom
      !cc = -1.469d0           ! mdyn / Angstrom
      !dc = 0.776d0            ! mdyn / Angstrom

      !charge(1) = qo
      !charge(2) = qh
      !charge(3) = qh

      ! convert parameters to atomic units

      epsilon = epsilon * Kelvin
      sigma = sigma * ang_to_au
      rhow = rhow / ang_to_au
      dw = dw * mdyntimesang_to_au
      roh_eq = roh_eq * ang_to_au
      theta_eq = theta_eq * deg_to_rad
      bc = bc * mdynoverang_to_au
      cc = cc * mdynoverang_to_au
      dc = dc * mdynoverang_to_au

      ! test
      !write(*,*)
      !write(*,*) 'SPC/F parameters in a.u.'
      !write(*,*) 'epsilon', epsilon
      !write(*,*) 'sigma', sigma
      !write(*,*) 'rhow^2 * Dw', rhow * rhow * dw
      !write(*,*) 'b/2', bc / 2.0
      !write(*,*) 'c', cc
      !write(*,*) 'd', dc
      !stop
      
      alj = four * epsilon * sigma ** twelve
      blj = four * epsilon * sigma ** six
      rhh_eq = two * roh_eq * dsin( theta_eq / two )
      kw = dw*rhow*rhow
      kb = bc/two
      kc = cc
      kd = dc

      ! test
      !print*, 'r_OH'
      !print*, sqrt( (rhh_eq/two)**two + (roh_eq*dcos(theta_eq/two))**two ) / ang_to_au
      !stop
      
end subroutine SPCF_setup

!*****************************************************************************
subroutine evaluate_energy_SPCF( nmols, position, energy, force )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Created on: Sep 2019 by Alfonso Gijón
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! Given the position of the atoms, calculate energy and force  
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  variables

      integer, intent (IN) :: nmols
      
      double precision, dimension(:,:), intent (IN) :: position

      double precision, intent (OUT) :: energy
      double precision, dimension(:,:), intent (OUT) :: force

!*****************************************************************************
!  local variables

      integer :: imol, jmol
      integer :: iox, ih1, ih2, jox

      integer :: k, l
      integer :: ik, il
      
      double precision :: uinter, uintra
      double precision :: ulj, ucoulomb
      double precision :: ustretch, ubend
      
      double precision :: rxij, ryij, rzij
      double precision :: rij, rij2, rij6, rij7, rij_inv
      double precision :: cosx, cosy, cosz
      double precision :: dV
      
      double precision :: qk, ql
      double precision :: rxkl, rykl, rzkl
      double precision :: rkl, rkl2, rkl_inv

      double precision :: rxoh1, ryoh1, rzoh1, rxoh2, ryoh2, rzoh2, rxhh, ryhh, rzhh
      double precision :: roh1, roh2, rhh, roh1_inv, roh2_inv, rhh_inv
      
      double precision :: vw1, dvw1, vw2, dvw2, vb, dvb, vc, dvc_droh, dvc_drhh
      double precision :: vd, dvd_droh1, dvd_droh2
      double precision :: deq, deq2

      double precision, dimension(3) :: rk, rl
      double precision, dimension(3) :: cosij, coskl, dVcos, cosoh1, cosoh2, coshh
      
!*****************************************************************************
!  start subroutine

      
      
      ! calculate energy and force

      ustretch = zero
      ubend = zero
      ulj = zero
      ucoulomb = zero
      force = zero

      !$OMP PARALLEL DO &
      !$OMP PRIVATE( imol, iox, ih1, ih2 ) &
      !$OMP PRIVATE( rxoh1, ryoh1, rzoh1, rxoh2, ryoh2, rzoh2, rxhh, ryhh, rzhh ) &
      !$OMP PRIVATE( roh1, roh2, rhh, roh1_inv, roh2_inv, rhh_inv ) &            
      !$OMP PRIVATE( cosoh1, cosoh2, deq, deq2, vw1, dvw1, vw2, dvw2 ) &
      !$OMP PRIVATE( coshh, vb, dvb, vc, dvc_droh, dvc_drhh, vd, dvd_droh1, dvd_droh2 ) &
      !$OMP PRIVATE( jmol, jox, rxij, ryij, rzij, rij2, rij, rij6, rij7, rij_inv ) &
      !$OMP PRIVATE( cosx, cosy, cosz, cosij, coskl, dV, dVcos ) &
      !$OMP PRIVATE( k, ik, rk, qk, l, il, rl, ql, rxkl, rykl, rzkl, rkl2, rkl, rkl_inv) &
      !$OMP REDUCTION( + : ustretch, ubend, ulj, ucoulomb, force ) 
      do imol = 1, nmols

         iox = 1 + 3*(imol-1)
         ih1 = iox + 1
         ih2 = iox + 2

         ! intra-molecular energy

         rxoh1 =  position(1,ih1) - position(1,iox)
         ryoh1 =  position(2,ih1) - position(2,iox)
         rzoh1 =  position(3,ih1) - position(3,iox)
         rxoh2 =  position(1,ih2) - position(1,iox)
         ryoh2 =  position(2,ih2) - position(2,iox)
         rzoh2 =  position(3,ih2) - position(3,iox)
         rxhh =  position(1,ih2) - position(1,ih1)
         ryhh =  position(2,ih2) - position(2,ih1)
         rzhh =  position(3,ih2) - position(3,ih1)         

         roh1 = sqrt( rxoh1*rxoh1 + ryoh1*ryoh1 + rzoh1*rzoh1 )
         roh2 = sqrt( rxoh2*rxoh2 + ryoh2*ryoh2 + rzoh2*rzoh2 )
         rhh = sqrt( rxhh*rxhh + ryhh*ryhh + rzhh*rzhh )
         roh1_inv = one / roh1
         roh2_inv = one / roh2
         rhh_inv = one / rhh
         
         ! O-H streching
         
         cosoh1 = (/ rxoh1 * roh1_inv, ryoh1 * roh1_inv, rzoh1 * roh1_inv /)
         cosoh2 = (/ rxoh2 * roh2_inv, ryoh2 * roh2_inv, rzoh2 * roh2_inv /)

         deq = roh1 - roh_eq
         deq2 = deq * deq
         vw1 = kw * deq2
         dvw1 = two * kw * deq
         
         deq = roh2 - roh_eq
         deq2 = deq * deq
         vw2 = kw * deq2
         dvw2 = two * kw * deq

         ustretch = ustretch + vw1 + vw2
         
         force(1:3,iox) = force(1:3,iox) + dvw1 * cosoh1 + dvw2 * cosoh2
         force(1:3,ih1) = force(1:3,ih1) - dvw1 * cosoh1
         force(1:3,ih2) = force(1:3,ih2) - dvw2 * cosoh2

         ! H-O-H bending is composed by three terms: Vb, Vc and Vd terms
         
         ! Vb term: H-H stretching
         
         coshh = (/ rxhh * rhh_inv, ryhh * rhh_inv, rzhh * rhh_inv /)
         
         deq = rhh - rhh_eq
         deq2 = deq * deq
         vb = kb * deq2
         dvb = two * kb * deq         

         ustretch = ustretch + vb

         force(1:3,ih1) = force(1:3,ih1) + dvb * coshh
         force(1:3,ih2) = force(1:3,ih2) - dvb * coshh

         ! Vc term

         vc = kc * ( roh1 + roh2 - two*roh_eq ) * ( rhh - rhh_eq )
         dvc_droh = kc * ( rhh - rhh_eq )
         dvc_drhh = kc * ( roh1 + roh2 - two*roh_eq )
         
         ubend = ubend + vc

         force(1:3,iox) = force(1:3,iox) + dvc_droh * cosoh1 + &
                                           dvc_droh * cosoh2 
                                           
         force(1:3,ih1) = force(1:3,ih1) - dvc_droh * cosoh1 + dvc_drhh * coshh
         force(1:3,ih2) = force(1:3,ih2) - dvc_droh * cosoh2 - dvc_drhh * coshh

         ! Vd term

         vd = kd * ( roh1 - roh_eq ) * ( roh2 - roh_eq )
         dvd_droh1 = kd * ( roh2 - roh_eq )
         dvd_droh2 = kd * ( roh1 - roh_eq )         
         
         ubend = ubend + vd

         force(1:3,iox) = force(1:3,iox) + dvd_droh1 * cosoh1 + &
                                           dvd_droh2 * cosoh2 
                                           
         force(1:3,ih1) = force(1:3,ih1) - dvd_droh1 * cosoh1
         force(1:3,ih2) = force(1:3,ih2) - dvd_droh2 * cosoh2
         
         ! inter-molecular energy
         
         do jmol = imol+1, nmols

            jox = 1 + 3*(jmol-1)
            
            rxij = position(1,jox) - position(1,iox)
            ryij = position(2,jox) - position(2,iox)
            rzij = position(3,jox) - position(3,iox)
            
            rij2 = rxij*rxij + ryij*ryij + rzij*rzij
            rij = sqrt(rij2)
            rij6 = rij2 * rij2 *rij2
            rij7 = rij6 * rij
            rij_inv = one / rij            
            
            cosx = rxij * rij_inv
            cosy = ryij * rij_inv
            cosz = rzij * rij_inv

            cosij = (/ cosx, cosy, cosz /)

            ! oxygen - oxygen LJ interaction
            
            ulj = ulj + ( alj/rij6 - blj ) / rij6

            dV = (-twelve * alj / rij6 + six * blj ) / rij7

            force(1:3,iox) = force(1:3,iox) + dV * cosij

            force(1:3,jox) = force(1:3,jox) - dV * cosij

            ! Coulomb interaction between charges k and l of imol and jmol
            
            do k = 1, 3 ! k-charge in imol

               ik = k + 3*(imol-1)
               rk = position( :, k+3*(imol-1) )

               !qk = charge(k)
               
               if ( k.eq.1 ) then ! oxygen
                  qk = qo
               else               ! hydrogen
                  qk = qh
               end if
               
               do l = 1, 3 ! l-charge in jmol

                  il = l + 3*(jmol-1)
                  rl = position( :, l+3*(jmol-1) )

                  !ql = charge(l)
                  
                  if ( l.eq.1 ) then ! oxygen
                     ql = qo
                  else               ! hydrogen
                     ql = qh
                  end if

                  rxkl = rl(1) - rk(1)
                  rykl = rl(2) - rk(2)
                  rzkl = rl(3) - rk(3)
                  
                  rkl2 = rxkl*rxkl + rykl*rykl + rzkl*rzkl
                  rkl = sqrt(rkl2)
                  rkl_inv = one / rkl

                  cosx = rxkl * rkl_inv
                  cosy = rykl * rkl_inv
                  cosz = rzkl * rkl_inv

                  coskl = (/ cosx, cosy, cosz /)

                  ucoulomb = ucoulomb + qk * ql / rkl

                  ! Coulomb forces
                  
                  dV = - qk * ql / rkl2

                  dVcos = dV * coskl

                  force(1:3,ik) = force(1:3,ik) + dVcos
                  force(1:3,il) = force(1:3,il) - dVcos
                     
               end do ! end l

            end do ! end k

         end do ! end jmol
         
      end do ! end imol
      !$OMP END PARALLEL DO      
      
      uinter = ulj + ucoulomb
      uintra = ustretch + ubend
      energy = uinter + uintra

      ! print*, 
      ! print*, 'Uintra, Uinter'
      ! print*, uintra, uinter
      ! stop

end subroutine evaluate_energy_SPCF

end module SPCF
