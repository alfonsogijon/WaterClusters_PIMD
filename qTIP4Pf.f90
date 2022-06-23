!*****************************************************************************
module qTIP4Pf
!*****************************************************************************
!
!  Project: MolecularDynamics
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
! This module provides a public subroutine, evaluate_energy_qTIP4Pf, which
! evaluates the energy and force of a water cluster system subject to a qTIP4Pf
! potential, given the positions of the atoms
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  private module variables

! constants to convert parameters of qTIP4Pf potential

      double precision, private, parameter :: zero = 0.0d0
      double precision, private, parameter :: one = 1.0d0
      double precision, private, parameter :: two = 2.0d0
      double precision, private, parameter :: three = 3.0d0
      double precision, private, parameter :: four = 4.0d0
      double precision, private, parameter :: five = 5.0d0
      double precision, private, parameter :: six = 6.0d0
      double precision, private, parameter :: seven = 7.0d0
      double precision, private, parameter :: eight = 8.0d0
      double precision, private, parameter :: nine = 9.0d0
      double precision, private, parameter :: ten = 10.0d0
      double precision, private, parameter :: eleven = 11.0d0
      double precision, private, parameter :: twelve = 12.0d0
      double precision, private, parameter :: half = 0.5d0
      double precision, private, parameter :: sixth = 1.0d0 / 6.0d0
      double precision, private, parameter :: pi = 3.14159265358979323846264338327950288419716939937510d0

      double precision, parameter :: kcalmol_to_au_p = one / 627.510d0
      double precision, parameter :: ang_to_au_p = one / 0.529177d0
      double precision, parameter :: deg_to_rad = pi / 180.0d0

      double precision, private :: epsilon, sigma, qh, qm
      double precision, private :: gammao, gammah, dr, alphar, req, ktheta, thetaeq
      double precision, private :: alj, blj, alphar2, alphar3, alphar4

!*****************************************************************************
!  private module subroutines

      private voh
      private dvoh

contains

!*****************************************************************************
subroutine qTIP4Pf_setup()
!*****************************************************************************
!
!  Project: MolecularDynamics    
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
! Define the parameters of the qTIP4Pf potential as in J. Chem. Phys 131, 024501 (2009)  
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
      
      epsilon = 0.1852d0       ! kcal/mol
      sigma = 3.1589d0         ! Angstrom
      qh = 0.5564d0            ! electron charge
      qm = -two * qh            ! electron charge
      gammao = 0.73612d0
      dr = 116.09d0            ! kcal/mol
      alphar = 2.287d0         ! 1 * Angstrom^-1
      req = 0.9419d0           ! Angstrom
      ktheta = 87.85d0         ! kcal/mol * rad^-2
      thetaeq = 107.4d0        ! º

      ! convert parameters to atomic units

      epsilon = epsilon * kcalmol_to_au_p
      sigma = sigma * ang_to_au_p
      dr = dr * kcalmol_to_au_p
      alphar = alphar / ang_to_au_p
      req = req * ang_to_au_p
      ktheta = ktheta* kcalmol_to_au_p
      thetaeq = thetaeq * deg_to_rad

      alj = four * epsilon * sigma ** twelve
      blj = four * epsilon * sigma ** six
      alphar2 = alphar * alphar
      alphar3 = alphar2 * alphar
      alphar4 = alphar2 * alphar2
      gammah = (one - gammao) / two
      
end subroutine qTIP4Pf_setup

!*****************************************************************************
subroutine evaluate_energy_qTIP4Pf( nmols, position, energy, force )
!*****************************************************************************
!
!  Project: MolecularDynamics
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
      integer :: iat, jat
      integer :: iox, ih1, ih2, jox

      integer :: k, l
      integer :: ik, il
      
      double precision :: uinter, uintra
      double precision :: ulj, ucoulomb
      double precision :: ustrech, ubend
      
      double precision :: rxij, ryij, rzij
      double precision :: rij, rij2, rij6, rij7
      double precision :: cosx, cosy, cosz
      double precision :: dV
      
      double precision :: qk, ql
      double precision :: rxkl, rykl, rzkl
      double precision :: rkl, rkl2

      double precision :: rxoh1, ryoh1, rzoh1, rxoh2, ryoh2, rzoh2
      double precision :: roh1, roh2
      double precision :: theta, dtheta, dtheta2
      double precision :: costheta, ftheta

      double precision, dimension(3) :: rk, rl
      double precision, dimension(3) :: cosij, coskl, dVcos, cosoh1, cosoh2
      double precision, allocatable, dimension(:,:) :: rM
      
!*****************************************************************************
!  start subroutine

      ! construct array with positions of M-points

      allocate( rM( 3,nmols ) )

      do imol = 1, nmols

         iox = 1 + 3 * (imol - 1)
         ih1 = iox + 1
         ih2 = iox + 2
         
         rM(:,imol) = gammao * position(:,iox) + &
              gammah * ( position(:,ih1) + position(:,ih2) )
         
      end do

      ! calculate energy and force

      ustrech = zero
      ubend = zero
      ulj = zero
      ucoulomb = zero
      force = zero

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

         roh1 = dsqrt( rxoh1*rxoh1 + ryoh1*ryoh1 + rzoh1*rzoh1 )
         roh2 = dsqrt( rxoh2*rxoh2 + ryoh2*ryoh2 + rzoh2*rzoh2 )

         ! streching

         cosoh1 = (/ rxoh1/roh1, ryoh1/roh1, rzoh1/roh1 /)
         cosoh2 = (/ rxoh2/roh2, ryoh2/roh2, rzoh2/roh2 /)

         ustrech = ustrech + voh(roh1) + voh(roh2)

         force(1:3,iox) = force(1:3,iox) + dvoh(roh1) * cosoh1 + dvoh(roh2) * cosoh2
         force(1:3,ih1) = force(1:3,ih1) - dvoh(roh1) * cosoh1
         force(1:3,ih2) = force(1:3,ih2) - dvoh(roh2) * cosoh2

         ! bending

         costheta = ( rxoh1*rxoh2 + ryoh1*ryoh2 + rzoh1*rzoh2 ) / roh1 / roh2
         theta = dacos( costheta  )
         dtheta = theta - thetaeq
         dtheta2 = dtheta * dtheta
         ftheta = ktheta * dtheta / dsin(theta)
         
         !goto 100
         ubend = ubend + half * ktheta * dtheta2
         
         force(1:3,iox) = force(1:3,iox) + ftheta * &
              ( ( -cosoh1 + costheta * cosoh2 ) / roh2 + &
                ( -cosoh2 + costheta * cosoh1 ) / roh1 )

         force(1:3,ih1) = force(1:3,ih1) + ftheta * &
              ( cosoh2 - costheta * cosoh1 ) / roh1

         force(1:3,ih2) = force(1:3,ih2) + ftheta * &
              ( cosoh1 - costheta * cosoh2 ) / roh2
         !100 continue
         
         ! inter-molecular energy

         do jmol = imol+1, nmols

            jox = 1 + 3*(jmol-1)
            
            rxij = position(1,jox) - position(1,iox)
            ryij = position(2,jox) - position(2,iox)
            rzij = position(3,jox) - position(3,iox)
            
            rij2 = rxij*rxij + ryij*ryij + rzij*rzij
            rij = dsqrt(rij2)
            rij6 = rij2 * rij2 *rij2
            rij7 = rij6 * rij

            cosx = rxij / rij
            cosy = ryij / rij
            cosz = rzij / rij

            cosij = (/ cosx, cosy, cosz /)

            ! oxygen - oxygen LJ interaction
            
            ulj = ulj + ( alj/rij6 - blj ) / rij6

            dV = (-twelve * alj / rij6 + six * blj ) / rij7

            force(1:3,iox) = force(1:3,iox) + dV * cosij

            force(1:3,jox) = force(1:3,jox) - dV * cosij 

            ! Coulomb interaction between charges k and l of imol and jmol

            do k = 1, 3 ! k-charge in imol

               ik = k + 3*(imol-1)
               
               if ( k.eq.1 ) then ! M-point
                  qk = qm
                  rk = rM(:,imol)
               else               ! hydrogen
                  qk = qh
                  rk = position( :, k+3*(imol-1) )
               end if
               
               do l = 1, 3 ! l-charge in jmol

                  il = l + 3*(jmol-1)
                  
                  if ( l.eq.1 ) then ! M-point
                     ql = qm
                     rl = rM(:,jmol)
                  else               ! hydrogen
                     ql = qh
                     rl = position( :, l+3*(jmol-1) )
                  end if

                  rxkl = rl(1) - rk(1)
                  rykl = rl(2) - rk(2)
                  rzkl = rl(3) - rk(3)
                  
                  rkl2 = rxkl*rxkl + rykl*rykl + rzkl*rzkl
                  rkl = dsqrt(rkl2)

                  cosx = rxkl / rkl
                  cosy = rykl / rkl
                  cosz = rzkl / rkl

                  coskl = (/ cosx, cosy, cosz /)
                  
                  ucoulomb = ucoulomb + qk * ql / rkl

                  ! Coulomb forces
                  
                  dV = - qk * ql / rkl2

                  dVcos = dV * coskl

                  if ( k.eq.1 ) then

                     ! forces on Ok, H1k, H2k

                     force(1:3,ik) = force(1:3,ik) + gammao * dVcos

                     force(1:3,ik+1) = force(1:3,ik+1) + gammah * dVcos

                     force(1:3,ik+2) = force(1:3,ik+2) + gammah * dVcos

                  else

                     ! forces on Hk

                     force(1:3,ik) = force(1:3,ik) + dVcos
                     
                  end if

                  if ( l.eq.1 ) then

                     ! forces on Ol, H1l, H2l

                     force(1:3,il) = force(1:3,il) - gammao * dVcos
                                          
                     force(1:3,il+1) = force(1:3,il+1) - gammah * dVcos

                     force(1:3,il+2) = force(1:3,il+2) - gammah * dVcos

                  else

                     ! forces on Hl

                     force(1:3,il) = force(1:3,il) - dVcos
                     
                  end if

               end do ! end l

            end do ! end k

         end do ! end jmol
         
      end do ! end imol

      uinter = ulj + ucoulomb
      uintra = ustrech + ubend
      energy = uinter + uintra

      deallocate(rM)

end subroutine evaluate_energy_qTIP4Pf

!*****************************************************************************
double precision function voh( roh )
!*****************************************************************************
!
!  Project: MolecularDynamics
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
! O-H stretching potential
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  variables

      double precision, intent (IN) :: roh

!*****************************************************************************
!  local variables

      double precision :: deq, deq2, deq3, deq4
      
!*****************************************************************************
!  start function

      deq = roh - req
      deq2 = deq * deq
      deq3 = deq2 * deq
      deq4 = deq2 * deq2
      
      voh = dr * ( alphar2 * deq2 - alphar3 * deq3 + seven/twelve * alphar4 * deq4 )
      
end function voh

!*****************************************************************************
double precision function dvoh( roh )
!*****************************************************************************
!
!  Project: MolecularDynamics
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
! derivative of the O-H stretching potential
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  variables

      double precision, intent (IN) :: roh

!*****************************************************************************
!  local variables

      double precision :: deq, deq2, deq3
      
!*****************************************************************************
!  start function

      deq = roh - req
      deq2 = deq * deq
      deq3 = deq2 * deq
      
      dvoh = dr * ( two * alphar2 * deq - three * alphar3 * deq2 + seven/three * alphar4 * deq3 )
      
end function dvoh    

end module qTIP4Pf
