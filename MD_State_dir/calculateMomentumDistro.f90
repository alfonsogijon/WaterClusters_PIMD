!*****************************************************************************
Subroutine calculateMomentumDistro( state_chain, temperature_target, &
     pdistro_o, pdistro_h, pdistro_h2o )
  
!*****************************************************************************
!
!  Project: WaterElectron_MD
!
!  Module: State (PUBLIC routine)
!
!  Created on: March 2020 by Alfonso Gijón
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This subroutine calculates the momentum distribution of the atoms O and H,
! and the intire molecule of H2O.
!
!*****************************************************************************
!  modules used

      use omp_lib  

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  variables

      type ( StateChain_type ), intent (INOUT) :: state_chain

      real(dp), intent (IN) :: temperature_target

      real(dp), optional, dimension(:), intent (OUT) :: pdistro_o

      real(dp), optional, dimension(:), intent (OUT) :: pdistro_h

      real(dp), optional, dimension(:), intent (OUT) :: pdistro_h2o


!*****************************************************************************
!  local variables

      integer i, j, k
      integer iat, imol

      integer nReplicas
      integer nbinp
      
      integer nAtoms, n_mols
      
      real(dp) :: sigma_o
      real(dp) :: sigma_h
      real(dp) :: sigma_h2o
      real(dp) :: mh2o

      real(dp) :: pmax_o, pmax_h, pmax_h2o, pmin_o, pmin_h, pmin_h2o
      real(dp) :: deltap_o, deltap_h, deltap_h2o
      real(dp) :: px, px_h2o
      
      real (dp), allocatable, dimension (:) :: mass
      double precision, allocatable, dimension (:,:) :: momentum      

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering calculateMomentumDistribution()'
      
      call getStateChain( state_chain, nReplicas = nReplicas, nAtoms = nAtoms )

      n_mols = nAtoms/3
      
      allocate( mass( nAtoms ) )
      allocate( momentum( 3, nAtoms ) )      

      call getStateChain( state_chain, mass = mass )
      
! initialise some quantities

      nbinp = size( pdistro_o(:) )
      pdistro_o = zero
      pdistro_h = zero
      pdistro_h2o = zero
      !sigma_o = dsqrt( boltzmann_k * temperature_target * mass(1) * massConversionFactor )      
      !sigma_h = dsqrt( boltzmann_k * temperature_target * mass(2) * massConversionFactor )      
      !mh2o = ( mass(1) + two*mass(2) ) * massConversionFactor
      sigma_o = dsqrt( boltzmann_k * temperature_target * mass(1) )
      sigma_h = dsqrt( boltzmann_k * temperature_target * mass(2) )
      mh2o = ( mass(1) + two*mass(2) ) 
      sigma_h2o = dsqrt( boltzmann_k * temperature_target * mh2o )      
      pmax_o = three * sigma_o
      pmax_h = three * sigma_h
      pmax_h2o = three * sigma_h2o
      pmin_o = -pmax_o
      pmin_h = -pmax_h
      pmin_h2o = - pmax_h2o
      deltap_o = (pmax_o-pmin_o)/dfloat(nbinp)
      deltap_h = (pmax_h-pmin_h)/dfloat(nbinp)
      deltap_h2o = (pmax_h2o-pmin_h2o)/dfloat(nbinp)

! compute px distribution of O and H atoms and H2O molecule

      ! These two groups of OMP directives are equivalent
      ! (!!!) is faster than (!!), but still slower than the calculation in serial
      
      !!$OMP PARALLEL &
      !!$OMP PRIVATE( i, momentum, imol, px_h2o, k, iat, px, j )
      !!$OMP DO

      !!!$OMP PARALLEL DO &
      !!!$OMP PRIVATE( i, momentum, imol, px_h2o, k, iat, px, j ) &
      !!!$OMP REDUCTION( + : pdistro_o, pdistro_h, pdistro_h2o )
      do i = 1, nReplicas

         call getStateChain( state_chain, iReplica = i, momentum = momentum )

         do imol = 1, n_mols
                 
            px_h2o = zero
                 
            do k = 1, 3
                     
               iat = k + 3*(imol-1)
               px = momentum(1,iat)
               px_h2o = px_h2o + px

               ! atomic momenta
               
               if ( k.eq.1 ) then ! oxygen atom

                  j = floor( ( px - pmin_o ) / deltap_o )
                  if ( 0.le.j .and. j.le.(nbinp-1) ) pdistro_o(j+1) = pdistro_o(j+1) + one
                  !!if ( 0.le.j .and. j.le.(nbinp-1) ) then
                     !!$OMP CRITICAL
                     !!pdistro_o(j+1) = pdistro_o(j+1) + one
                     !!$OMP END CRITICAL
                  !!end if
                  
               else               ! hydrogen atoms 

                  j = floor( ( px - pmin_h ) / deltap_h )
                  if ( 0.le.j .and. j.le.(nbinp-1) ) pdistro_h(j+1) = pdistro_h(j+1) + one
                  !!if ( 0.le.j .and. j.le.(nbinp-1) ) then
                     !!$OMP CRITICAL
                     !!pdistro_h(j+1) = pdistro_h(j+1) + one
                     !!$OMP END CRITICAL
                  !!end if
                  
               end if

            end do

            ! water molecule momentum

            j = floor( ( px_h2o - pmin_h2o ) / deltap_h2o )
            if ( 0.le.j .and. j.le.(nbinp-1) ) pdistro_h2o(j+1) = pdistro_h2o(j+1) + one
            !!if ( 0.le.j .and. j.le.(nbinp-1) ) then
               !!$OMP CRITICAL
               !!pdistro_h2o(j+1) = pdistro_h2o(j+1) + one
               !!$OMP END CRITICAL
            !!end if
                  
         end do ! end of imol
               
      end do ! end of replicas
      !!!$OMP END PARALLEL DO
      
      
      !!$OMP END DO
      !!$OMP END PARALLEL
      
      deallocate( mass )
      deallocate( momentum )

      if ( debug ) write(*,*) 'Exiting calculateMomentumDistribution()'
      
end subroutine calculateMomentumDistro
