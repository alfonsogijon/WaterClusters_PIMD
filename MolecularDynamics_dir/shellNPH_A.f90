!*****************************************************************************
subroutine shellNPH_A( state )
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: MolecularDynamics (PRIVATE routine)
!
!  Created on: Thu 21 Mar 9:35:05 2019  by Eduardo R. Hernandez
!
!  Last Modified: Thu  4 Jul 17:24:04 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! Here we perform NPT dynamics in a different way, without using a barostat, 
! but rather by taking a number of atoms (typically 3), designated as "shell"
! atoms, from whose coordinates we define the periodic cell The
! shell-particle method uses a number of so-called shell particles, that define
! the volume and shape of the simulation box; during the dynamics the simulation
! box is updated from the positions of the shell particles, that are otherwise
! normal particles of the simulation. This algorithm does not rely on a barostat
! to impose the external pressure, i.e. the volume and shape are not dynamical
! variables with associated momenta and masses; rather the cell volume and shape
! result from the dynamics of the shell particles, that respond to the imposed
! external pressure. For details of the general idea, see Uline and Corti, 
! J. Chem. Phys. vol. 123, 164101 (2005), and for the practical details of my
! fully-flexible-cell see Hernandez, J. Chem. Phys. (in preparation). The 
! isothermal part of the dynamics is achieved by coupling the system to a 
! recursive chain of Nosé-Poincaré thermostats. 
!
!*****************************************************************************
!  modules used

      use LU

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  variables

      type ( state_type ), intent (INOUT) :: state

!*****************************************************************************
!  local variables

      logical :: updateReciprocal

      integer :: i

      integer, parameter :: nIterationsMax = 20

      real (dp) :: factor, factor_1, factor_2
      real (dp) :: kinetic
      real (dp) :: shellKinetic
      real (dp) :: stress_energy

      real (dp), parameter :: tolerance = 1.0e-7_dp

      real (dp), allocatable, dimension (:,:) :: lattice_velocity_1
      real (dp), allocatable, dimension (:,:) :: lattice_velocity_2
      real (dp), dimension (3) :: vector
      real (dp), dimension (9) :: flatStress, flatMetric

      type ( representation_type ) :: representation

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering shellNPH_A'

      representation % position = .false.     ! in this routine we need a lattice rep.
      representation % force = .false.
      representation % momentum = .false.

!  the following array will be useful for updating ionic positions

      allocate( lattice_velocity_1( 3, nAtoms ) )
      allocate( lattice_velocity_2( 3, nAtoms ) )

!  make sure coordinate, momenta and forces are in lattice representation

      call get( state, nThermostats = nThermostats,                              &
                position = position, force = force, momentum = momentum,         &
                barostat = barostat, energy = energy,   &
                H0 = H0, Cell = cell, reciprocalCell = reciprocalCell,           &
                stress = stress, representation = representation )

!  we assume that the shell atoms are the last three in the list of atoms; these
!  will be treated slightly differently...

!  and calculate the kinetic energy and kinetic stress

      kinetic = zero
      kinetic_stress = zero

      do i=1, nAtoms - 3

         momentum(:,i) = momentum(:,i) + half_delta_t * force(:,i)

         vector = matmul( reciprocalCell % metricTensor, momentum(1:3,i) ) / mass(i)

         lattice_velocity_1(:,i) = vector

         kinetic = kinetic +  dot_product( vector, momentum(1:3,i) )

         kinetic_stress(1,1) = kinetic_stress(1,1) + mass(i) * vector(1) * vector(1)
         kinetic_stress(1,2) = kinetic_stress(1,2) + mass(i) * vector(1) * vector(2)
         kinetic_stress(1,3) = kinetic_stress(1,3) + mass(i) * vector(1) * vector(3)

         kinetic_stress(2,1) = kinetic_stress(2,1) + mass(i) * vector(2) * vector(1)
         kinetic_stress(2,2) = kinetic_stress(2,2) + mass(i) * vector(2) * vector(2)
         kinetic_stress(2,3) = kinetic_stress(2,3) + mass(i) * vector(2) * vector(3)

         kinetic_stress(3,1) = kinetic_stress(3,1) + mass(i) * vector(3) * vector(1)
         kinetic_stress(3,2) = kinetic_stress(3,2) + mass(i) * vector(3) * vector(2)
         kinetic_stress(3,3) = kinetic_stress(3,3) + mass(i) * vector(3) * vector(3)

      end do

      factor = half 

      energy % kinetic = factor * kinetic

      kinetic_stress = factor * kinetic_stress

!  now we need to advance the Cartesian momenta of the shell atoms to half step

      call advanceShellAtomMomenta

!  after updating the metric tensor components, their contribution to the kinetic energy 
!  has to be accumulated onto that of the standard atoms

      energy % kinetic = energy % kinetic + shellKinetic

!  now we can advance the Cartesian positions of the shell atoms, and the 
!  metric tensor components that depend on them

      call advanceShellAtomPositions

!  update the cell according to the new metric tensor components
!  we also need the new reciprocal cell

!     updateReciprocal = .true.

!     call updateCell( cell, cell % metricTensor, updateReciprocal )

!  advance the lattice components of the atomic positions
!  for this we also need the new reciprocal metric tensor (see equations in paper)

      do i=1, nAtoms - 3

         vector = matmul( reciprocalCell % metricTensor, momentum(1:3,i) ) / mass(i)

         lattice_velocity_2(:,i) = vector

      end do

      do i=1, nAtoms - 3

         position(:,i) = position(:,i) + half_delta_t * ( lattice_velocity_1(:,i) + &
                                         lattice_velocity_2(:,i) )

      end do

      flatStress = pack( simulation % stress % lattice, .true. )
      flatMetric = pack( cell % metricTensor, .true. )

      stress_energy = half * dot_product( flatStress, flatMetric )

      energy % barostat = simulation % pressure * cell % volume + &
                stress_energy

!  before returning, store updated variables in state

      call set( state, position = position, momentum = momentum,       &
                barostat = barostat, energy = energy, Cell = cell,     &
                representation = representation )

!  before leaving, deallocate lattice_velocity storage

      deallocate( lattice_velocity_1 )
      deallocate( lattice_velocity_2 )

      if ( debug ) write(*,*) 'Exiting shellNPH_A'

contains     ! the following routines are specific to NPSM_A, so we keep them here rather than as module routines

!*****************************************************************************
subroutine advanceShellAtomPositions
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Module: MolecularDynamics (PRIVATE routine)
!
!  Created on: Tue  9 Apr 14:50:18 2019  by Eduardo R. Hernandez
!
!  Last Modified: Tue  9 Apr 15:20:46 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
!     Purpose:          This subroutine advances the shell-atom positions; 
!                  once new shell-atom positions are available, the routine
!                  calculates the strain tensor that results from the displacement
!                  from their previous position to the current one, and this is
!                  then used to update the cell lattice vectors and metric
!                  tensor.               
!
!*****************************************************************************
!  No implicits please!
   
     implicit none

!*****************************************************************************
!  shared variables

!*****************************************************************************
!  local variables

      integer :: i, j
      integer :: d, rc

      integer, dimension (9) :: indx

      real (dp) :: factor

      real (dp), dimension (3,3) :: CartesianPosition
      real (dp), dimension (3,3) :: CartesianPositionOld
      real (dp), dimension (3,3) :: CartesianVelocity
      real (dp), dimension (6,6) :: CoeffMatrix
      real (dp), dimension (3,3) :: metricTensor
      real (dp), dimension (3,3) :: H, Ht
      real (dp), dimension (3,3) :: matrix
      real (dp), dimension (6) :: rhs
      real (dp), dimension (3) :: p, r, s, t
      real (dp), dimension (3) :: d21, d31, d32
      real (dp), dimension (3) :: s1, s2, s3

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering advanceShellAtomPositions()'

!  transform the shell-atom positions and momenta to Cartesian representation

      r = (/ cell % a(1), cell % b(1), cell % c(1) /)
      s = (/ cell % a(2), cell % b(2), cell % c(2) /)
      t = (/ cell % a(3), cell % b(3), cell % c(3) /)

      j = 1

      do i=nAtoms-2, nAtoms

         CartesianPositionOld(1,j) = dot_product( r, position(:,i) )
         CartesianPositionOld(2,j) = dot_product( s, position(:,i) )
         CartesianPositionOld(3,j) = dot_product( t, position(:,i) )

         p = matmul( reciprocalCell % MetricTensor, momentum(:,i) )

         CartesianVelocity(1,j) = dot_product( r, p ) / mass(i)
         CartesianVelocity(2,j) = dot_product( s, p ) / mass(i)
         CartesianVelocity(3,j) = dot_product( t, p ) / mass(i)

         j = j + 1

      end do

!  now advance the shell-atom positions

      factor = delta_t

      CartesianPosition = CartesianPositionOld + factor * CartesianVelocity

! TEST TEST TEST convert Cartesian positions of shell atoms to lattice

!     j = 1

!     do i=nAtoms-2, nAtoms
!        
!        position(1,i) = dot_product( reciprocalCell % a, CartesianPosition(:,j) )
!        position(2,i) = dot_product( reciprocalCell % b, CartesianPosition(:,j) )
!        position(3,i) = dot_product( reciprocalCell % c, CartesianPosition(:,j) )

!        j = j + 1

!     end do

! now, with the new Cartesian positions we can obtain the metric tensor 

! the index correspondence is as follows:
! (1,1)->1, (2,2)->2, (3,3)->3, (1,2)->4, (1,3)->5, 
! (2,3)->6, (2,1)->7, (3,1)->8, (3,2)->9

      CoeffMatrix = zero
      rhs = zero

      p = CartesianPosition(:,2) - CartesianPosition(:,1)
      rhs(1) = dot_product( p, p )

      p = CartesianPosition(:,3) - CartesianPosition(:,1)
      rhs(2) = dot_product( p, p )

      p = CartesianPosition(:,3) - CartesianPosition(:,2)
      rhs(3) = dot_product( p, p )

      rhs(4) = dot_product( CartesianPosition(:,1), CartesianPosition(:,1) )
      rhs(5) = dot_product( CartesianPosition(:,2), CartesianPosition(:,2) )
      rhs(6) = dot_product( CartesianPosition(:,3), CartesianPosition(:,3) )

      s1 = position(:,nAtoms-2)
      s2 = position(:,nAtoms-1)
      s3 = position(:,nAtoms)

      d21 = s2 - s1
      d31 = s3 - s1
      d32 = s3 - s2

      CoeffMatrix(1,1) =       d21(1) * d21(1)
      CoeffMatrix(1,2) = two * d21(1) * d21(2)
      CoeffMatrix(1,3) = two * d21(1) * d21(3)
      CoeffMatrix(1,4) =       d21(2) * d21(2)
      CoeffMatrix(1,5) = two * d21(2) * d21(3)
      CoeffMatrix(1,6) =       d21(3) * d21(3)
   
      CoeffMatrix(2,1) =       d31(1) * d31(1)
      CoeffMatrix(2,2) = two * d31(1) * d31(2)
      CoeffMatrix(2,3) = two * d31(1) * d31(3)
      CoeffMatrix(2,4) =       d31(2) * d31(2)
      CoeffMatrix(2,5) = two * d31(2) * d31(3)
      CoeffMatrix(2,6) =       d31(3) * d31(3)

      CoeffMatrix(3,1) =       d32(1) * d32(1)
      CoeffMatrix(3,2) = two * d32(1) * d32(2)
      CoeffMatrix(3,3) = two * d32(1) * d32(3)
      CoeffMatrix(3,4) =       d32(2) * d32(2)
      CoeffMatrix(3,5) = two * d32(2) * d32(3)
      CoeffMatrix(3,6) =       d32(3) * d32(3)

      CoeffMatrix(4,1) =       s1(1) * s1(1)
      CoeffMatrix(4,2) = two * s1(1) * s1(2)
      CoeffMatrix(4,3) = two * s1(1) * s1(3)
      CoeffMatrix(4,4) =       s1(2) * s1(2)
      CoeffMatrix(4,5) = two * s1(2) * s1(3)
      CoeffMatrix(4,6) =       s1(3) * s1(3)

      CoeffMatrix(5,1) =       s2(1) * s2(1)
      CoeffMatrix(5,2) = two * s2(1) * s2(2)
      CoeffMatrix(5,3) = two * s2(1) * s2(3)
      CoeffMatrix(5,4) =       s2(2) * s2(2)
      CoeffMatrix(5,5) = two * s2(2) * s2(3)
      CoeffMatrix(5,6) =       s2(3) * s2(3)

      CoeffMatrix(6,1) =       s3(1) * s3(1)
      CoeffMatrix(6,2) = two * s3(1) * s3(2)
      CoeffMatrix(6,3) = two * s3(1) * s3(3)
      CoeffMatrix(6,4) =       s3(2) * s3(2)
      CoeffMatrix(6,5) = two * s3(2) * s3(3)
      CoeffMatrix(6,6) =       s3(3) * s3(3)

! now we can solve this linear system and obtain the matrix of lattice vectors
! we do this by LU-decomposition of CoeffMatrix

      call LUDCMP( CoeffMatrix, 6, indx, d, rc )

      if ( rc .eq. 0 ) then ! if everything went well

         call LUBKSB( CoeffMatrix, 6, indx, rhs )

      end if

!     write(*,*) 'rc = ', rc

      metricTensor(1,1) = rhs(1)
      metricTensor(1,2) = rhs(2)
      metricTensor(2,1) = rhs(2)
      metricTensor(1,3) = rhs(3)
      metricTensor(3,1) = rhs(3)
      metricTensor(2,2) = rhs(4)
      metricTensor(2,3) = rhs(5)
      metricTensor(3,2) = rhs(5)
      metricTensor(3,3) = rhs(6)

      call updateCell( cell, metricTensor, .true. )

! now we calculate the derivatives of G, the metric tensor, with respect to the positions of shell atoms

! first with respect to x of first shell atom

      rhs = zero

      rhs(1) = -two * ( CartesianPosition(1,2) - CartesianPosition(1,1) )
      rhs(2) = -two * ( CartesianPosition(1,3) - CartesianPosition(1,1) )
      rhs(4) =  two * CartesianPosition(1,1)
 
      call  LUBKSB( CoeffMatrix, 6, indx, rhs )

      matrix(1,1) = rhs(1)
      matrix(1,2) = rhs(2)
      matrix(2,1) = rhs(2)
      matrix(1,3) = rhs(3)
      matrix(3,1) = rhs(3)
      matrix(2,2) = rhs(4)
      matrix(2,3) = rhs(5)
      matrix(3,2) = rhs(5)
      matrix(3,3) = rhs(6)

      dGdr(:,1,1) = pack( matrix, .true. )

! then with respect to y of first shell atom

      rhs = zero

      rhs(1) = -two * ( CartesianPosition(2,2) - CartesianPosition(2,1) )
      rhs(2) = -two * ( CartesianPosition(2,3) - CartesianPosition(2,1) )
      rhs(4) =  two * CartesianPosition(2,1)

      call  LUBKSB( CoeffMatrix, 6, indx, rhs )

      matrix(1,1) = rhs(1)
      matrix(1,2) = rhs(2)
      matrix(2,1) = rhs(2)
      matrix(1,3) = rhs(3)
      matrix(3,1) = rhs(3)
      matrix(2,2) = rhs(4)
      matrix(2,3) = rhs(5)
      matrix(3,2) = rhs(5)
      matrix(3,3) = rhs(6)

      dGdr(:,2,1) = pack( matrix, .true. )

! then with respect to z of first shell atom

      rhs = zero

      rhs(1) = -two * ( CartesianPosition(3,2) - CartesianPosition(3,1) )
      rhs(2) = -two * ( CartesianPosition(3,3) - CartesianPosition(3,1) )
      rhs(4) =  two * CartesianPosition(3,1)
 
      call  LUBKSB( CoeffMatrix, 6, indx, rhs )

      matrix(1,1) = rhs(1)
      matrix(1,2) = rhs(2)
      matrix(2,1) = rhs(2)
      matrix(1,3) = rhs(3)
      matrix(3,1) = rhs(3)
      matrix(2,2) = rhs(4)
      matrix(2,3) = rhs(5)
      matrix(3,2) = rhs(5)
      matrix(3,3) = rhs(6)

      dGdr(:,3,1) = pack( matrix, .true. )

! then with respect to x of second shell atom

      rhs = zero

      rhs(1) =  two * ( CartesianPosition(1,2) - CartesianPosition(1,1) )
      rhs(3) = -two * ( CartesianPosition(1,3) - CartesianPosition(1,2) )
      rhs(5) =  two * CartesianPosition(1,2)

      call  LUBKSB( CoeffMatrix, 6, indx, rhs )

      matrix(1,1) = rhs(1)
      matrix(1,2) = rhs(2)
      matrix(2,1) = rhs(2)
      matrix(1,3) = rhs(3)
      matrix(3,1) = rhs(3)
      matrix(2,2) = rhs(4)
      matrix(2,3) = rhs(5)
      matrix(3,2) = rhs(5)
      matrix(3,3) = rhs(6)

      dGdr(:,1,2) = pack( matrix, .true. )

! then with respect to y of second shell atom

      rhs = zero

      rhs(1) =  two * ( CartesianPosition(2,2) - CartesianPosition(2,1) )
      rhs(3) = -two * ( CartesianPosition(2,3) - CartesianPosition(2,2) )
      rhs(5) =  two * CartesianPosition(2,2)
 
      call  LUBKSB( CoeffMatrix, 6, indx, rhs )

      matrix(1,1) = rhs(1)
      matrix(1,2) = rhs(2)
      matrix(2,1) = rhs(2)
      matrix(1,3) = rhs(3)
      matrix(3,1) = rhs(3)
      matrix(2,2) = rhs(4)
      matrix(2,3) = rhs(5)
      matrix(3,2) = rhs(5)
      matrix(3,3) = rhs(6)

      dGdr(:,2,2) = pack( matrix, .true. )

! then with respect to z of second shell atom

      rhs = zero

      rhs(1) =  two * ( CartesianPosition(3,2) - CartesianPosition(3,1) )
      rhs(3) = -two * ( CartesianPosition(3,3) - CartesianPosition(3,2) )
      rhs(5) =  two * CartesianPosition(3,2)
 
      call  LUBKSB( CoeffMatrix, 6, indx, rhs )

      matrix(1,1) = rhs(1)
      matrix(1,2) = rhs(2)
      matrix(2,1) = rhs(2)
      matrix(1,3) = rhs(3)
      matrix(3,1) = rhs(3)
      matrix(2,2) = rhs(4)
      matrix(2,3) = rhs(5)
      matrix(3,2) = rhs(5)
      matrix(3,3) = rhs(6)

      dGdr(:,3,2) = pack( matrix, .true. )

! then with respect to x of third shell atom

      rhs = zero

      rhs(2) =  two * ( CartesianPosition(1,3) - CartesianPosition(1,1) )
      rhs(3) =  two * ( CartesianPosition(1,3) - CartesianPosition(1,2) )
      rhs(6) =  two * CartesianPosition(1,3)

      call  LUBKSB( CoeffMatrix, 6, indx, rhs )

      matrix(1,1) = rhs(1)
      matrix(1,2) = rhs(2)
      matrix(2,1) = rhs(2)
      matrix(1,3) = rhs(3)
      matrix(3,1) = rhs(3)
      matrix(2,2) = rhs(4)
      matrix(2,3) = rhs(5)
      matrix(3,2) = rhs(5)
      matrix(3,3) = rhs(6)

      dGdr(:,1,3) = pack( matrix, .true. )

! then with respect to y of third shell atom

      rhs = zero

      rhs(2) =  two * ( CartesianPosition(2,3) - CartesianPosition(2,1) )
      rhs(3) =  two * ( CartesianPosition(2,3) - CartesianPosition(2,2) )
      rhs(6) =  two * CartesianPosition(2,3)
 
      call  LUBKSB( CoeffMatrix, 6, indx, rhs )

      matrix(1,1) = rhs(1)
      matrix(1,2) = rhs(2)
      matrix(2,1) = rhs(2)
      matrix(1,3) = rhs(3)
      matrix(3,1) = rhs(3)
      matrix(2,2) = rhs(4)
      matrix(2,3) = rhs(5)
      matrix(3,2) = rhs(5)
      matrix(3,3) = rhs(6)

      dGdr(:,2,3) = pack( matrix, .true. )

! then finally with respect to z of third shell atom

      rhs = zero

      rhs(2) =  two * ( CartesianPosition(3,3) - CartesianPosition(3,1) )
      rhs(3) =  two * ( CartesianPosition(3,3) - CartesianPosition(3,2) )
      rhs(6) =  two * CartesianPosition(3,3)
 
      call  LUBKSB( CoeffMatrix, 6, indx, rhs )

      matrix(1,1) = rhs(1)
      matrix(1,2) = rhs(2)
      matrix(2,1) = rhs(2)
      matrix(1,3) = rhs(3)
      matrix(3,1) = rhs(3)
      matrix(2,2) = rhs(4)
      matrix(2,3) = rhs(5)
      matrix(3,2) = rhs(5)
      matrix(3,3) = rhs(6)

      dGdr(:,3,3) = pack( matrix, .true. )

! TEST TEST TEST convert Cartesian positions of shell atoms to lattice

!     j = 1

!     do i=nAtoms-2, nAtoms
!        
!        position(1,i) = dot_product( reciprocalCell % a, CartesianPosition(:,j) )
!        position(2,i) = dot_product( reciprocalCell % b, CartesianPosition(:,j) )
!        position(3,i) = dot_product( reciprocalCell % c, CartesianPosition(:,j) )

!        j = j + 1

!     end do

      if ( debug ) write(*,*) 'Exiting advanceShellAtomPositions()'

end subroutine advanceShellAtomPositions

!*****************************************************************************
subroutine advanceShellAtomMomenta
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Created on: Tue  9 Apr 09:43:03 2019  by Eduardo R. Hernandez
!
!  Last Modified: Tue  9 Apr 14:31:50 2019
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
!                       This subroutine updates the momenta of the metric 
!                       tensor components to half-step. The equations that 
!                       appear in the Souza-Martins algorithm for constant
!                       pressure molecular dynamics can be integrated using
!                       the Generalised Leap-Frog (GLF) integration method.
!                       When this is done, the equations that appear for 
!                       the numerical integration of the momenta of the 
!                       metric tensor components are implicit, and therefore
!                       we must use a Newton-Raphson procedure to solve them.
!                       This is what this subroutine does.
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  shared variables

!*****************************************************************************
!  local variables

      integer :: i, j

      real (dp), dimension (3,3) :: CartesianForce
      real (dp), dimension (3,3) :: CartesianMomentum
      real (dp), dimension (3,3) :: cellForce
      real (dp), dimension (9) :: difference
      real (dp), dimension (9) :: flatCellForce
      real (dp), dimension (3) :: p, r, s, t, vector

!*****************************************************************************
!  start of subroutine

      if ( debug ) write(*,*) 'Entering advanceShellAtomMomenta()'

! convert the shell-atom momenta to Cartesian; likewise with the forces

      ! this transforms m_i G_ab \dot{s_ib} first to m_i \dot{s_ia}
      ! and then to m_i H_ab \dot_{s_ib} = m_i v_ia

      r = (/ cell % a(1), cell % b(1), cell % c(1) /)
      s = (/ cell % a(2), cell % b(2), cell % c(2) /)
      t = (/ cell % a(3), cell % b(3), cell % c(3) /)

      j = 1

      do i=nAtoms-2, nAtoms

         p = matmul( reciprocalCell % metricTensor, momentum(1:3,i) )

         CartesianMomentum(1,j) = dot_product( r, p )
         CartesianMomentum(2,j) = dot_product( s, p )
         CartesianMomentum(3,j) = dot_product( t, p )

         j = j + 1

      end do 

! now convert the forces

      j = 1

      do i=nAtoms-2, nAtoms

         CartesianForce(1,j) = dot_product( reciprocalCell % a, force(1:3,i) )
         CartesianForce(2,j) = dot_product( reciprocalCell % b, force(1:3,i) )
         CartesianForce(3,j) = dot_product( reciprocalCell % c, force(1:3,i) )

         j = j + 1

      end do

! now calculate the contributions to the shell-atom forces coming from the cell

      cellForce = -( stress % lattice + kinetic_stress +                    &
                       half * simulation % pressure * cell % volume *       &
                       reciprocalCell % metricTensor +                      &
                       half * simulation % stress % lattice )

! with this we have everything we need to advance the shell-atom momenta to half step

      flatCellForce = pack( cellForce, .true. )

      do i=1, 3

         CartesianMomentum(1,i) = CartesianMomentum(1,i) +                  &
                half_delta_t * ( CartesianForce(1,i) +  &
                      dot_product( flatCellForce, dGdr(:,1,i) ) )

         CartesianMomentum(2,i) = CartesianMomentum(2,i) +                  &
                half_delta_t * ( CartesianForce(2,i) +  &
                      dot_product( flatCellForce, dGdr(:,2,i) ) )

         CartesianMomentum(3,i) = CartesianMomentum(3,i) +                  &
                half_delta_t * ( CartesianForce(3,i) +  &
                      dot_product( flatCellForce, dGdr(:,3,i) ) )

      end do

!     write(*,'(3f17.7)') CartesianMomentum(1:3,1)
!     write(*,'(3f17.7)') CartesianMomentum(1:3,2)
!     write(*,'(3f17.7)') CartesianMomentum(1:3,3)

! finally, transform the shell-atom momenta back to lattice reppresentation, 
! and in passing, calculate the shell-atom contribution to the kinetic energy

      shellKinetic = zero

      j = 1

      do i=nAtoms-2, nAtoms

         shellKinetic = shellKinetic +                                      &
           dot_product( CartesianMomentum(:,j), CartesianMomentum(:,j) ) /  &
                          mass(i)

         p(1) = dot_product( reciprocalCell % a,         &
                                       CartesianMomentum(1:3,j) )
         p(2) = dot_product( reciprocalCell % b,         &
                                       CartesianMomentum(1:3,j) )
         p(3) = dot_product( reciprocalCell % c,         &
                                       CartesianMomentum(1:3,j) )

         momentum(:,i) = matmul( cell % metricTensor, p )

         j = j + 1

      end do

      shellKinetic = half * shellKinetic 

      if ( debug ) write(*,*) 'Exiting advanceShellAtomMomenta()'

end subroutine advanceShellAtomMomenta

end subroutine shellNPH_A

