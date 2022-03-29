!*****************************************************************************
program PIMD_Water
!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Created on: Jan 2020 by Alfonso Gijón
!
! This means to serve as an example of a client program that runs a
! multireplica simulation using the Path Integral interpretation of quantum
! statistical mechanics. A Langevin dynamics is followed by each atom of
! of the extended system, to ensure that the NVT ensemble is correctly sampled
!
! The physical system is a cluster of water molecules, interacting via a SPC/F
! potential. 
!
!*****************************************************************************
!  modules used

      use MD_State_module
      use MolecularDynamics
      use SPCF

      use omp_lib

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  local variables
      
      integer i, j, k
      integer i1, i2, ir
      integer n, n0
      integer nStep
      integer iat, imol, jmol, nmol
      integer iblock

      integer length
      integer energ
      integer dyntype
      integer nSteps
      integer nSteps_equi
      integer nSteps_prod            
      integer n_atoms
      integer n_mols
      integer nReplicas
      integer nblocks
      
      integer jat, iox, jox, n_ox
      
      integer istart
      integer noutputs, nperiod, nperiod_density

      integer nthreads, nthreads_max

      integer nThermostats
      double precision, dimension(:), allocatable :: tMass
      
      logical fixCellShape
      logical periodic
      logical start

      double precision, parameter :: ang_to_au = one / 0.529177d0
      double precision, parameter :: kcalmol_to_au = one / 627.510d0
      double precision, parameter :: massConversionFactor = 1822.8885302342505_dp

      double precision :: gamma
      double precision :: temperature_target
      double precision :: delta_t      

      double precision :: tstart_cpu, tfinish_cpu, cputime
      integer :: count_0, count_1, count_rate, count_max      
      double precision :: walltime
      
      double precision :: kinetic_energy
      double precision :: potential_energy
      double precision :: thermostat_energy
      double precision :: total_energy
      double precision :: temperature
      double precision :: beta, beta2, kT
      double precision :: qKineticEnergy
      double precision :: qKineticEnergy_v2
      double precision :: qPotentialEnergy
      double precision :: qTotalEnergy
      double precision :: qTotalEnergy_v2

      double precision :: temperature_avrg, kinetic_avrg, potential_avrg, totalenergy_avrg
      double precision :: qkinetic_avrg, qkinetic_v2_avrg, qpotential_avrg
      double precision :: qtotalenergy_avrg, qtotalenergy_v2_avrg
      double precision :: sigma_kinetic, sigma_potential, sigma_totalenergy, sigma_qtotalenergy_v2
      double precision :: sigma_qkinetic, sigma_qkinetic_v2, sigma_qpotential, sigma_qtotalenergy

      double precision, allocatable, dimension(:) :: gyration2, gyration2_avrg
      double precision, allocatable, dimension(:,:) :: g2_avrg
      double precision :: radius2, radius2_avrg, radmax2, radmax2_avrg
      double precision :: gyration2_o, gyration2_h, gyration2_o_avrg, gyration2_h_avrg
      double precision :: sigma_radius2, sigma_gyration2_o, sigma_gyration2_h, sigma_radmax2
      
      double precision:: msd_ox_avrg, msd_hy_avrg, sigma_msd_ox, sigma_msd_hy
      double precision:: msd_ox, msd_hy

      double precision:: msd_w

      double precision :: roh2, theta
      double precision :: roh2_avrg, theta_avrg
      double precision :: sigma_roh2, sigma_theta

      integer, parameter :: nblocks_max = 20
      double precision, dimension (nblocks_max) :: kin_avrg, pot_avrg, tote_avrg, qkin_avrg
      double precision, dimension (nblocks_max) :: qkin_v2_avrg, qpot_avrg
      double precision, dimension(nblocks_max) :: r2_avrg, qtote_avrg, qtote_v2_avrg
      double precision, dimension(nblocks_max) :: rmax2_avrg
      double precision, dimension(nblocks_max) :: doh2_avrg, hoh_avrg
      double precision :: dblock

      ! integer, parameter :: nbinp = 100
      ! double precision :: pmin_o, pmax_o, pmin_h, pmax_h, deltap_o, deltap_h
      ! double precision :: pmin_h2o, pmax_h2o, deltap_h2o, mh2o, px_h2o
      ! double precision :: norm_ph, norm_po, norm_ph2o
      ! double precision :: pl, pr, pc, px, gaussian, sigma
      ! double precision :: sigma_o, sigma_h, sigma_h2o
      ! double precision, dimension (nbinp) :: pdistro_h_avrg
      ! double precision, dimension (nbinp) :: pdistro_o_avrg
      ! double precision, dimension (nbinp) :: pdistro_h2o_avrg
      ! double precision, dimension (nbinp) :: pdistro_h
      ! double precision, dimension (nbinp) :: pdistro_o
      ! double precision, dimension (nbinp) :: pdistro_h2o
      ! double precision, allocatable, dimension (:,:) :: momentum

      integer, parameter :: nbinx = 200
      double precision :: x1c, x2c, xmax, xmin, dx
      double precision :: dox_xy, dox_yz, dox_zx, dhy_xy, dhy_yz, dhy_zx
      double precision, dimension(3,nbinx,nbinx) :: density_o
      double precision, dimension(3,nbinx,nbinx) :: density_h
      double precision, dimension(3,nbinx,nbinx) :: density_o_avrg      
      double precision, dimension(3,nbinx,nbinx) :: density_h_avrg

      integer:: nbinr
      double precision :: rc, rmax_rho, dr_rho, rmax_g, dr_g, rmax_histo
      double precision, allocatable, dimension(:) :: rhor_o
      double precision, allocatable, dimension(:) :: rhor_h
      double precision, allocatable, dimension(:) :: rhor_h2o      
      double precision, allocatable, dimension(:) :: rhor_o_avrg
      double precision, allocatable, dimension(:) :: rhor_h_avrg
      double precision, allocatable, dimension(:) :: rhor_h2o_avrg
      double precision, allocatable, dimension(:) :: goo, ghh, goh
      double precision, allocatable, dimension(:) :: goo_avrg, ghh_avrg, goh_avrg
      
      double precision, dimension(3) :: a, b, c
      double precision, dimension(3) :: rcm
      double precision, allocatable, dimension (:,:) :: position0
      double precision, allocatable, dimension (:,:) :: centroid
      double precision, allocatable, dimension (:,:) :: centroid0              
      double precision, allocatable, dimension (:,:) :: centroid_avrg
      double precision, allocatable, dimension (:,:,:) :: ctr_avrg
      double precision, allocatable, dimension (:,:) :: ctr2_avrg      
      double precision, allocatable, dimension (:) :: centroid2_avrg
      double precision, allocatable, dimension (:) :: centroid2
      double precision, allocatable, dimension (:,:,:) :: polymer
      double precision, allocatable, dimension (:) :: mass
      
      double precision, allocatable, dimension (:,:) :: rij, rij2
      double precision, allocatable, dimension (:,:,:) :: rij_block, rij2_block
      double precision, dimension(3) :: vij
      double precision :: dij, dij2, lindemann, linde, sigma_linde

      type ( StateChain_type ) :: state_chain
      
      character ( len = 30 ) :: restartFile
      character ( len = 30 ) :: newRestartFile
      character ( len = 30 ) :: logfile, xyzfile
      character ( len = 30 ) :: idat
      logical exist

!*****************************************************************************
! begin program

      call cpu_time(tstart_cpu)
      call system_clock( count_0, count_rate, count_max )
      !walltime = omp_get_wtime()
      
! simulation info
      
      restartFile = 'restartPI'
      newRestartFile = 'restartPI_2'

      write(*,*) 'Reading input'
      write(*,*)

      read(5,*)
      read(5,*) logfile
      logfile = trim(logfile)
      write(*,*) '# logfile'
      write(*,*) logfile
      !logfile = "logfile.dat"
      inquire( file = logfile, exist = exist )
      if (exist) then
         open( unit = 20, file = logfile, status = "old", position = "append", action = "write")
      else
         open( unit = 20, file = logfile)
      end if
      
      read(5,*)
      read(5,*) nSteps_equi, nSteps_prod
      write(*,*) '# Number of timesteps of equilibration and production'
      write(*,*) nSteps_equi, nSteps_prod
      write(20,*) '# Number of timesteps of equilibration and production'
      write(20,*) nSteps_equi, nSteps_prod      
      !nSteps = 100000
      nSteps = nSteps_equi + nSteps_prod

      read(5,*)
      read(5,*) noutputs
      write(*,*) '# Number of outputs in logfile'
      write(*,*) noutputs
      write(20,*) '# Number of outputs in logfile'
      write(20,*) noutputs      
      !noutputs = 1000
      nperiod = nSteps_prod / noutputs

      read(5,*)
      read(5,*) nperiod_density, nbinr, rmax_histo
      write(*,*) '# Frequency to compute spatial density, number of bins and rmax (in Angstroms)'
      write(*,*) nperiod_density, nbinr, rmax_histo
      write(20,*) '# Frequency to compute spatial density, number of bins and rmax (in Angstroms)'
      write(20,*) nperiod_density, nbinr, rmax_histo
      !nperiod_density = 1
      !nbinr = 500

      dyntype = 4  ! Bussi-Parrinello dynamics for each degree of freedom

      read(5,*)
      read(5,*) istart
      write(*,*) '# Restart label (1 if restart, 0 if not)'
      write(*,*) istart
      write(20,*) '# Restart label (1 if restart, 0 if not)'
      write(20,*) istart      
      if (istart .eq. 1) then
         start = .false.
      else
         start = .true.
      end if
      !start = .true.    ! this is not a re-start simulation
      !start = .false.   ! restart simulation
      
      periodic = .false. ! we do not assume periodic-boundary conditions
      fixCellShape = .false. ! we do not impose a fixed cell shape
      
      length = 0 ! length units of client are Bohr (atomic units)
      energ = 0 ! energy units are Hartree (atomic units)
      call setUnits( Bohr, Hartree )

      ! LATTICE VECTORS ARE NOT USED BY THE PROGRAM, AS WE ARE
      ! CONSIDERING FINITE WATER CLUSTERS, NOT PERIODIC SYSTEMS. 
      ! a = zero
      ! b = zero
      ! c = zero
      ! read(5,*)
      ! read(5,*) a(1), b(2), c(3)
      ! write(*,*) '# Lattice vectors'
      ! write(*,*) a
      ! write(*,*) b
      ! write(*,*) c
      ! write(20,*) '# Lattice vectors'
      ! write(20,*) a
      ! write(20,*) b
      ! write(20,*) c            
      a = (/60.55_dp, 0.0_dp, 0.0_dp/)
      b = (/0.0_dp, 60.55_dp, 0.0_dp/)
      c = (/0.0_dp, 0.0_dp, 60.55_dp/)

      if ( length .eq. 1 ) then ! cell parameters and atomic positions 
                                ! must be in client units (Angstrom in this case)
         a = a / Angstrom
         b = b / Angstrom
         c = c / Angstrom

      end if

      nThermostats = 0
      !tMass = 1.0_dp

      read(5,*)
      read(5,*) n_mols
      write(*,*) '# Number of water molecules of the cluster'
      write(*,*) n_mols
      write(20,*) '# Number of water molecules of the cluster'
      write(20,*) n_mols      
      !n_mols = 2
      n_atoms = 3 * n_mols
      
      read(5,*)
      read(5,*) nReplicas
      write(*,*) '# Number of replicas'
      write(*,*) nReplicas
      write(20,*) '# Number of replicas'
      write(20,*) nReplicas      
      !nReplicas = 120

      read(5,*)
      read(5,*) gamma
      write(*,*) '# Friction parameter of Langevin dynamics (NVE dynamics if zero)'
      write(*,*) gamma
      write(20,*) '# Friction parameter of Langevin dynamics (NVE dynamics if zero)'
      write(20,*) gamma      
      !gamma = 0.001_dp   ! friction parameter of Langevin dynamics
      !gamma = zero      ! if gamma = zero, NVE dynamics is followed

      read(5,*)
      read(5,*) temperature_target
      write(*,*) '# Target temperature (K)'
      write(*,*) temperature_target
      write(20,*) '# Target temperature (K)'
      write(20,*) temperature_target      
      !temperature_target = 50.0_dp
      beta = one / ( boltzmann_k * temperature_target )
      beta2 = beta*beta
      kT = one / beta

      read(5,*)
      read(5,*) delta_t
      write(*,*) '# Timestep (fs)'
      write(*,*) delta_t
      write(20,*) '# Timestep (fs)'
      write(20,*) delta_t      
      !delta_t = 1.0_dp ! 1 fs

      ! write some information about units, number of atoms...
      
      if ( length .eq. 0 ) then
         write(20,*) 'Length units are Bohr (atomic units)'
      else if ( length .eq. 1 ) then
         write(20,*) 'Length units are Angstrom (puaj!!!)' 
      else
         write(20,*) 'Unrecognised length units!'
         stop
      end if

      if ( energ .eq. 0 ) then
         write(20,*) 'Energy units are Hartree (atomic units)'
      else if ( energ .eq. 1 ) then
         write(20,*) 'Energy units are Rydberg (tutut!)' 
      else if ( energ .eq. 2 ) then
         write(20,*) 'Energy units are eV (puaj!!!!)' 
      else
         write(20,*) 'Unrecognised energy units!'
         stop
      end if

      write(20,*) 'n_atoms = ', n_atoms

! allocate some arrays
      
      allocate( mass( n_atoms ) )
      allocate( position0( 3, n_atoms ) )
      !allocate( momentum( 3, n_atoms ) )
      allocate( centroid( 3, n_atoms ) )
      allocate( centroid0( 3, n_atoms ) )      
      allocate( centroid_avrg( 3, n_atoms ) )
      allocate( centroid2( n_atoms ) )            
      allocate( centroid2_avrg( n_atoms ) )      
      allocate( polymer( 3, n_atoms, nReplicas ) )
      allocate( gyration2( n_atoms ) )
      allocate( gyration2_avrg( n_atoms ) )

      n_ox = n_mols
      allocate( rij( n_ox, n_ox ) )
      allocate( rij2( n_ox, n_ox ) )      

      allocate( rhor_o(nbinr) )
      allocate( rhor_h(nbinr) )
      allocate( rhor_h2o(nbinr) )
      allocate( rhor_o_avrg(nbinr) )
      allocate( rhor_h_avrg(nbinr) )
      allocate( rhor_h2o_avrg(nbinr) )
      allocate( goo(nbinr) )
      allocate( ghh(nbinr) )
      allocate( goh(nbinr) )
      allocate( goo_avrg(nbinr) )
      allocate( ghh_avrg(nbinr) )
      allocate( goh_avrg(nbinr) )

! parallelization settings

      ! Choose the number of threads

      !nthreads = omp_get_num_threads ( )      ! OMP_NUM_THREADS enviroment variable
      nthreads = omp_get_num_procs( )         ! number of processors of my machine

      nthreads_max = omp_get_max_threads ( )
      call omp_set_num_threads(nthreads)

      write(*,*)
      write(*,*) 'Number of OpenMP threads'
      write(*,*) nthreads
      write(*,*)

      write(20,*)
      write(20,*) 'Number of OpenMP threads'
      write(20,*) nthreads
      write(20,*)                        

! state info

      read(5,*)
      read(5,*) xyzfile
      xyzfile = trim(xyzfile)
      write(*,*) 'Structures file'
      write(*,*) xyzfile
      write(20,*) 'Structures file'
      write(20,*) xyzfile      
      close(5)
      
      if ( start ) then

         ! read positions of atoms
      
         open( unit = 10, file = xyzfile )

         read(10,*)
         read(10,*)
         do imol = 1, n_mols
            do j = 1, 3
               iat = j + 3 * (imol-1)
               read(10,*) idat, position0(1:3,iat) 
               if ( j .eq. 1) then
                  mass(iat) = 15.999_dp ! oxygen mass
               else
                  mass(iat) = 1.00784_dp ! hydrogen mass
               end if
            end do
         end do
         
         close(10)

         ! Translate the center of masses of the cluster to the origin

         rcm = zero
         do imol = 1, n_mols
            do j = 1, 3
               iat = j + 3 * (imol-1)
               rcm = rcm + mass(iat) * position0(:,iat)
            end do
         end do
         rcm = rcm / sum( mass )
         do iat = 1, n_atoms
            position0(:,iat) = position0(:,iat) - rcm(:)
         end do
         position0 = position0 * ang_to_au ! convert Angstrom to atomic units

         ! now create a chain of replicas (with same positions in each replica but random momenta)
         
         call setupStateChain( state_chain, nReplicas, temperature_target, n_atoms, &
              nThermostats, tMass, a, b, c, mass, position0, periodic = periodic )
         n0 = 0

      else
         
         call readRestartChain( state_chain, restartFile, n0 )
         call getMass( state_chain, mass )
         
      end if

! setup SPCF water potential and compute initial energy and forces

      call SPCF_setup()

      call computeForces( state_chain, temperature_target, gyration2, centroid )

! compute initial radius to set size of the spatial mesh
      
      call CentroidProperties( n_atoms, mass, centroid, radius2, radmax2, roh2, theta )
      centroid0 = centroid
      open(unit=88,file='msd.dat')
      
      call MSD(n_atoms, mass, centroid0, centroid, msd_w)
      write(88,*) int(0), msd_w / ang_to_au

! before entering the MD loop, set-up MD and compute initial classical energy

      call MDChain_setup( state_chain, dyntype, delta_t, temperature_target, friction = gamma )

      call inquireStateChain(state_chain, temperature = temperature, &
           kineticEnergy = kinetic_energy, potentialEnergy = potential_energy)

      total_energy = kinetic_energy + potential_energy

      write(*,*) 'INITIAL STEP'
      write(*,'(i10,10f17.8)') n0, temperature, kinetic_energy, potential_energy, &
           total_energy, qKineticEnergy, qKineticEnergy_v2, qPotentialEnergy, &
           qTotalEnergy, qTotalEnergy_v2
      write(*,*)
      
      write(20,*) 'INITIAL STEP'
      write(20,'(i10,10f17.8)') n0, temperature, kinetic_energy, potential_energy, &
           total_energy, qKineticEnergy, qKineticEnergy_v2, qPotentialEnergy, &
           qTotalEnergy, qTotalEnergy_v2
      write(20,*)

! initialise some quantities      
      
      ! initialise momentum distributions
      
      ! pdistro_o_avrg = zero
      ! pdistro_h_avrg = zero
      ! pdistro_h2o_avrg = zero
      ! sigma_o = dsqrt( boltzmann_k * temperature_target * mass(1) * massConversionFactor )
      ! sigma_h = dsqrt( boltzmann_k * temperature_target * mass(2) * massConversionFactor )      
      ! mh2o = ( mass(1) + two*mass(2) ) * massConversionFactor
      ! sigma_h2o = dsqrt( boltzmann_k * temperature_target * mh2o )      
      ! pmax_o = three * sigma_o
      ! pmax_h = three * sigma_h
      ! pmax_h2o = three * sigma_h2o
      ! pmin_o = -pmax_o
      ! pmin_h = -pmax_h
      ! pmin_h2o = - pmax_h2o
      ! deltap_o = (pmax_o-pmin_o)/dfloat(nbinp)
      ! deltap_h = (pmax_h-pmin_h)/dfloat(nbinp)
      ! deltap_h2o = (pmax_h2o-pmin_h2o)/dfloat(nbinp)

      ! initialise averages
      
      temperature_avrg = zero
      kinetic_avrg = zero
      potential_avrg = zero
      totalenergy_avrg = zero
      qtotalenergy_avrg = zero
      qtotalenergy_v2_avrg = zero      
      qkinetic_avrg = zero
      qkinetic_v2_avrg = zero      
      qpotential_avrg = zero

      kin_avrg = zero
      pot_avrg = zero
      tote_avrg = zero
      qkin_avrg = zero
      qkin_v2_avrg = zero
      qpot_avrg = zero
      qtote_avrg = zero
      qtote_v2_avrg = zero
      
      if( nSteps_prod .gt. nblocks_max ) then
         nblocks = nblocks_max
      else
         nblocks = nSteps_prod
      end if
      dblock = dfloat(nSteps_prod)/dfloat(nblocks)

      allocate( ctr_avrg( 3, n_atoms, nblocks ) )
      allocate( ctr2_avrg( n_atoms, nblocks ) )
      allocate( g2_avrg( n_atoms, nblocks ) )

      allocate(rij_block(n_ox,n_ox,nblocks) )
      allocate(rij2_block(n_ox,n_ox,nblocks) )      

      radius2_avrg = zero
      radmax2_avrg = zero            
      gyration2_avrg = zero
      r2_avrg = zero
      rmax2_avrg = zero
      centroid_avrg = zero
      centroid2_avrg = zero
      ctr_avrg = zero
      ctr2_avrg = zero
      g2_avrg = zero

      rij_block = zero
      rij2_block = zero      

      roh2_avrg = zero
      theta_avrg = zero
      doh2_avrg = zero
      hoh_avrg = zero

      density_o_avrg = zero      
      density_h_avrg = zero

      rhor_o_avrg = zero
      rhor_h_avrg = zero
      rhor_h2o_avrg = zero

      goo_avrg = zero
      ghh_avrg = zero
      goh_avrg = zero

      rij = zero
      rij2 = zero
      
! launch dynamics

      do n = 1, nSteps

         nStep = n0 + n

         ! advances positions of every replica by a time step, and momenta by half a time step

         call BussiParrinello_A( state_chain )
         !call velocityVerlet_StateChain_A( state_chain )

         ! update forces

         call computeForces( state_chain, temperature_target, gyration2, centroid )

         ! with the new forces, uptade the remaining half a time step of momenta for every replica

         call BussiParrinello_B( state_chain )
         !call velocityVerlet_StateChain_B( state_chain )

         ! Compute centroid properties and momentum distribution

         if ( n.gt. nSteps_equi ) then
         
            if ( mod(n, nperiod_density) .eq. 0 ) then

               call CentroidProperties( n_atoms, mass, centroid, radius2, radmax2, roh2, theta, &
                    rmax_histo = rmax_histo, density_o = density_o, density_h = density_h, &
                    rhor_o = rhor_o, rhor_h = rhor_h, rhor_h2o = rhor_h2o, &
                    goo = goo, ghh = ghh, goh = goh)

               call MSD(n_atoms, mass, centroid0, centroid, msd_w)
               write(88,*) n, msd_w / ang_to_au               

               ! call CentroidProperties( n_atoms, mass, centroid, radius2, roh2, theta, &
               !      rmax_histo = rmax_histo, &
               !      rhor_o = rhor_o, rhor_h = rhor_h, rhor_h2o = rhor_h2o, &
               !      goo = goo, ghh = ghh, goh = goh)
            
               ! call calculateMomentumDistro( state_chain, temperature_target, &
               !      pdistro_o, pdistro_h, pdistro_h2o )

               ! pdistro_o_avrg = pdistro_o_avrg + pdistro_o
               ! pdistro_h_avrg = pdistro_h_avrg + pdistro_h
               ! pdistro_h2o_avrg = pdistro_h2o_avrg + pdistro_h2o         
         
               density_o_avrg = density_o_avrg + density_o
               density_h_avrg = density_h_avrg + density_h
            
               rhor_o_avrg = rhor_o_avrg + rhor_o
               rhor_h_avrg = rhor_h_avrg + rhor_h
               rhor_h2o_avrg = rhor_h2o_avrg + rhor_h2o

               goo_avrg = goo_avrg + goo
               ghh_avrg = ghh_avrg + ghh
               goh_avrg = goh_avrg + goh
            
            else
            
               call CentroidProperties( n_atoms, mass, centroid, radius2, radmax2, roh2, theta )

            end if

         end if ! end if n > nSteps_equi

         ! instantaneous State information
         
         call inquireStateChain(state_chain, temperature = temperature, &
              kineticEnergy = kinetic_energy, potentialEnergy = potential_energy, &
              qKineticEnergy = qKineticEnergy, qPotentialEnergy = qPotentialEnergy, &
              qKineticEnergy_v2 = qKineticEnergy_v2)

         ! accumulate averages

         if( n.gt.nSteps_equi ) then

            gyration2_avrg = gyration2_avrg + gyration2
            radius2_avrg = radius2_avrg + radius2
            radmax2_avrg = radmax2_avrg + radmax2            
            centroid_avrg = centroid_avrg + centroid
            roh2_avrg = roh2_avrg + roh2
            theta_avrg = theta_avrg + theta

            do iat = 1, n_atoms
               centroid2(iat) = dot_product( centroid(:,iat), centroid(:,iat) )
            end do
            centroid2_avrg = centroid2_avrg + centroid2

            do iox = 1, n_ox
               iat = 1 + 3*(iox-1)
               do jox = iox + 1, n_ox

                  ! Total average
                  
                  jat = 1 + 3*(jox-1)
                  vij = centroid(:,jat) - centroid(:,iat)
                  dij2 = dot_product( vij, vij )
                  dij = dsqrt(dij2)
                  rij(iox,jox) = rij(iox,jox) + dij
                  rij2(iox,jox) = rij2(iox,jox) + dij2

                  ! Block average
                  iblock = int( (n-nSteps_equi-1) / dblock ) + 1
                  rij_block(iox,jox,iblock) = rij_block(iox,jox,iblock) + dij
                  rij2_block(iox,jox,iblock) = rij2_block(iox,jox,iblock) + dij2
                  
               end do
            end do
            
            total_energy = kinetic_energy + potential_energy
            qTotalEnergy = qKineticEnergy + qPotentialEnergy
            qTotalEnergy_v2 = qKineticEnergy_v2 + qPotentialEnergy

            temperature_avrg = temperature_avrg + temperature
            kinetic_avrg = kinetic_avrg + kinetic_energy
            potential_avrg = potential_avrg + potential_energy
            totalenergy_avrg = totalenergy_avrg + total_energy
            qtotalenergy_avrg = qtotalenergy_avrg + qTotalEnergy
            qtotalenergy_v2_avrg = qtotalenergy_v2_avrg + qTotalEnergy_v2         
            qkinetic_avrg = qkinetic_avrg + qKineticEnergy
            qkinetic_v2_avrg = qkinetic_v2_avrg + qKineticEnergy_v2
            qpotential_avrg = qpotential_avrg + qPotentialEnergy

            iblock = int( (n-nSteps_equi-1) / dblock ) + 1
            kin_avrg(iblock) = kin_avrg(iblock) + kinetic_energy
            pot_avrg(iblock) = pot_avrg(iblock) + potential_energy
            tote_avrg(iblock) = tote_avrg(iblock) + total_energy
            qkin_avrg(iblock) = qkin_avrg(iblock) + qKineticEnergy
            qkin_v2_avrg(iblock) = qkin_v2_avrg(iblock) + qKineticEnergy_v2
            qpot_avrg(iblock) = qpot_avrg(iblock) + qPotentialEnergy
            qtote_avrg(iblock) = qtote_avrg(iblock) + qTotalEnergy
            qtote_v2_avrg(iblock) = qtote_v2_avrg(iblock) + qTotalEnergy_v2
            r2_avrg(iblock) = r2_avrg(iblock) + radius2
            rmax2_avrg(iblock) = rmax2_avrg(iblock) + radmax2
            doh2_avrg(iblock) = doh2_avrg(iblock) + roh2
            hoh_avrg(iblock) = hoh_avrg(iblock) + theta
            ctr_avrg(:,:,iblock) = ctr_avrg(:,:,iblock) + centroid(:,:)
            ctr2_avrg(:,iblock) = ctr2_avrg(:,iblock) + centroid2(:)
            g2_avrg(:,iblock) = g2_avrg(:,iblock) + gyration2(:)
            
            if ( mod(n, nperiod) .eq. 0 ) then

               ! current values
            
               !write(*,'(i10,10f17.8)') nStep, temperature, kinetic_energy, potential_energy, &
               !     total_energy, qKineticEnergy, qKineticEnergy_v2, qPotentialEnergy, &
               !     qTotalEnergy, qTotalEnergy_v2

               ! accumulated averages

               write(20,'(i10,10f17.8)') nStep, &
                    temperature_avrg / dfloat(n-nSteps_equi), &
                    kinetic_avrg / dfloat(n-nSteps_equi), &
                    potential_avrg / dfloat(n-nSteps_equi), &
                    totalenergy_avrg / dfloat(n-nSteps_equi), &
                    qkinetic_avrg / dfloat(n-nSteps_equi), &
                    qkinetic_v2_avrg / dfloat(n-nSteps_equi), &
                    qpotential_avrg / dfloat(n-nSteps_equi), &
                    qtotalenergy_avrg / dfloat(n-nSteps_equi), &
                    qtotalenergy_v2_avrg / dfloat(n-nSteps_equi)
               
            end if

         end if ! end if n > nSteps_prod

         if ( mod(n,200000) .eq. 0 ) then
            call writeRestartChain( state_chain, newrestartFile, nStep )
         end if

      end do ! end of dynamics

      close(88)
      
! normalize momenta distribution and averages
      
      ! norm_po = zero
      ! norm_ph = zero
      ! norm_ph2o = zero
      ! do j = 1, nbinp
      !    norm_po = norm_po + pdistro_o_avrg(j)*deltap_o
      !    norm_ph = norm_ph + pdistro_h_avrg(j)*deltap_h
      !    norm_ph2o = norm_ph2o + pdistro_h2o_avrg(j)*deltap_h2o
      ! enddo
      ! pdistro_o_avrg = pdistro_o_avrg / norm_po
      ! pdistro_h_avrg = pdistro_h_avrg / norm_ph
      ! pdistro_h2o_avrg = pdistro_h2o_avrg / norm_ph2o

      gyration2_avrg = gyration2_avrg / dfloat(nSteps_prod)
      radius2_avrg = radius2_avrg / dfloat(nSteps_prod)
      radmax2_avrg = radmax2_avrg / dfloat(nSteps_prod)      
      centroid_avrg = centroid_avrg / dfloat(nSteps_prod)
      centroid2_avrg = centroid2_avrg / dfloat(nSteps_prod)
      roh2_avrg = roh2_avrg / dfloat(nSteps_prod)
      theta_avrg = theta_avrg / dfloat(nSteps_prod)

      rij = rij / dfloat(nSteps_prod)
      rij2 = rij2 / dfloat(nSteps_prod)
      
      temperature_avrg = temperature_avrg / dfloat(nSteps_prod)
      kinetic_avrg = kinetic_avrg / dfloat(nSteps_prod)
      potential_avrg = potential_avrg / dfloat(nSteps_prod)
      totalenergy_avrg = totalenergy_avrg / dfloat(nSteps_prod)
      qtotalenergy_avrg = qtotalenergy_avrg / dfloat(nSteps_prod)
      qtotalenergy_v2_avrg = qtotalenergy_v2_avrg / dfloat(nSteps_prod)
      qkinetic_v2_avrg = qkinetic_v2_avrg / dfloat(nSteps_prod)                  
      qkinetic_avrg = qkinetic_avrg / dfloat(nSteps_prod)
      qpotential_avrg = qpotential_avrg / dfloat(nSteps_prod)

      kin_avrg = kin_avrg / dblock
      pot_avrg = pot_avrg / dblock
      tote_avrg = tote_avrg / dblock
      qkin_avrg = qkin_avrg / dblock
      qkin_v2_avrg = qkin_v2_avrg / dblock
      qpot_avrg = qpot_avrg / dblock
      qtote_avrg = qtote_avrg / dblock
      qtote_v2_avrg = qtote_v2_avrg / dblock
      r2_avrg = r2_avrg / dblock
      rmax2_avrg = rmax2_avrg / dblock            
      doh2_avrg = doh2_avrg / dblock
      hoh_avrg = hoh_avrg / dblock
      ctr_avrg = ctr_avrg / dblock
      ctr2_avrg = ctr2_avrg / dblock
      g2_avrg = g2_avrg / dblock

      rij_block = rij_block / dblock
      rij2_block = rij2_block / dblock      

      do k = 1, 3
         density_o_avrg(k,:,:) = density_o_avrg(k,:,:) / sum( density_o_avrg(k,:,:) )
         density_h_avrg(k,:,:) = density_h_avrg(k,:,:) / sum( density_h_avrg(k,:,:) )
      end do

      rhor_o_avrg(:) = rhor_o_avrg(:) / sum( rhor_o_avrg(:) )
      rhor_h_avrg(:) = rhor_h_avrg(:) / sum( rhor_h_avrg(:) )
      rhor_h2o_avrg(:) = rhor_h2o_avrg(:) / sum( rhor_h2o_avrg(:) )

      goo_avrg(:) = goo_avrg(:) / sum( goo_avrg(:) )
      ghh_avrg(:) = ghh_avrg(:) / sum( ghh_avrg(:) )
      goh_avrg(:) = goh_avrg(:) / sum( goh_avrg(:) )      

! compute statistical errors

      sigma_kinetic = zero
      sigma_potential = zero
      sigma_totalenergy = zero
      sigma_qkinetic = zero
      sigma_qkinetic_v2 = zero
      sigma_qpotential = zero
      sigma_qtotalenergy = zero
      sigma_qtotalenergy_v2 = zero
      sigma_radius2 = zero
      sigma_radmax2 = zero      
      sigma_roh2 = zero
      sigma_theta = zero
      do iblock = 1, nblocks
         sigma_kinetic = sigma_kinetic + ( kin_avrg(iblock) - kinetic_avrg )**two
         sigma_potential = sigma_potential + ( pot_avrg(iblock) - potential_avrg )**two
         sigma_totalenergy = sigma_totalenergy + ( tote_avrg(iblock) - totalenergy_avrg )**two
         sigma_qkinetic = sigma_qkinetic + ( qkin_avrg(iblock) - qkinetic_avrg )**two
         sigma_qkinetic_v2 = sigma_qkinetic_v2 + ( qkin_v2_avrg(iblock) - qkinetic_v2_avrg )**two
         sigma_qpotential = sigma_qpotential + ( qpot_avrg(iblock) - qpotential_avrg )**two
         sigma_qtotalenergy = sigma_qtotalenergy + ( qtote_avrg(iblock) - qtotalenergy_avrg )**two
         sigma_qtotalenergy_v2 = sigma_qtotalenergy_v2 + &
              ( qtote_v2_avrg(iblock) - qtotalenergy_v2_avrg )**two
         sigma_radius2 = sigma_radius2 + ( r2_avrg(iblock) - radius2_avrg )**two
         sigma_radmax2 = sigma_radmax2 + ( rmax2_avrg(iblock) - radmax2_avrg )**two         
         sigma_roh2 = sigma_roh2 + ( doh2_avrg(iblock) - roh2_avrg )**two
         sigma_theta = sigma_theta + ( hoh_avrg(iblock) - theta_avrg )**two
      end do
      sigma_kinetic = dsqrt( sigma_kinetic / dfloat(nblocks) ) 
      sigma_potential = dsqrt( sigma_potential / dfloat(nblocks) )
      sigma_totalenergy = dsqrt( sigma_totalenergy / dfloat(nblocks) )
      sigma_qkinetic = dsqrt( sigma_qkinetic / dfloat(nblocks) )
      sigma_qkinetic_v2 = dsqrt( sigma_qkinetic_v2 / dfloat(nblocks) )
      sigma_qpotential = dsqrt( sigma_qpotential / dfloat(nblocks) )
      sigma_qtotalenergy = dsqrt( sigma_qtotalenergy / dfloat(nblocks) )
      sigma_qtotalenergy_v2 = dsqrt( sigma_qtotalenergy_v2 / dfloat(nblocks) )
      sigma_radius2 = dsqrt( sigma_radius2 / dfloat(nblocks) )
      sigma_radmax2 = dsqrt( sigma_radmax2 / dfloat(nblocks) )      
      sigma_roh2 = dsqrt( sigma_roh2 / dfloat(nblocks) )
      sigma_theta = dsqrt( sigma_theta / dfloat(nblocks) )
      
! outputs

      ! ! momenta distribution

      ! open(unit=21, file="pdistro-o.dat")
      ! open(unit=22, file="pdistro-h.dat")
      ! open(unit=23, file="pdistro-h2o.dat")

      ! ! p distribution of oxygens
      ! do j = 1, nbinp
      !    pc = pmin_o + (dfloat(j)-half)*deltap_o
      !    gaussian = one / dsqrt( two * pi ) / sigma_o * dexp( -pc*pc /two /sigma_o**two )
      !    write(21,*) pc, pdistro_o_avrg(j), gaussian
      ! end do

      ! ! p distribution of hydrogens
      ! do j = 1, nbinp
      !    pc = pmin_h + (dfloat(j)-half)*deltap_h
      !    gaussian = one / dsqrt( two * pi ) / sigma_h * dexp( -pc*pc /two /sigma_h**two )
      !    write(22,*) pc, pdistro_h_avrg(j), gaussian
      ! end do
      
      ! ! p distribution of water molecules
      ! do j = 1, nbinp
      !    pc = pmin_h2o + (dfloat(j)-half)*deltap_h2o
      !    gaussian = one / dsqrt( two * pi ) / sigma_h2o * dexp( -pc*pc /two /sigma_h2o**two )
      !    write(23,*) pc, pdistro_h2o_avrg(j), gaussian
      ! end do
         
      ! close(21)
      ! close(22)
      ! close(23)

      ! Final averages

      write(20,*)
      write(20,*) 'Temperature average:', temperature_avrg / temperature_target
      !write(*,*) 'Kinetic energy average:', kinetic_avrg / ( three / two * kT ) / dfloat(nReplicas)
      !write(*,*) 'Potential energy average:', potential_avrg / ( three / two * kT ) / &
      !     dfloat(nReplicas)

      if (nReplicas.gt.1) then
         write(20,*)
         write(20,*) 'Quantum energies per molecule (in kcalmol)'
         write(20,*) 'Total:', qtotalenergy_avrg / dfloat(n_mols) / kcalmol_to_au, &
              sigma_qtotalenergy / dfloat(n_mols) / kcalmol_to_au
         write(20,*) 'Total_v2:', qtotalenergy_v2_avrg / dfloat(n_mols) / kcalmol_to_au, &
              sigma_qtotalenergy_v2 / dfloat(n_mols) / kcalmol_to_au         
         write(20,*) 'Kinetic:', qkinetic_avrg / dfloat(n_mols) / kcalmol_to_au, &
              sigma_qkinetic / dfloat(n_mols) / kcalmol_to_au
         write(20,*) 'Kinetic_v2:', qkinetic_v2_avrg / dfloat(n_mols) / kcalmol_to_au, &
              sigma_qkinetic_v2 / dfloat(n_mols) / kcalmol_to_au
         write(20,*) 'Potential:', qpotential_avrg / dfloat(n_mols) / kcalmol_to_au, &
              sigma_qpotential / dfloat(n_mols) / kcalmol_to_au
      end if
      write(20,*)
      write(20,*) 'Extended system classical energies per molecule (in kcalmol)'
      write(20,*) 'Total:', totalenergy_avrg / dfloat(n_mols) / kcalmol_to_au, &
           sigma_totalenergy / dfloat(n_mols) / kcalmol_to_au
      write(20,*) 'Kinetic:', kinetic_avrg / dfloat(n_mols) / kcalmol_to_au, &
           sigma_kinetic / dfloat(n_mols) / kcalmol_to_au
      write(20,*) 'Potential:', potential_avrg / dfloat(n_mols) / kcalmol_to_au, &
           sigma_potential / dfloat(n_mols) / kcalmol_to_au

      ! squared mean radius and gyration radius
      
      write(20,*)
      write(20,*) 'Squared mean radius (in Angstrom)'
      write(20,*) dsqrt(radius2_avrg) / ang_to_au, &
           half / dsqrt(radius2_avrg) * sigma_radius2 / ang_to_au
      write(20,*) 'Root mean square max radius of the cluster (in Angstrom)'
      write(20,*) dsqrt(radmax2_avrg) / ang_to_au, &
           half / dsqrt(radmax2_avrg) * sigma_radmax2 / ang_to_au            

      if ( nReplicas.gt.1 ) then

         ! Mean gyration radius
         
         gyration2_o_avrg = zero
         gyration2_h_avrg = zero
         do imol = 1, n_mols
            do j = 1, 3
               iat = j + 3*(imol-1)
               if (j.eq.1) then
                  gyration2_o_avrg = gyration2_o_avrg + gyration2_avrg(iat)
               else
                  gyration2_h_avrg = gyration2_h_avrg + gyration2_avrg(iat)
               end if
            end do
         end do
         gyration2_o_avrg = gyration2_o_avrg / dfloat( n_mols )
         gyration2_h_avrg = gyration2_h_avrg / ( two*dfloat(n_mols) )

         ! Uncertainty
         
         sigma_gyration2_o = zero
         sigma_gyration2_h = zero

         do iblock = 1, nblocks

            ! Gyration radius of each block

            gyration2_o = zero
            gyration2_h = zero            
            do imol = 1, n_mols
               do j = 1, 3
                  iat = j + 3*(imol-1)
                  if (j.eq.1) then
                     gyration2_o = gyration2_o + g2_avrg(iat,iblock)                     
                  else
                     gyration2_h = gyration2_h + g2_avrg(iat,iblock)                  
                  end if
               end do
            end do
            gyration2_o = gyration2_o / dfloat( n_mols )
            gyration2_h = gyration2_h / ( two*dfloat(n_mols) )
            
            ! Accumulate deviation

            sigma_gyration2_o = sigma_gyration2_o + ( gyration2_o - gyration2_o_avrg )**two
            sigma_gyration2_h = sigma_gyration2_h + ( gyration2_h - gyration2_h_avrg )**two         
            
         end do
         sigma_gyration2_o = dsqrt( sigma_gyration2_o / dfloat(nblocks) )
         sigma_gyration2_h = dsqrt( sigma_gyration2_h / dfloat(nblocks) )                  
         
         ! Error(R) = 1 / (2 * R ) * Error(R^2)
         ! Error( sqrt(R^2) ) = Error(R) = 1 / (2 * R ) * Error(R^2)
         write(20,*) 'Gyration radius of O and H (in Angstrom)'
         write(20,*) dsqrt(gyration2_o_avrg) / ang_to_au, &
              half / dsqrt(gyration2_o_avrg) * sigma_gyration2_o / ang_to_au, &
              dsqrt(gyration2_h_avrg) / ang_to_au, &
              half / dsqrt(gyration2_h_avrg) * sigma_gyration2_h / ang_to_au

      end if

      ! Mean squared deviation from the averaged centroid position

      msd_ox_avrg = zero
      msd_hy_avrg = zero

      do imol = 1, n_mols
         do j = 1,3
            iat = j + 3*(imol-1)
            if (j.eq.1) then
               msd_ox_avrg = msd_ox_avrg + centroid2_avrg(iat) - &
                    dot_product( centroid_avrg(:,iat), centroid_avrg(:,iat) )
            else
               msd_hy_avrg = msd_hy_avrg + centroid2_avrg(iat) - &
                    dot_product( centroid_avrg(:,iat), centroid_avrg(:,iat) )               
            end if
         end do
      end do
      msd_ox_avrg = msd_ox_avrg / dfloat(n_mols)
      msd_hy_avrg = msd_hy_avrg / dfloat(2*n_mols)

      sigma_msd_ox = zero
      sigma_msd_hy = zero      
      do iblock = 1, nblocks

         ! Compute the MSD of each block

         msd_ox = zero
         msd_hy = zero
         do imol = 1, n_mols
            do j = 1, 3
               iat = j + 3*(imol-1)
               if (j.eq.1) then
                  msd_ox = msd_ox + ctr2_avrg(iat,iblock) - &
                       dot_product( ctr_avrg(:,iat,iblock), ctr_avrg(:,iat,iblock) )
               else
                  msd_hy = msd_hy + ctr2_avrg(iat,iblock) - &
                       dot_product( ctr_avrg(:,iat,iblock), ctr_avrg(:,iat,iblock) )        
               end if
            end do
         end do
         msd_ox = msd_ox / dfloat(n_mols)
         msd_hy = msd_hy / dfloat(2*n_mols)

         ! Accumulate the deviation

         sigma_msd_ox = sigma_msd_ox + ( msd_ox - msd_ox_avrg )**two
         sigma_msd_hy = sigma_msd_hy + ( msd_hy - msd_hy_avrg )**two         

      end do
      sigma_msd_ox = dsqrt( sigma_msd_ox / dfloat(nblocks) )
      sigma_msd_hy = dsqrt( sigma_msd_hy / dfloat(nblocks) )            
      
      ! Error(MSD) = 1 / (2 * MSD ) * Error(MSD^2)
      ! Error( sqrt(MSD^2) ) = Error(MSD) = 1 / (2 * MSD ) * Error(MSD^2)      
      write(20,*)
      write(20,*) 'Root Mean Squared Deviation of O and H (in Angstrom)'
      write(20,*) dsqrt(msd_ox_avrg) / ang_to_au, &
           half / dsqrt(msd_ox_avrg) * sigma_msd_ox / ang_to_au, &
           dsqrt(msd_hy_avrg) / ang_to_au, &
           half / dsqrt(msd_hy_avrg) * sigma_msd_hy / ang_to_au 

      ! rOH distance, HOH angle and dipole moment
      
      write(20,*)
      write(20,*) 'Mean rOH distance (in Angstrom)'
      write(20,*) dsqrt(roh2_avrg) / ang_to_au, &
           half / dsqrt(roh2_avrg) * sigma_roh2 / ang_to_au

      write(20,*) 'Mean HOH angle (in degrees)'
      write(20,*) theta_avrg * rad_to_deg, sigma_theta * rad_to_deg

      ! Lindemann index

      lindemann = zero
      do iox = 1, n_ox
         do jox = iox+1, n_ox
            dij = rij(iox,jox)
            dij2 = rij2(iox,jox)
            lindemann = lindemann + dsqrt( dij2 - dij*dij ) / dij
         end do
      end do
      lindemann = lindemann * two / n_ox / (n_ox -1)

      sigma_linde = zero
      do iblock = 1, nblocks

         ! lindemann of iblock
         linde = zero
         do iox = 1, n_ox
            do jox = iox+1, n_ox
               dij = rij_block(iox,jox,iblock)
               dij2 = rij2_block(iox,jox,iblock)
               linde = linde + dsqrt( dij2 - dij*dij ) / dij
            end do
         end do
         linde = linde * two / n_ox / (n_ox -1)
         
         ! Accumulate the deviation
         sigma_linde = sigma_linde + ( linde - lindemann )**two         

      end do
      sigma_linde = dsqrt( sigma_linde / dfloat(nblocks) )                  
      
      write(20,*)
      write(20,*) 'Lindemann index'
      write(20,*) lindemann, sigma_linde
      
      ! Averaged centroid 
      
      open(unit=30, file='centroid.xyz')
      write(30,*) n_atoms
      write(30,*)
      do iat = 1, n_atoms
         if (mod(iat,3).eq.1) then
            !write(30,*) 'O', centroid_avrg(1:3,iat) / ang_to_au ! convert a.u. to Angstrom
            write(30,*) 'O', centroid(1:3,iat) / ang_to_au ! convert a.u. to Angstrom
         else
            !write(30,*) 'H', centroid_avrg(1:3,iat) / ang_to_au ! convert a.u. to Angstrom
            write(30,*) 'H', centroid(1:3,iat) / ang_to_au ! convert a.u. to Angstrom
         end if
      end do
      close(30)

      ! Polymer positions of atoms over different replicas

      if (nReplicas.gt.1) then
      
         call PolymerDistribution( state_chain, polymer )     

         open(unit=31, file='polymer.xyz')
         write(31,*) n_atoms*nReplicas
         write(31,*)
         do i = 1, nReplicas
            do iat = 1, n_atoms
               if (mod(iat,3).eq.1) then
                  write(31,*) 'O', polymer(1:3,iat,i) / ang_to_au ! convert a.u. to Angstrom
               else
                  write(31,*) 'H', polymer(1:3,iat,i) / ang_to_au ! convert a.u. to Angstrom
               end if
            end do
         end do
         close(31)
         
      end if

      ! Spatial distribution of the centroid

      xmax = rmax_histo / dsqrt(three) * ang_to_au
      xmin = -xmax
      dx = (xmax-xmin)/dfloat(nbinx)
      open(unit=32, file='xdistro-o.dat')
      open(unit=33, file='xdistro-h.dat')
      do i1 = 0, nbinx-1
         x1c = ( xmin + ( dfloat(i1)+half )*dx ) / ang_to_au
         do i2 = 0, nbinx-1
            x2c = ( xmin + ( dfloat(i2)+half )*dx ) / ang_to_au
            write(32,'(5f17.8)') x1c, x2c, density_o_avrg(1,i1+1,i2+1), &
                 density_o_avrg(2,i1+1,i2+1), density_o_avrg(3,i1+1,i2+1)
            write(33,'(5f17.8)') x1c, x2c, density_h_avrg(1,i1+1,i2+1), &
                 density_h_avrg(2,i1+1,i2+1), density_h_avrg(3,i1+1,i2+1)
         end do ! end iy

         write(32,*)
         write(33,*)
         
      end do ! end ix
      close(32)
      close(33)

      rmax_rho = rmax_histo * ang_to_au
      dr_rho = rmax_rho / dfloat(nbinr)
      open(unit=34, file='rdistro.dat')
      do ir = 0, nbinr-1
         rc = ( dfloat(ir) + half ) * dr_rho / ang_to_au
         write(34,'(4f17.8)') rc, rhor_o_avrg(ir+1), rhor_h_avrg(ir+1), rhor_h2o_avrg(ir+1)
      end do
      close(34)

      rmax_g = rmax_histo * ang_to_au
      dr_g = rmax_g / dfloat(nbinr)
      open(unit=35, file='pairdistro.dat')
      do ir = 0, nbinr-1
         rc = ( dfloat(ir) + half ) * dr_g / ang_to_au
         write(35,'(4f17.8)') rc, goo_avrg(ir+1), ghh_avrg(ir+1), goh_avrg(ir+1)
      end do
      close(35)   
      
! write restart file
      
      call writeRestartChain( state_chain, newrestartFile, nStep )
      
! deallocate

      deallocate( mass )
      deallocate( position0 )
      !deallocate( momentum )
      deallocate( centroid )
      deallocate( centroid0 )
      deallocate( centroid2 )      
      deallocate( centroid_avrg )
      deallocate( centroid2_avrg )
      deallocate( polymer )      
      deallocate( gyration2 )
      deallocate( gyration2_avrg )
      deallocate( ctr_avrg )
      deallocate( ctr2_avrg )
      deallocate( g2_avrg )
      deallocate( rij )
      deallocate( rij2 )
      deallocate( rij_block )
      deallocate( rij2_block )      

      deallocate( rhor_o )
      deallocate( rhor_h )
      deallocate( rhor_h2o )
      deallocate( rhor_o_avrg )
      deallocate( rhor_h_avrg )
      deallocate( rhor_h2o_avrg )
      deallocate( goo )
      deallocate( ghh )
      deallocate( goh )
      deallocate( goo_avrg )
      deallocate( ghh_avrg )
      deallocate( goh_avrg )
      
      call deleteStateChain( state_chain )

! end program

      call cpu_time(tfinish_cpu)
      cputime = tfinish_cpu - tstart_cpu
      call system_clock( count_1, count_rate, count_max )
      walltime = dfloat(count_1)/dfloat(count_rate) - dfloat(count_0)/dfloat(count_rate)
      !walltime = omp_get_wtime() - walltime
      write(*,*)
      write(*,*) 'CPU time (seconds, minutes, hours):'
      write(*,'(3f10.3)') cputime, cputime/60.0_dp, cputime/3600.0_dp
      write(20,*)
      write(20,*) 'CPU time (seconds, minutes, hours):'
      write(20,'(3f10.3)') cputime, cputime/60.0_dp, cputime/3600.0_dp
      write(*,*)
      write(*,*) 'Wall clock time (seconds, minutes, hours):'
      write(*,'(3f10.3)') walltime, walltime/60.0_dp, walltime/3600.0_dp
      write(20,*)
      write(20,*) 'Wall clock time (seconds, minutes, hours):'
      write(20,'(3f10.3)') walltime, walltime/60.0_dp, walltime/3600.0_dp            

      close(20) ! close logfile            

contains
  
!*****************************************************************************
subroutine computeForces( state_chain, temperature_target, gyration2, centroid )

!*****************************************************************************
!
!  Project: MolecularDynamics
!
!  Created on: March 2019 by Alfonso Gijón
!
! ---
! Copyright (C) 2018       Eduardo R. Hernandez - ICMM
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! ---
!
! This subroutine receives a chain of states to compute intra and inter-replica
! energies and forces, and update them.
!
! Also compute the gyration radius and the centroid of the atoms
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
      implicit none

!*****************************************************************************
!  arguments

      type (StateChain_type), intent (INOUT) :: state_chain

      double precision, intent(IN) :: temperature_target

      double precision, dimension(:), intent (OUT) :: gyration2

      double precision, dimension (:,:), intent(OUT) :: centroid

!*****************************************************************************
!  local variables

      integer i, j, inext
      integer imol, iat

      integer nAtoms
      integer nReplicas
      integer nmols

      double precision :: beta, beta2

      double precision potential_energy
      double precision qPotentialEnergy
      double precision qKineticEnergy
      double precision qKineticEnergy_v2
      
      double precision vintra_total
      double precision vspring, vspring_total
      double precision vspring_v2, vspring_total_v2
      double precision :: kfactor, kspring
      
      double precision, dimension (3) :: rspring
      double precision, dimension (3) :: fspring
      double precision, dimension (3) :: deltar

      double precision, dimension (3) :: a, b, c
      double precision, allocatable, dimension (:,:,:) :: position
      double precision, allocatable, dimension (:,:,:) :: force
      double precision, allocatable, dimension (:) :: mass
      double precision, allocatable, dimension (:) :: vintra      

!*****************************************************************************
!  start subroutine

      call getStateChain( state_chain, nReplicas = nReplicas, nAtoms = nAtoms )

      nmols = nAtoms/3
      
      allocate( position( 3, nAtoms, nReplicas ) )
      allocate( force( 3, nAtoms, nReplicas ) )
      allocate( mass( nAtoms ) )
      allocate( vintra( nReplicas ) )      

      call getStateChain( state_chain, mass = mass )
      
! compute intra-replica energy and forces
      
      vintra = zero
      position = zero
      force = zero

      !$OMP PARALLEL &
      !$OMP PRIVATE( i )
      !$OMP DO
      do i = 1, nReplicas

         ! get position from each replica of the chain

         call getStateChain( state_chain, iReplica = i, position = position(:,:,i) )

         call evaluate_energy_SPCF( nmols, position(:,:,i), vintra(i), &
              force(:,:,i) )

      end do
      !$OMP END DO
      !$OMP END PARALLEL

      ! compute centroid position and total intra-replica energy

      vintra_total = zero
      centroid = zero         
      
      do i = 1, nReplicas

         vintra_total = vintra_total + vintra(i)
         centroid(:,:) = centroid(:,:) + position(:,:,i)
         
      end do
      
      centroid = centroid / dfloat(nReplicas)

! compute inter-replica energy and forces
      
      ! Calculate vspring, using the internal forces and the path integral virial theorem

      vspring_total_v2 = zero
      
      gyration2 = zero

      do i = 1, nReplicas

         vspring_v2 = zero

         do j = 1, nAtoms

            !deltar = position(:,j,i) - position(:,j,1)
            deltar = position(:,j,i) - centroid(:,j)
            vspring_v2 = vspring_v2 + dot_product( deltar, force(:,j,i) )
            
            ! gyration radius squared
            gyration2(j) = gyration2(j) + dot_product( deltar, deltar )
            
         end do

         vspring_total_v2 = vspring_total_v2 + vspring_v2
         
      end do

      vspring_total_v2 = vspring_total_v2 * half / dfloat(nReplicas)

      gyration2 = gyration2 / dfloat(nReplicas)

      ! In the the polymer Hamiltonian, the internal potential is reduced by a 1/nReplicas factor
      
      vintra_total = vintra_total / dfloat(nReplicas)
      force = force / dfloat(nReplicas)
      
      ! add inter-replica energy and forces (in a.u.)

      beta = one / ( boltzmann_k * temperature_target )
      beta2 = beta*beta

      kfactor = dfloat(nReplicas) / beta2

      vspring_total = zero

      do i = 1, nReplicas

         inext = i + 1
         if (inext == nReplicas+1) inext = 1

         vspring = zero

         do j = 1, nAtoms

            kspring = mass(j) * kfactor ! it depends on mass, number of replicas and temperature
            rspring = position(:,j,inext) - position(:,j,i)
            vspring = vspring + half * kspring * dot_product( rspring, rspring )
            
            fspring = kspring * rspring ! factor 2 is compensated with factor 1/2
            force(:,j,i) = force(:,j,i) + fspring
            force(:,j,inext) = force(:,j,inext) - fspring
                  
         end do

         vspring_total = vspring_total + vspring

      end do

      do i = 1, nReplicas

         call setStateChain( state_chain, iReplica = i, force = force(:,:,i) )
         
      end do

      ! update classical potential energy of the chain

      potential_energy = vintra_total + vspring_total
            
      call setStateChain( state_chain, potentialEnergy = potential_energy )
      
      ! update quantum energies of the chain

      qPotentialEnergy = vintra_total
      qKineticEnergy = three * half * dfloat(nReplicas*nAtoms) / beta - vspring_total
      qKineticEnergy_v2 = three * half * dfloat(nAtoms) / beta - vspring_total_v2
      
      call setStateChain( state_chain, qPotentialEnergy = qPotentialEnergy, &
           qKineticEnergy = qKineticEnergy, qKineticEnergy_v2 = qKineticEnergy_v2 )
      
      ! deallocate local variables
      
      deallocate( position )
      deallocate( force )
      deallocate( mass )
      deallocate( vintra )

end subroutine computeForces    

end program PIMD_Water
