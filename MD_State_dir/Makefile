
MODULE = ../MD_State.f90

# This Makefile creates ../MD_State.f90 by catenating the following files
# The first files correspond to public routines in MD_State.f90 module in alphabetic order;
# then the private ones (also in alphabetic order) follow

${MODULE} : MD_State_header \
            createState.f90 \
	    calculateKEChain.f90 \
	    calculateMomentumDistro.f90 \
	    CentroidProperties.f90 \
	    MSD.f90 \
	    CentroidProperties_dipole.f90 \
            deleteState.f90 \
	    deleteStateChain.f90 \
            get.f90 \
            getConfiguration.f90 \
            getMass.f90 \
	    getStateChain.f90 \
            initialiseState.f90 \
            inquire.f90 \
            inquireStateChain.f90 \
            PolymerDistribution.f90 \
            PolymerProperties.f90 \
            set.f90 \
            setForce.f90 \
	    setStateChain.f90 \
            setUnits.f90 \
	    setupStateChain.f90 \
            readRestart.f90 \
            readRestartChain.f90 \
            writeRestart.f90 \
            writeRestartChain.f90 \
            private.txt \
            calculateInternalPressure.f90 \
            calculateKineticEnergy.f90 \
            generateMomenta.f90 \
            initialiseCell.f90 \
            random.f90 \
            setReciprocalCell.f90 \
            transformCoordinates2Lattice.f90 \
            transformCoordinates2Cartesian.f90 \
            transformForces2Lattice.f90 \
            transformForces2Cartesian.f90 \
            transformMomenta2Lattice.f90 \
            transformMomenta2Cartesian.f90 \
            MD_State_tail
	@cat $^ > $@
