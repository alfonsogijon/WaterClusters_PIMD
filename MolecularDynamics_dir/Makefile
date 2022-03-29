
MODULE = ../MolecularDynamics.f90

# This Makefile creates ../MolecularDynamics.f90 by catenating the following files
# The first files correspond to public routines in MolecularDynamics.f90 module in alphabetic order;
# then the private ones (also in alphabetic order) follow

${MODULE} : MolecularDynamics_header \
	    BussiParrinello_A.f90 \
	    BussiParrinello_B.f90 \
            MD_setUp.f90 \
            MD_cleanUp.f90 \
	    MDChain_setup.f90 \
            MD_stepA.f90 \
            MD_stepB.f90 \
	    NPChain_StateChain_A.f90 \
	    NPChain_StateChain_B.f90 \
            velocityVerlet_StateChain_A.f90 \
            velocityVerlet_StateChain_B.f90 \
            private.txt \
            calculateThermostatEnergy.f90 \
            NosePoincareChain_A.f90 \
            NosePoincareChain_B.f90 \
            NPSM_A.f90 \
            NPSM_B.f90 \
	    random.f90 \
            velocityVerlet_A.f90 \
            velocityVerlet_B.f90 \
            updateCell.f90 \
            MolecularDynamics_tail
	@cat $^ > $@
