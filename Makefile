
NAME1 = PIMD_Water
NAME2 = PIMD_Water_qTIP4Pf

all: modules ${NAME1} ${NAME2}

FC = gfortran
OFLAGS = -O2
DFLAGS = #-pg -fcheck='all' #-Wall -Wextra -g -pedantic-errors -Wsurprising -Wunderflow -Wintrinsic-shadow -fbacktrace -fimplicit-none
OMPFLAGS = -fopenmp
#LFLAGS = 

MDOFILES = SPCF.o \
           qTIP4Pf.o \
           MD_precision_module.o \
           MD_constants_module.o \
           MD_State.o \
           MolecularDynamics.o

${NAME1} : ${MDOFILES} PIMD_Water.o
	${FC} ${OFLAGS} ${OMPFLAGS} -o PIMD_Water.x ${DFLAGS} PIMD_Water.o ${MDOFILES}

${NAME2} : ${MDOFILES} PIMD_Water_qTIP4Pf.o
	${FC} ${OFLAGS} ${OMPFLAGS} -o PIMD_Water_qTIP4Pf.x ${DFLAGS} PIMD_Water_qTIP4Pf.o ${MDOFILES}

MolecularDynamics.o: MolecularDynamics.f90 \
                     MD_State.o
	${FC} -c ${OFLAGS} ${DFLAGS} ${OMPFLAGS} $<

SPCF.o: SPCF.f90 
	${FC} -c ${OFLAGS} ${DFLAGS} ${OMPFLAGS} $<

qTIP4Pf.o: qTIP4Pf.f90 
	${FC} -c ${OFLAGS} ${DFLAGS} ${OMPFLAGS} $<

MD_State.o: MD_State.f90 \
            MD_precision_module.o \
            MD_constants_module.o
	${FC} -c ${OFLAGS} ${DFLAGS} ${OMPFLAGS} $<

PIMD_Water.o: PIMD_Water.f90 \
           SPCF.o \
           MolecularDynamics.o \
           MD_State.o
	${FC} -c ${OFLAGS} ${DFLAGS} ${OMPFLAGS} $<

PIMD_Water_qTIP4Pf.o: PIMD_Water_qTIP4Pf.f90 \
           qTIP4Pf.o \
           MolecularDynamics.o \
           MD_State.o
	${FC} -c ${OFLAGS} ${DFLAGS} ${OMPFLAGS} $<

MD_precision_module.o: MD_precision_module.f90
	${FC} -c ${OFLAGS} ${DFLAGS} $<

MD_constants_module.o: MD_constants_module.f90 \
                       MD_precision_module.o
	${FC} -c ${OFLAGS} ${DFLAGS} $<

MODULE_DIRS = MD_State_dir \
              MolecularDynamics_dir \

.PHONY: module_dirs ${MODULE_DIRS}

modules:
	for dir in ${MODULE_DIRS}; do \
	    $(MAKE) -s -C $$dir; \
	done

clean:
	rm *.o *.mod
