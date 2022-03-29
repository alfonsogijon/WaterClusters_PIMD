!
! ======================================================
! =         Metodes  Numerics per a la Fisica          =
! =          R. Guardiola  E. Higon  J. Ros            =
! =          Universitat de Valencia (1995)            =
! ======================================================
!
  real(dp) function randini(dseed)
! Generador pres de Carlson, codi modificat
! Crida per a iniciar el generador

    implicit none

! arguments

    real(dp), intent(INOUT) :: dseed

! local variables

    integer i
    real(dp) :: pa1
    
! start function

    iset = 0
    
    ! En cas d'entrada erronia, assegurem sequencia
    if (dseed .GT. vmodul1 .OR. dseed .LT. 1.d0) dseed=91063
    DO i=1,55
       dseed=MOD(vmultp*dseed,vmodul1)
       pa1=dseed*vmodul1
       dseed=MOD(vmultp*dseed,vmodul1)
       y(i)=(pa1+dseed)
    end do
    jang=24
    kang=55
    randini=0
    return
  end function randini
!----------------------------------------------------------------------------


  REAL*8 FUNCTION RANDU()
! Crida per a generar un nombre a l'atzar
!    

    implicit none

    y(kang)=MOD(y(kang)+y(jang),vmodul)
    RANDU=y(kang)/vmodul
    jang=jang-1
    kang=kang-1
    IF (jang .EQ. 0) jang=55
    IF (kang .EQ. 0) kang=55
    return
  END FUNCTION RANDU
!-----------------------------------------------------------------------------

  
  real(dp) FUNCTION GASDEV()

    implicit none

    real(dp) :: rr1, rr2, v1, v2, rsq, fac
    
    IF (ISET.EQ.0) THEN
1      RR1=RANDU()
       RR2=RANDU()
       V1=2.0*RR1-1.0
       V2=2.0*RR2-1.0
       RSQ=V1**2.+V2**2.
       IF(RSQ.GE.1.0.OR.RSQ.EQ.0.0) GO TO 1
       FAC=SQRT(-2.0*LOG(RSQ)/RSQ)
       GSET=V1*FAC
       GASDEV=V2*FAC
       ISET=1
    ELSE
       GASDEV=GSET
       ISET=0
    ENDIF
!
    RETURN
  END FUNCTION GASDEV
