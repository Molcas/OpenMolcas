************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE FOPAB(FIFA,NFIFA,IBRA,IKET,FOPEL)
      use definitions, only: iwp, wp
      use constants, only: Zero, One, Two
      use gugx, only: SGS, L2ACT, EXS, CIS
      use caspt2_global, only: LUCIEX, IDCIEX
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NSYM,NORB,NISH,ISCF,NCONF,STSYM,NASH,
     &                         NAES
      IMPLICIT None

      integer(kind=iwp), intent(in):: NFIFA, IBRA, IKET
      real(kind=wp), intent(in):: FIFA(NFIFA)
      real(kind=wp), intent(out):: FOPEL

* Purely local array, offsets:
      integer(kind=iwp) IOFF(8)
      integer(kind=iwp) :: nLev
      real(kind=wp), allocatable:: BRA(:), KET(:), SGM(:)
      integer(kind=iwp) IOF,ISYM,IFTEST,IJ,I,ID,II,ISCR,IST,ISU,IT,
     &                  ITABS,ITTOT,ITUTOT,IU,IUABS,IUTOT,J,LEVT,LEVU,
     &                  NI
      real(kind=wp) ESUM,OCC,EINACT,FTU,TRC
      real(kind=wp), external:: DDot_

      nLev = SGS%nLev

* Procedure for computing one matrix element of the Fock matrix in the
* basis of the CASSCF states: <BRA|FOP|KET>
* In: The (possibly average) Fock matrix, active indices only, over the
* original CASSCF orbitals and the indices of the two states


* Offset table for accessing FIFA array:
      IOF=0
      DO ISYM=1,NSYM
        IOFF(ISYM)=IOF
        IOF=IOF+(NORB(ISYM)*(NORB(ISYM)+1))/2
      END DO

      IFTEST=0
      IF (IFTEST.GT.0) THEN
        WRITE(6,*)' The FIFA array:'
        DO ISYM=1,NSYM
          IJ=IOFF(ISYM)+1
          DO I=1,NORB(ISYM)
            WRITE(6,'(1x,5F16.8)')(FIFA(IJ+J),J=0,I-1)
            IJ=IJ+I
          END DO
        END DO
      END IF

* Specialized code for Closed-shell or Hi-spin HF:
* Sum up diagonal elements of FIFA times occ. number
* FIXME: This only works for diagonal elements, thus
* in the case of XMS this will not work...
      IF (ISCF.EQ.1 .OR. ISCF.EQ.2) THEN
        ESUM=Zero
        IF (IBRA.EQ.IKET) THEN
          DO ISYM=1,NSYM
            OCC=Two
            DO I=1,NISH(ISYM)
              ESUM=ESUM+OCC*FIFA(IOFF(ISYM)+(I*(I+1))/2)
            END DO
            IF (ISCF.eq.2) OCC=One
            DO J=1,NASH(ISYM)
              I=NISH(ISYM)+J
              ESUM=ESUM+OCC*FIFA(IOFF(ISYM)+(I*(I+1))/2)
            END DO
          END DO
        ELSE
          WRITE(6,*) ' Warning: neglecting the off-diagonal entries'
          WRITE(6,*) ' of H0, XMS will be equal to MS!'
        END IF
        FOPEL=ESUM
        RETURN
      END IF

* General CASSCF or RASSCF case:
* Sum up trace of FIFA over inactive orbitals only:
      TRC=Zero
      DO ISYM=1,NSYM
        DO I=1,NISH(ISYM)
          II=IOFF(ISYM)+(I*(I+1))/2
          TRC=TRC+FIFA(II)
        END DO
      END DO
* Contribution from inactive orbitals:
      EINACT=Two*TRC

      IF (IFTEST.GT.0) THEN
        WRITE(6,*)' Energy contrib from inactive orbitals:',EINACT
      END IF

* Allocate arrays for ket and bra wave functions
      CALL mma_allocate(BRA,NCONF,Label='BRA')
      CALL mma_allocate(KET,NCONF,Label='KET')
* Allocate array for sigma = Fock operator acting on ket:
      CALL mma_allocate(SGM,NCONF,LABEL='SGM')

* Load ket wave function
      ID=IDCIEX
      DO I=1,IKET-1
        CALL DDAFILE(LUCIEX,0,KET,NCONF,ID)
      END DO
      CALL DDAFILE(LUCIEX,2,KET,NCONF,ID)

      IF (IFTEST.GT.0) THEN
        WRITE(6,*)' IKET:',IKET
        WRITE(6,*)' Ket CI array:'
        ISCR=MIN(NCONF,20)
        WRITE(6,'(1x,5F16.8)')(KET(I),I=1,ISCR)
      END IF

* Compute (lowering part of) FIFA operator acting on
* the ket wave function.
      CALL DCOPY_(NCONF,[Zero],0,SGM,1)
      DO LEVU=1,NLEV
        IUABS=L2ACT(LEVU)
        ISU=SGS%ISM(LEVU)
        IU=IUABS-NAES(ISU)
        NI=NISH(ISU)
        IUTOT=NI+IU
        DO LEVT= 1,LEVU
          IF(SGS%ISM(LEVT).NE.ISU) GOTO 10
          ITABS=L2ACT(LEVT)
          IST=ISU
          IT=ITABS-NAES(IST)
          ITTOT=NI+IT
          ITUTOT=(IUTOT*(IUTOT-1))/2+ITTOT
          IF (ITTOT.GT.IUTOT) ITUTOT=(ITTOT*(ITTOT-1))/2+IUTOT
          FTU=FIFA(IOFF(ISU)+ITUTOT)
          IF(ABS(FTU).LT.1.0D-16) GOTO 10
          CALL SIGMA1(SGS,CIS,EXS,
     &                LEVT,LEVU,FTU,STSYM,KET,SGM)
  10      CONTINUE
        END DO
      END DO
* Add contribution from inactive part:
      CALL DAXPY_(NCONF,EINACT,KET,1,SGM,1)

      IF (IFTEST.GT.0) THEN
        WRITE(6,*)' SGM array from (lowering F)|KET>:'
        WRITE(6,'(1x,5F16.8)')(SGM(I),I=1,NCONF)
      END IF

* Load bra wave function
      ID=IDCIEX
      DO I=1,IBRA-1
        CALL DDAFILE(LUCIEX,0,BRA,NCONF,ID)
      END DO
      CALL DDAFILE(LUCIEX,2,BRA,NCONF,ID)

* Put matrix element into FOPEL:
      FOPEL=DDOT_(NCONF,BRA,1,SGM,1)

      IF (IFTEST.GT.0) THEN
        WRITE(6,*)' FOPEL is now:'
        WRITE(6,'(1x,5f16.8)')(FOPEL)
      END IF

* Compute (strictly lowering part of) FIFA operator acting on |BRA>.
* We are computing contributions <KET|Etu|BRA> with t<u, then
* using them as <BRA|Eut|KET>
* Note that I already have BRA in memory
      CALL DCOPY_(NCONF,[Zero],0,SGM,1)
      DO LEVU=2,NLEV
        IUABS=L2ACT(LEVU)
        ISU=SGS%ISM(LEVU)
        IU=IUABS-NAES(ISU)
        NI=NISH(ISU)
        IUTOT=NI+IU
        DO LEVT= 1,LEVU-1
          IF(SGS%ISM(LEVT).NE.ISU) GOTO 20
          ITABS=L2ACT(LEVT)
          IST=ISU
          IT=ITABS-NAES(IST)
          ITTOT=NI+IT
          ITUTOT=(IUTOT*(IUTOT-1))/2+ITTOT
          IF (ITTOT.GT.IUTOT) ITUTOT=(ITTOT*(ITTOT-1))/2+IUTOT
          FTU=FIFA(IOFF(ISU)+ITUTOT)
          IF(ABS(FTU).LT.1.0D-16) GOTO 20
          CALL SIGMA1(SGS,CIS,EXS,
     &                LEVT,LEVU,FTU,STSYM,BRA,SGM)
  20      CONTINUE
        END DO
      END DO

      IF (IFTEST.GT.0) THEN
        WRITE(6,*)' SGM array from (strictly lowering F)|BRA>:'
        WRITE(6,'(1x,5F16.8)')(SGM(I),I=1,NCONF)
      END IF

* Load ket wave function
      ID=IDCIEX
      DO I=1,IKET-1
        CALL DDAFILE(LUCIEX,0,KET,NCONF,ID)
      END DO
      CALL DDAFILE(LUCIEX,2,KET,NCONF,ID)

* Add contribution to matrix element FOPEL
      FOPEL=FOPEL+DDOT_(NCONF,KET,1,SGM,1)

      IF (IFTEST.GT.0) THEN
        WRITE(6,*)' FOPEL is now:'
        WRITE(6,'(1x,5f16.8)')(FOPEL)
      END IF

      CALL mma_deallocate(SGM)
      CALL mma_deallocate(BRA)
      CALL mma_deallocate(KET)

      END SUBROUTINE FOPAB

