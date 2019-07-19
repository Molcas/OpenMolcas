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
      SUBROUTINE MKFOP(FIFA,NGRP,JSTATE_OFF,FOPXMS)
      IMPLICIT REAL*8 (a-h,o-z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "SysDef.fh"
#include "WrkSpc.fh"
#include "pt2_guga.fh"
      DIMENSION FIFA(NFIFA)
      DIMENSION FOPXMS(NGRP,NGRP)
* Purely local array, offsets:
      DIMENSION IOFF(8)

* Procedure for computing the matrix FOPXMS, which is the matrix
* of the average Fock operator in the basis of input CASSCF w.f.
* In: The average Fock matrix, active indices only, over the original
* CASSCF MO-orbitals.
* Compute the matrix elements of this Fock operator over the basis
* of CASSCF states that are treated together in XMS fashion.

      CALL QENTER('MKFOP')

      IFTEST=0

* Offset table for accessing FIFA array:
      IOF=0
      DO ISYM=1,NSYM
       IOFF(ISYM)=IOF
       IOF=IOF+(NORB(ISYM)*(NORB(ISYM)+1))/2
      END DO

      IF (IFTEST.gt.0) THEN
       write(6,*)' The FIFA array:'
       DO ISYM=1,NSYM
        IJ=IOFF(ISYM)+1
        DO I=1,NORB(ISYM)
         WRITE(6,'(1x,5F16.8)')(FIFA(IJ+J),J=0,I-1)
         IJ=IJ+I
        END DO
       END DO
       CALL XFLUSH(6)
      END IF

      IF (ISCF.eq.1 .or. ISCF.eq.2) THEN
* Specialized code for Closed-shell or Hi-spin HF:
* Sum up diagonal elements of FIFA times occ. number
       ESUM=0.0D0
       DO ISYM=1,NSYM
        OCC=2.0D0
        DO I=1,NISH(ISYM)
         ESUM=ESUM+OCC*FIFA(IOFF(ISYM)+(I*(I+1))/2)
        END DO
        IF (ISCF.eq.2) OCC=1.0D0
        DO J=1,NASH(ISYM)
         I=NISH(ISYM)+J
         ESUM=ESUM+OCC*FIFA(IOFF(ISYM)+(I*(I+1))/2)
        END DO
       END DO
       FOPXMS(1,1)=ESUM
       CALL QEXIT('MKFOP')
       RETURN
      END IF

* General CASSCF or RASSCF case:
* Sum up trace of FIFA over inactive orbitals only:
      TRC=0.0D0
      DO ISYM=1,NSYM
       DO I=1,NISH(ISYM)
        II=IOFF(ISYM)+(I*(I+1))/2
        TRC=TRC+FIFA(II)
       END DO
      END DO
* Contribution from inactive orbitals:
      EINACT=2.0D0*TRC
      IF (IFTEST.gt.0) THEN
       write(6,*)' Energy contrib from inactive orbitals:',EINACT
       CALL XFLUSH(6)
      END IF

* Allocate arrays for ket and bra wave functions
      CALL GETMEM('LBRA','ALLO','REAL',LBRA,NCONF)
      CALL GETMEM('LKET','ALLO','REAL',LKET,NCONF)
* Allocate array for sigma, = average Fock operator acting on ket:
      CALL GETMEM('SGM','ALLO','REAL',LSGM,NCONF)
      CALL DCOPY_(NGRP**2,[0.0D0],0,FOPXMS,1)

* Loop over ket functions:
      DO J2=1,NGRP
       IKET=JSTATE_OFF+J2
       ID=IDCIEX
       DO I=1,IKET-1
        CALL DDAFILE(LUCIEX,0,WORK(LKET),NCONF,ID)
       END DO
       CALL DDAFILE(LUCIEX,2,WORK(LKET),NCONF,ID)

      IF (IFTEST.gt.0) THEN
       write(6,*)' IKET:',IKET
       write(6,*)' Ket CI array:'
       ISCR=MIN(NCONF,20)
       write(6,'(1x,5F16.8)')(WORK(LKET+I),I=0,ISCR-1)
       call xflush(6)
      END IF

* Compute (lowering part of) FIFA operator acting on
* the ket wave function.
       CALL DCOPY_(NCONF,[0.0D0],0,WORK(LSGM),1)
       DO LEVU=1,NLEV
        IUABS=L2ACT(LEVU)
        ISU=ISM(LEVU)
        IU=IUABS-NAES(ISU)
        NI=NISH(ISU)
        IUTOT=NI+IU
        DO LEVT= 1,LEVU
         IF(ISM(LEVT).NE.ISU) GOTO 10
         ITABS=L2ACT(LEVT)
         IST=ISU
         IT=ITABS-NAES(IST)
         ITTOT=NI+IT
         ITUTOT=(IUTOT*(IUTOT-1))/2+ITTOT
         IF (ITTOT.GT.IUTOT) ITUTOT=(ITTOT*(ITTOT-1))/2+IUTOT
         FTU=FIFA(IOFF(ISU)+ITUTOT)
         IF(ABS(FTU).lt.1.0D-16) GOTO 10
         CALL SIGMA1_CP2(LEVT,LEVU,FTU,LSYM,WORK(LKET),WORK(LSGM),
     &     IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &     IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &     WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
  10     CONTINUE
        END DO
       END DO
* Add contribution from inactive part:
       CALL DAXPY_(NCONF,EINACT,WORK(LKET),1,WORK(LSGM),1)

      IF (IFTEST.gt.0) THEN
       write(6,*)' SGM array from (lowering F)|KET>:'
       write(6,'(1x,5F16.8)')(WORK(LSGM+I),I=0,NCONF-1)
       call xflush(6)
      END IF

* Loop over all bra functions:
       DO J1=1,NGRP
        IBRA=JSTATE_OFF+J1
        ID=IDCIEX
        DO I=1,IBRA-1
         CALL DDAFILE(LUCIEX,0,WORK(LBRA),NCONF,ID)
        END DO
        CALL DDAFILE(LUCIEX,2,WORK(LBRA),NCONF,ID)
* Put matrix element into FOPXMS:
        FOPXMS(J1,J2)=DDOT_(NCONF,WORK(LBRA),1,WORK(LSGM),1)
       END DO
      END DO

      IF (IFTEST.gt.0) THEN
       write(6,*)' FOPXMS is now:'
       do J1=1,NGRP
        write(6,'(1x,5f16.8)')(FOPXMS(J1,J2),J2=1,NGRP)
       end do
       call xflush(6)
      END IF

* Repeat, now switch role of bra/ket:
* Loop over bra functions:
      DO J1=1,NGRP
       IBRA=JSTATE_OFF+J1
       ID=IDCIEX
       DO I=1,IBRA-1
        CALL DDAFILE(LUCIEX,0,WORK(LBRA),NCONF,ID)
       END DO
       CALL DDAFILE(LUCIEX,2,WORK(LBRA),NCONF,ID)

      IF (IFTEST.gt.0) THEN
       write(6,*)' IBRA:',IBRA
       write(6,*)' Bra CI array:'
       write(6,'(1x,5F16.8)')(WORK(LBRA+I),I=0,NCONF-1)
       call xflush(6)
      END IF

* Compute (strictly lowering part of) FIFA operator acting on |BRA>.
*  We are computing contributions <KET|Etu|BRA> with t<u, then
* using them as <BRA|Eut|KET>
       CALL DCOPY_(NCONF,[0.0D0],0,WORK(LSGM),1)
       DO LEVU=2,NLEV
        IUABS=L2ACT(LEVU)
        ISU=ISM(LEVU)
        IU=IUABS-NAES(ISU)
        NI=NISH(ISU)
        IUTOT=NI+IU
        DO LEVT= 1,LEVU-1
         IF(ISM(LEVT).NE.ISU) GOTO 20
         ITABS=L2ACT(LEVT)
         IST=ISU
         IT=ITABS-NAES(IST)
         ITTOT=NI+IT
         ITUTOT=(IUTOT*(IUTOT-1))/2+ITTOT
         IF (ITTOT.GT.IUTOT) ITUTOT=(ITTOT*(ITTOT-1))/2+IUTOT
         FTU=FIFA(IOFF(ISU)+ITUTOT)
         IF(ABS(FTU).lt.1.0D-16) GOTO 20
         CALL SIGMA1_CP2(LEVT,LEVU,FTU,LSYM,WORK(LBRA),WORK(LSGM),
     &     IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &     IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &     WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
  20     CONTINUE
        END DO
       END DO

      IF (IFTEST.gt.0) THEN
       write(6,*)' SGM array from (strictly lowering F)|BRA>:'
       write(6,'(1x,5F16.8)')(WORK(LSGM+I),I=0,NCONF-1)
       call xflush(6)
      END IF

* Loop over all ket functions:
       DO J2=1,NGRP
        IKET=JSTATE_OFF+J2
        ID=IDCIEX
        DO I=1,IKET-1
         CALL DDAFILE(LUCIEX,0,WORK(LKET),NCONF,ID)
        END DO
        CALL DDAFILE(LUCIEX,2,WORK(LKET),NCONF,ID)
* Put matrix element into FOPXMS:
        FOPXMS(J1,J2)=FOPXMS(J1,J2)+
     &         DDOT_(NCONF,WORK(LKET),1,WORK(LSGM),1)
       END DO
      END DO

      IF (IFTEST.gt.0) THEN
       write(6,*)' FOPXMS is now:'
       do J1=1,NGRP
        write(6,'(1x,5f16.8)')(FOPXMS(J1,J2),J2=1,NGRP)
       end do
       call xflush(6)
      END IF

      CALL GETMEM('SGM','FREE','REAL',LSGM,NCONF)
      CALL GETMEM('LBRA','FREE','REAL',LBRA,NCONF)
      CALL GETMEM('LKET','FREE','REAL',LKET,NCONF)
      CALL QEXIT('MKFOP')
      RETURN
      END
