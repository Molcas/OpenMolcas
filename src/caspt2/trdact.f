************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1998, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1998  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE TRDACT(IVEC,JVEC,DTU)
      IMPLICIT NONE

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "pt2_guga.fh"
#include "stdalloc.fh"
#include "SysDef.fh"

      INTEGER IVEC, JVEC
      REAL*8 DTU(NASHT,NASHT)
C Local array:
      INTEGER IATOG(MXLEV)

      REAL*8, ALLOCATABLE :: TRDOP1(:), TRDOP2(:), TRDOP3(:)
      INTEGER :: NOP1, NOP2, NOP3

      REAL*8, ALLOCATABLE :: TRDTMP(:), TRDCI(:), TRDSGM(:)
      INTEGER :: NTMP, NSGM

      INTEGER :: I, J, ID
      INTEGER :: ISYM, ISYMT
      INTEGER :: ITABS, ITLEV, IU, IUABS, IULEV
      REAL*8 :: OP0, OCCNUM, SCP, DDOT_

C Add to the active-active block of transition density matrix,
C    D(t,u) = Add <IVEC| E(t,u) |JVEC> = <0| W1T E(t,u) W2 |0>
C where t,u are active indices. IVEC and JVEC are integer labels,
C denoting sets of coefficients stored on LUSOLV. These are assumed
C to be contravariant representations of the wave operators W1 and W2,
C in the notation of the comments.

      CALL QENTER('TRDACT')

CTEST      WRITE(*,*)' Memory list, TRDACT point A:'
CTEST      call getmem('TRDACT A','LIST','REAL',LDUM,NDUM)

C (1): Compute a representation of the operator PCAS*W1T*W2
      NOP1=NASHT**2
      NOP2=(NOP1*(NOP1+1))/2
      NOP3=(NOP2*(NOP1+2))/3
      CALL MMA_ALLOCATE(TRDOP1,NOP1)
      CALL MMA_ALLOCATE(TRDOP2,NOP2)
      CALL MMA_ALLOCATE(TRDOP3,NOP3)
      CALL MKWWOP(IVEC,JVEC,OP0,TRDOP1,NOP2,TRDOP2,NOP3,TRDOP3)
CTEST      WRITE(*,*)' TRDACT, back from MKWWOP.'
CTEST      WRITE(*,*)' OP0:'
CTEST      WRITE(*,'(1x,5f16.8)') OP0
CTEST      WRITE(*,*)' OP1:'
CTEST      WRITE(*,'(1x,5f16.8)')(TRDOP1(I),I=1,NOP1)
CTEST      WRITE(*,*)' OP2:'
CTEST      WRITE(*,'(1x,5f16.8)')(TRDOP2(I),I=1,NOP2)
CTEST      WRITE(*,*)' OP3:'
CTEST      WRITE(*,'(1x,5f16.8)')(TRDOP3(I),I=1,NOP3)

CTEST      WRITE(*,*)' Memory check, TRDACT point B:'
CTEST      call getmem('TRDACT B','CHEC','REAL',LDUM,NDUM)

C (2): Compute the state vector |Temp> = (PCAS*W1T*W2) |0>
C First modify the coefficients, see subroutine MODOP.
      CALL MODOP(TRDOP1,NOP2,TRDOP2,NOP3,TRDOP3)
CTEST      WRITE(*,*)' TRDACT, back from MODOP.'
CTEST      WRITE(*,*)' OP0:'
CTEST      WRITE(*,'(1x,5f16.8)') OP0
CTEST      WRITE(*,*)' OP1:'
CTEST      WRITE(*,'(1x,5f16.8)')(TRDOP1(I),I=1,NOP1)
CTEST      WRITE(*,*)' OP2:'
CTEST      WRITE(*,'(1x,5f16.8)')(TRDOP2(I),I=1,NOP2)
CTEST      WRITE(*,*)' OP3:'
CTEST      WRITE(*,'(1x,5f16.8)')(TRDOP3(I),I=1,NOP3)
      NTMP=NCONF
      CALL MMA_ALLOCATE(TRDTMP,NTMP)
      CALL MMA_ALLOCATE(TRDCI,NCONF)
      IF(ISCF.EQ.0) THEN
*PAM07 Eliminate the unsafe IPOSFILE call
*        ID=IDTCEX+iPosFile(NCONF)*(JSTATE-1)
*PAM07 Use instead dummy operations:
        ID=IDTCEX
        DO J=1,JSTATE-1
         CALL DDAFILE(LUCIEX,0,TRDCI,NCONF,ID)
        END DO
        CALL DDAFILE(LUCIEX,2,TRDCI,NCONF,ID)
CTEST      WRITE(*,*)' Before HAM3, the CI array:'
CTEST      CALL PRWF_CP2(LSYM,NCONF,TRDCI,0.0D0)
CTEST      WRITE(*,*)' Norm:',DNRM2_(NCONF,TRDCI,1)
      ELSE
        TRDCI=1.0D0
      END IF
      CALL DCOPY_(NTMP,[0.0D0],0,TRDTMP,1)
      CALL HAM3(OP0,TRDOP1,NOP2,TRDOP2,NOP3,TRDOP3,LSYM,TRDCI,TRDTMP)
CTEST      WRITE(*,*)' After HAM3, the wave function W(+)*W*(PSI0):'
CTEST      CALL PRWF_CP2(LSYM,NCONF,TRDTMP,0.0D0)
CTEST      WRITE(*,*)' Ovlp with CI:',DDOT_(NCONF,TRDCI,1,TRDTMP,1)
C No more need for the operators:
      CALL MMA_DEALLOCATE(TRDOP1)
      CALL MMA_DEALLOCATE(TRDOP2)
      CALL MMA_DEALLOCATE(TRDOP3)

CTEST      WRITE(*,*)' Memory check, TRDACT point C:'
CTEST      call getmem('TRDACT C','CHEC','REAL',LDUM,NDUM)

C (3) compute <0| E(t,u) W1T W2 |0> as <Sigma_ut| Temp>, and add to DTU

      IF(ISCF.EQ.0) THEN
C Create reorder table giving the GUGA level, i.e. CI-coupling
C ordinal number of each active orbital.
        ITABS=0
        DO ISYM=1,NSYM
          DO I=1,NLEV
            IF(ISM(I).EQ.ISYM) THEN
              ITABS=ITABS+1
              IATOG(ITABS)=I
            END IF
          END DO
        END DO
CTEST      WRITE(*,*)' At point 3 in TRDACT. The scalar products:'
        NSGM=NCONF
        CALL MMA_ALLOCATE(TRDSGM,MXCI)
        DO ITABS=1,NASHT
          ISYMT=IASYM(ITABS)
          ITLEV=IATOG(ITABS)
          DO IU=1,NASH(ISYMT)
            IUABS=NAES(ISYMT)+IU
            IULEV=IATOG(IUABS)
CPAM00          CALL GETSGM(IULEV,ITLEV,IDEX,TRDSGM)
CPAM00 GETSGM replaced by GETSGM2
            CALL GETSGM2(IULEV,ITLEV,LSYM,TRDCI,TRDSGM)
            SCP=DDOT_(NCONF,TRDSGM,1,TRDTMP,1)
CTEST      WRITE(*,*)' ITABS,IUABS,SCP:',ITABS,IUABS,SCP
            DTU(ITABS,IUABS)=DTU(ITABS,IUABS)+SCP
          END DO
        END DO
        CALL MMA_DEALLOCATE(TRDSGM)
      ELSE
        OCCNUM=2.0D0
        IF(ISCF.EQ.2) OCCNUM=1.0D0
        DO ITABS=1,NASHT
          DTU(ITABS,ITABS)=DTU(ITABS,ITABS)+OCCNUM
        END DO
      END IF
CPAM00 No more need for CI array
      CALL MMA_DEALLOCATE(TRDCI)
C No more need for the TMP state vector
      CALL MMA_DEALLOCATE(TRDTMP)

CTEST      WRITE(*,*)' After point 3 in TRDACT, the array DTU is:'
CTEST      do i=1,nasht
CTEST        WRITE(*,'(1x,6f12.6)')(dtu(i,j),j=1,nasht)
CTEST      end do

CTEST      WRITE(*,*)' Memory check, TRDACT point D:'
CTEST      call getmem('TRDACT D','CHEC','REAL',LDUM,NDUM)

C (4): Add the correction <0| [W1T,E(tu)] W2 |0>.
      CALL COMMWEW(IVEC,JVEC,DTU)

CTEST      WRITE(*,*)' After point 4 in TRDACT, the array DTU is:'
CTEST      do i=1,nasht
CTEST        WRITE(*,'(1x,6f12.6)')(dtu(i,j),j=1,nasht)
CTEST      end do

CTEST      WRITE(*,*)' Memory check, TRDACT point E:'
CTEST      call getmem('TRDACT E','CHEC','REAL',LDUM,NDUM)

      CALL QEXIT('TRDACT')
      RETURN
      END
