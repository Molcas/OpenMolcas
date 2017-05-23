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
      SUBROUTINE SSOTRA(ISGS,ICIS,IXS,ISYM,LSM,NA,NO,TRA,NCO,CI,TMP)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='SSOTRA')
      DIMENSION TRA(NO,NO),CI(NCO),TMP(NCO)
#include "rassi.fh"
#include "Struct.fh"
#include "WrkSpc.fh"
      Dimension ISGS(nSGSize)
      Dimension ICIS(nCISize)
      Dimension IXS (nXSize)



      CALL QENTER(ROUTINE)

C Dereference ISGS to get at ISM table:
      LISM=ISGS(3)
C ILEV(IORB)=GUGA LEVEL CORRESPONDING TO A SPECIFIC ACTIVE ORBITAL
C OF SYMMETRY ISYM.
      CALL GETMEM('ILEV','ALLO','INTE',LILEV,NA)
      NI=NO-NA
      IL=0
      DO 10 IP=1,NA
5       IL=IL+1
        IF(IWORK(LISM-1+IL).NE.ISYM) GOTO 5
        IWORK(LILEV-1+IP)=IL
10    CONTINUE
CTEST      write(*,*)' Check prints in SSOTRA.'
CTEST      write(*,*)' ISYM:',ISYM
      DO 100 IK=1,NA
        IKLEV=IWORK(LILEV-1+IK)
        CALL DCOPY_(NCO,0.0D0,0,TMP,1)
        DO 50 IP=1,NA
          IPLEV=IWORK(LILEV-1+IP)
          CPK=TRA(NI+IP,NI+IK)
          IF(IP.EQ.IK) CPK=CPK-1.0D00
          X=0.5D0*CPK
CTEST          write(*,*)' IP,IK,X:',IP,IK,X
          IF(ABS(X).LT.1.0D-14) GOTO 50
CPAM98          CALL SIGMA_1(IPLEV,IKLEV,X,LSM,CI,TMP,IWORK(LNOCSF),
CPAM98     *                 IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
CPAM98     *                 IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
CPAM98     *                 WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
CPAM98 Replaced by more modern routine:
          CALL SGMONE(ISGS,ICIS,IXS,IPLEV,IKLEV,X,LSM,CI,TMP)
50      CONTINUE
        CKK=TRA(NI+IK,NI+IK)
        X= 3.0D00-CKK
        CALL DAXPY_(NCO,X,TMP,1,CI,1)
        DO 60 IP=1,NA
          IPLEV=IWORK(LILEV-1+IP)
          CPK=TRA(NI+IP,NI+IK)
          IF(IP.EQ.IK) CPK=CPK-1.0D00
          IF(ABS(CPK).LT.1.0D-14) GOTO 60
CPAM98          CALL SIGMA_1(IPLEV,IKLEV,CPK,LSM,TMP,CI,IWORK(LNOCSF),
CPAM98     *                 IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
CPAM98     *                 IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
CPAM98     *                 WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          CALL SGMONE(ISGS,ICIS,IXS,IPLEV,IKLEV,CPK,LSM,TMP,CI)

60      CONTINUE
100   CONTINUE
      CALL GETMEM('ILEV','FREE','INTE',LILEV,NA)

      CALL QEXIT(ROUTINE)
      RETURN
      END
