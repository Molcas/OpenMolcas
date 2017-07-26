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
* Copyright (C) 1989,1998, Per Ake Malmqvist                           *
************************************************************************
*  PROGRAM RASSI        PER-AAKE MALMQVIST
*  SUBROUTINE CITRA     IBM-3090 RELEASE 89 01 31
*  USE THE COEFFICIENTS FOR A SEQUENCE OF SINGLE-ORBITAL TRANSFOR-
*  MATION, TRA, TO TRANSFORM THE CI EXPANSION COEFFICIENTS
*  IN-PLACE TO A NEW NON-ON ORBITAL BASIS.
*  NEW VERSION 981122, using arrays ISGS,ICIS,IXS.
************************************************************************
*  CITRA
*
*> @brief
*>   Recompute a CI coefficient array to use another orbital basis
*> @author P. &Aring;. Malmqvist
*>
*> @details
*> For a given linear transformation of the orbitals, and a
*> CI array where the CSF basis is built from the old
*> orbitals, compute the CI array using the new orbitals
*> instead. The orbitals are transformed sequentially, and
*> for each active orbital, a call to ::SSOTRA performs the
*> single-orbital transformation.
*>
*> @param[in]     WFTP Wave function Type Name
*> @param[in]     ISGS Split Graph Structure Array
*> @param[in]     ICIS CI Structure Array
*> @param[in]     IXS  Excitation operator Structure Array
*> @param[in]     LSM  Wave function Symmetry Label
*> @param[in]     TRA  Transformation Matrix
*> @param[in]     NCO  Number of Configuration Functions
*> @param[in,out] CI   CI Array
************************************************************************
      SUBROUTINE CITRA(WFTP,ISGS,ICIS,IXS,LSM,TRA,NCO,CI)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='CITRA')
      DIMENSION TRA(NTRA),CI(NCO)
#include "WrkSpc.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Struct.fh"
      CHARACTER*8 WFTP
      DIMENSION ISGS(NSGSIZE),ICIS(NCISIZE),IXS(NXSIZE)

      CALL QENTER(ROUTINE)

CTEST      write(*,*)' Entering CITRA. TRA='
CTEST      write(*,'(1x,5f16.8)')(TRA(I),I=1,NTRA)
C TRA contains square matrices, one per symmetry
C  FIRST TRANSFORM THE INACTIVE ORBITALS:
      FAC=1.0D00
      ISTA=1
      DO ISYM=1,NSYM
        NO=NOSH(ISYM)
        DO I=1,NISH(ISYM)
          II=ISTA+(NO+1)*(I-1)
          CKK=TRA(II)
          FAC=FAC*CKK
        END DO
        ISTA=ISTA+NO**2
      END DO
      FAC=FAC**2
      CALL DSCAL_(NCO,FAC,CI,1)
C  THEN THE ACTIVE ONES:
      IF(WFTP.EQ.'EMPTY   ') GOTO 100
* The HISPIN case may be buggy and is not presently used.
      IF(WFTP.EQ.'HISPIN  '.or.WFTP.EQ.'CLOSED  ') THEN
        ISTA=1
        DO ISYM=1,NSYM
          NI=NISH(ISYM)
          NA=NASH(ISYM)
          NO=NOSH(ISYM)
          DO I=NI+1,NI+NA
            II=ISTA+(NO+1)*(I-1)
            CKK=TRA(II)
            FAC=FAC*CKK
          END DO
          ISTA=ISTA+NO**2
        END DO
        IF(WFTP.EQ.'CLOSED  ') FAC=FAC**2
        CALL DSCAL_(NCO,FAC,CI,1)
      ELSE
C The general case:
        CALL GETMEM('TMP   ','ALLO','REAL',LTMP,NCO)
        ISTA=1
        DO ISYM=1,NSYM
          NA=NASH(ISYM)
          NO=NOSH(ISYM)
          IF(NA.NE.0) THEN
            CALL SSOTRA(ISGS,ICIS,IXS,ISYM,LSM,NA,NO,
     *                TRA(ISTA),NCO,CI,WORK(LTMP))
          END IF
          ISTA=ISTA+NO**2
        END DO
        CALL GETMEM('TMP   ','FREE','REAL',LTMP,NCO)
      END IF

 100  CONTINUE
      CALL QEXIT(ROUTINE)
      RETURN
      END
