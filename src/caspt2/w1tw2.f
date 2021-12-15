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
      SUBROUTINE W1TW2(IVEC,JVEC,CI,SGM)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"

#include "SysDef.fh"
      dimension ci(*),sgm(*)

C Given contravariant indices of two wave operators W1 and W2,
C in the vectors numbered IVEC and JVEC on file (unit LUSOLV),
C compute the vector in CAS space
C   | SGM > := | SGM > + (W1 conj)*(W2)*| CI >


C (1): Compute a representation of the operator PCAS*W1T*W2
      NOP1=NASHT**2
      NOP2=(NOP1*(NOP1+1))/2
      NOP3=(NOP2*(NOP1+2))/3
      CALL GETMEM('TRDOP1','ALLO','REAL',LOP1,NOP1)
      CALL GETMEM('TRDOP2','ALLO','REAL',LOP2,NOP2)
      CALL GETMEM('TRDOP3','ALLO','REAL',LOP3,NOP3)
      CALL MKWWOP(IVEC,JVEC,OP0,WORK(LOP1),NOP2,WORK(LOP2),
     &                                           NOP3,WORK(LOP3))
C Modify the coefficients, see subroutine MODOP.
      CALL MODOP(WORK(LOP1),NOP2,WORK(LOP2),NOP3,WORK(LOP3))

CTEST      WRITE(*,*)' W1TW2 after MODOP.'
CTEST      WRITE(*,*)' OP0:'
CTEST      WRITE(*,'(1x,5f16.8)') OP0
CTEST      WRITE(*,*)' OP1:'
CTEST      WRITE(*,'(1x,5f16.8)')(WORK(LOP1-1+I),I=1,NOP1)
CTEST      WRITE(*,*)' OP2:'
CTEST      WRITE(*,'(1x,5f16.8)')(WORK(LOP2-1+I),I=1,NOP2)
CTEST      WRITE(*,*)' OP3:'
CTEST      WRITE(*,'(1x,5f16.8)')(WORK(LOP3-1+I),I=1,NOP3)

C (2) Apply the operators:
CTEST      WRITE(*,*)' In W1TW2, before HAM3, the CI array:'
CTEST      CALL PRWF_CP2(LSYM,NCONF,CI,0.05D0)
      CALL HAM3(OP0,WORK(LOP1),NOP2,WORK(LOP2),NOP3,WORK(LOP3),
     &                           LSYM,CI,SGM)
CTEST      WRITE(*,*)' In W1TW2,  after HAM3, the SGM array:'
CTEST      CALL PRWF_CP2(LSYM,NCONF,SGM,0.05D0)

      CALL GETMEM('TRDOP1','FREE','REAL',LOP1,NOP1)
      CALL GETMEM('TRDOP2','FREE','REAL',LOP2,NOP2)
      CALL GETMEM('TRDOP3','FREE','REAL',LOP3,NOP3)
      RETURN
      END
