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
* Copyright (C) 1992,1994, Per Ake Malmqvist                           *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE POLY1(CI)
      IMPLICIT NONE
* PER-AAKE MALMQUIST, 92-12-07
* THIS PROGRAM CALCULATES THE 1-EL DENSITY
* MATRIX FOR A CASSCF WAVE FUNCTION.
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "pt2_guga.fh"
#include "SysDef.fh"

      REAL*8, INTENT(IN) :: CI(NCONF)

      INTEGER LSGM1,LG1TMP

      INTEGER I
      REAL*8, EXTERNAL :: DDOT_,DNRM2_

      CALL QENTER('POLY1')

      IF(NLEV.GT.0) THEN
        CALL GETMEM('LSGM1','ALLO','REAL',LSGM1 ,MXCI)
        CALL GETMEM('LG1TMP','ALLO','REAL',LG1TMP,NG1)
        CALL DENS1_RPT2(CI,WORK(LSGM1),WORK(LG1TMP))
      END IF

* REINITIALIZE USE OF DMAT.
* The fields IADR10 and CLAB10 are kept in common included from pt2_guga.fh
* CLAB10 replaces older field called LABEL.
      DO I=1,64
        IADR10(I,1)=-1
        IADR10(I,2)=0
        CLAB10(I)='   EMPTY'
      END DO
      IADR10(1,1)=0
* HENCEFORTH, THE CALL PUT(NSIZE,LABEL,ARRAY) WILL ENTER AN
* ARRAY ON LUDMAT AND UPDATE THE TOC.
      IF(NLEV.GT.0) THEN
        CALL PT2_PUT(NG1,' GAMMA1',WORK(LG1TMP))

        CALL GETMEM('LSGM1','FREE','REAL',LSGM1 ,MXCI)
        CALL GETMEM('LG1TMP','FREE','REAL',LG1TMP,NG1)
      END IF

      CALL QEXIT('POLY1')

      RETURN
      END

