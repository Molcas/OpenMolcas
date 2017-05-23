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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE SPEC1C(IFC,FACT,ISYM,X,Y)
      USE SUPERINDEX
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
      INTEGER IFC,ISYM
      REAL*8 FACT,X(*),Y(*)
      INTEGER NAS,NT,NA,IT,ITQ,IUQ,ITUU
C If IFC=0, compute
C X(tuu,a) <- X(tuu,a)+FACT*Y(t,a), else
C the conjugate expression (summing into Y, values from X).

      NAS=NTUV(ISYM)
      NT=NASH(ISYM)
      NA=NSSH(ISYM)
      DO 10 IT=1,NT
        ITQ=IT+NAES(ISYM)
        DO 10 IUQ=1,NASHT
          ITUU=KTUV(ITQ,IUQ,IUQ)-NTUVES(ISYM)
          IF(IFC.EQ.0) THEN
            CALL DAXPY_(NA,FACT,Y(IT),NT,X(ITUU),NAS)
          ELSE
            CALL DAXPY_(NA,FACT,X(ITUU),NAS,Y(IT),NT)
          END IF
  10  CONTINUE
      RETURN
      END
