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
      SUBROUTINE SPEC1A(IFC,FACT,ISYM,X,nX,Y,nY)
      USE SUPERINDEX, only: KTUV
      use caspt2_module, only: nTUV, nAsh, nIsh, nAES, nTUVES, nAshT
      use definitions, only: iwp, wp
      IMPLICIT NONE
      INTEGER(kind=iwp), intent(in):: IFC,ISYM, nX, nY
      REAL(kind=wp), intent(in):: FACT
      REAL(kind=wp), intent(inout):: X(nX),Y(nY)

      INTEGER(kind=iwp) NAS,NT,NI,IT,ITQ,IUQ,ITUU
C If IFC=0, compute
C X(tuu,i) <- X(tuu,i)+FACT*Y(t,i), else
C the conjugate expression (summing into Y, values from X).
      NI=NISH(ISYM)
      IF (NI<=0) RETURN
      NAS=NTUV(ISYM)
      NT=NASH(ISYM)
      DO IT=1,NT
        ITQ=IT+NAES(ISYM)
        DO IUQ=1,NASHT
          ITUU=KTUV(ITQ,IUQ,IUQ)-NTUVES(ISYM)
          IF(IFC==0) THEN
            CALL DAXPY_(NI,FACT,Y(IT),NT,X(ITUU),NAS)
          ELSE
            CALL DAXPY_(NI,FACT,X(ITUU),NAS,Y(IT),NT)
          END IF
        END DO
      END DO
      END SUBROUTINE SPEC1A
