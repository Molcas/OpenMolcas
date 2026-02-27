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
      SUBROUTINE SPEC1D(IFC,FACT,X,nX,Y,nY)
      USE SUPERINDEX, only: KTU
      USE caspt2_module, only: nAshT, nASup, nISup
      use constants, only: Zero
      use definitions, only: iwp, wp
      IMPLICIT NONE
      INTEGER(kind=iwp), intent(in):: IFC, nX, nY
      REAL(kind=wp), intent(in):: FACT
      REAL(kind=wp), intent(inout):: X(nX),Y(nY)

      INTEGER(kind=iwp) NAS,NIS,ITQ,ITT
C If IFC=0, compute
C X(tt1,ia) <- X(tt1,ia)+FACT*Y(ia), else
C the conjugate expression (summing into Y, values from X).

      NIS=NISUP(1,5)
      IF(NIS==0) RETURN
      NAS=NASUP(1,5)
      IF(IFC==0) THEN
        DO ITQ=1,NASHT
          ITT=KTU(ITQ,ITQ)
          CALL DAXPY_(NIS,FACT,Y,1,X(ITT),NAS)
        END DO
      ELSE
        Y(1:NIS)=Zero
        DO ITQ=1,NASHT
          ITT=KTU(ITQ,ITQ)
          CALL DAXPY_(NIS,FACT,X(ITT),NAS,Y,1)
        END DO
      END IF
      END SUBROUTINE SPEC1D
