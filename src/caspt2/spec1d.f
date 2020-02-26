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
      SUBROUTINE SPEC1D(IFC,FACT,X,Y)
      USE SUPERINDEX
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
      INTEGER IFC
      REAL*8 FACT,X(*),Y(*)
      INTEGER NAS,NIS,ITQ,ITT
C If IFC=0, compute
C X(tt1,ia) <- X(tt1,ia)+FACT*Y(ia), else
C the conjugate expression (summing into Y, values from X).

      NAS=NASUP(1,5)
      NIS=NISUP(1,5)
      IF(NIS*NAS.EQ.0) RETURN
      IF(IFC.EQ.0) THEN
        DO 20 ITQ=1,NASHT
          ITT=KTU(ITQ,ITQ)
          CALL DAXPY_(NIS,FACT,Y,1,X(ITT),NAS)
  20    CONTINUE
      ELSE
        CALL DCOPY_(NIS,[0.0D0],0,Y,1)
        DO 30 ITQ=1,NASHT
          ITT=KTU(ITQ,ITQ)
          CALL DAXPY_(NIS,FACT,X(ITT),NAS,Y,1)
  30    CONTINUE
      END IF
      RETURN
      END
