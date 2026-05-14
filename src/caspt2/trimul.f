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
      SUBROUTINE TRIMUL(N,M,ALPHA,ASYM,X,LDX,Y,LDY)
      use constants, only: One
      use definitions, only: iwp, wp
      IMPLICIT NONE
      INTEGER(kind=iwp), intent(in):: N,M,LDX,LDY
      REAL(kind=wp), intent(in):: ALPHA,ASYM((N*(N+1))/2),X(LDX,M)
      REAL(kind=wp), intent(inout):: Y(LDY,M)

      INTEGER(kind=iwp) I
C Multiply symmetric matrix ASYM with matrix X.
C Scale result with ALPHA and add it to matrix Y.
      DO I=1,M
!       CALL DSLMX(N,ALPHA,ASYM,X(1,I),1,Y(1,I),1)
        CALL DSPMV_('U',N,ALPHA,ASYM,X(1,I),1,One,Y(1,I),1)
      End Do

      END SUBROUTINE TRIMUL
