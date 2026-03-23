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
* Copyright (C) 2006, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 2006  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE MKDREF_RPT2(N,G1,DREF,nDREF)
      use definitions, only: iwp, wp
      IMPLICIT NONE
      INTEGER(kind=iwp), INTENT(IN) :: N, NDREF
      REAL(kind=wp), INTENT(IN) :: G1(N,N)
      REAL(kind=wp), INTENT(OUT) :: DREF(NDREF)

      INTEGER(kind=iwp) I,J,IJ

C Compute DREF(PQ) = <0| Epq |0>
C from G1(P,Q) = <0| Epq |0>
C Storage differs: DREF is triangular.

      DO I=1,N
       DO J=1,I
        IJ=(I*(I-1))/2+J
        DREF(IJ)=G1(I,J)
       END DO
      END DO

      END SUBROUTINE MKDREF_RPT2
