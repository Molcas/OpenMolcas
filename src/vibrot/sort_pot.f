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
* Copyright (C) 1981, Per-Olof Widmark                                 *
************************************************************************
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  THIS ROUTINE SORTS THE VECTORS X AND Y, WHERE THE VECTOR X IS       C
C  USED AS THE KEYS.                                                   C
C                                                                      C
C  DATE: 81 04 14                                                      C
C                                                                      C
C  AUTHOR: P-O WIDMARK                                                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SORT_POT(X,Y,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(N),Y(N)
      DO 10 I=2,N
         IND=I-1
         SMALL=X(IND)
         DO 20 J=I,N
            IF(X(J).GT.SMALL) GOTO 20
            IND=J
            SMALL=X(IND)
20       CONTINUE
         Z=X(I-1)
         X(I-1)=X(IND)
         X(IND)=Z
         Z=Y(I-1)
         Y(I-1)=Y(IND)
         Y(IND)=Z
10    CONTINUE
      RETURN
      END
