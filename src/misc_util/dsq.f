************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE DSQ(A,B,ICB,IRB,NROW)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),B(*)
      IND=0
      DO 10 IROW=0,NROW-1
      DO 10 ICOL=0,IROW
         IND=IND+1
         B(1+IROW*ICB+ICOL*IRB)=0.5D0*A(IND)
         B(1+ICOL*ICB+IROW*IRB)=0.5D0*A(IND)
10    CONTINUE
      DO 20 IROW=0,NROW-1
         B(1+IROW*(ICB+IRB))=2.0D0*B(1+IROW*(ICB+IRB))
20    CONTINUE
      RETURN
      END
