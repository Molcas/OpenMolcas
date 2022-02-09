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
      SUBROUTINE MXMT(A,ICA,IRA, B,ICB,IRB, C, NROW,NSUM)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),B(*),C(*)
#define _ALTERNATIVE_
#ifdef _ALTERNATIVE_
      Character(LEN=1) OPA, OPB
      Integer :: LDA, LDB

      If (ICA.eq.1) Then
        OPA='N'
        LDA=IRA
      Else If (IRA.eq.1) Then
        OPA='T'
        LDA=ICA
      Else
        Call AbENd()
      End If

      If (ICB.eq.1) Then
        OPB='N'
        LDB=IRB
      Else If (IRB.eq.1) Then
        OPB='T'
        LDB=ICB
      Else
        Call AbENd()
      End If

      Call DGEMM_Tri(OPA,OPB,NROW,NROW,NSUM,
     &               1.0D0,A,LDA,
     &                     B,LDB,
     &               0.0D0,C,NROW)
#else
      IND=0
      DO 100 IROW=0,NROW-1
      DO 101 ICOL=0,IROW
         SUM=0.0D0
         DO 110 ISUM=0,NSUM-1
            AB=A(1+IROW*ICA+ISUM*IRA) *
     &         B(1+ISUM*ICB+ICOL*IRB)
            SUM=SUM+AB
110      CONTINUE
         IND=IND+1
         C(IND)=SUM
101   CONTINUE
100   CONTINUE
#endif
      RETURN
      END
