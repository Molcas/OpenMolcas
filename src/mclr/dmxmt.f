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
      SUBROUTINE DMXMT(A,LDA,NTA,B,LDB,NTB,
     &                C,NRC,NSUM,NCC)
      IMPLICIT NONE
      INTEGER IROW,ICOL,LDA,LDB,NRC,NSUM,NCC,IND,ISUM
      CHARACTER NTA,NTB
      REAL*8 A(LDA,*),B(LDB,*),C(*),SUM
      IND=0
      IF (NTA.eq.'N'.and.NTB.eq.'N') Then
      DO ICol=1,NRC
       DO IRow=icol,nrc
         SUM=0.0D0
         DO ISUM=1,NSUM
            SUM=SUM +
     &          A(IROW,ISUM)*
     &          B(ISUM,iCOL)
         End Do
         IND=IND+1
         C(IND)=SUM
       End Do
      End Do
      Else
       Call SysHalt('dmxmt')
      End If
*     IF (NTA.eq.'T'.and.NTB.eq.'N') Then
*     DO IROW=1,NRC
*      DO ICOL=0,IROW
*        SUM=0.0D0
*        DO ISUM=0,NSUM-1
*           SUM=SUM +
*    &          A(ISUM,IROW)*
*    &          B(ISUM,iCOL)
*        End Do
*        IND=IND+1
*        C(IND)=SUM
*      End Do
*     End Do
*     Else IF (NTA.eq.'N'.and.NTB.eq.'T') Then
*     DO IROW=1,NRC
*      DO ICOL=0,IROW
*        SUM=0.0D0
*        DO ISUM=0,NSUM-1
*           SUM=SUM +
*    &          A(IROW,ISUM)*
*    &          B(ICOL,iSUM)
*        End Do
*        IND=IND+1
*        C(IND)=SUM
*      End Do
*     End Do
*     Else IF (NTA.eq.'T'.and.NTB.eq.'T') Then
*     DO IROW=1,NRC
*      DO ICOL=0,IROW
*        SUM=0.0D0
*        DO ISUM=0,NSUM-1
*           SUM=SUM +
*    &          A(ISUM,IROW)*
*    &          B(ICOL,iSUM)
*        End Do
*        IND=IND+1
*        C(IND)=SUM
*      End Do
*     End Do
*     End If
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(NCC)
      END
