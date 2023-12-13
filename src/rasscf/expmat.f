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
* Copyright (C) 2022, Jie J. Bao                                       *
************************************************************************
* ****************************************************************
* history:                                                       *
* Jie J. Bao, on Feb. 03, 2022, created this file.               *
* ****************************************************************
* Purpose: Return to the exponent of a skew-symmetric matrix
      Subroutine ExpMat(M,A,DimM,LenA)
*     Input
      INTEGER DimM, LenA
      Real*8 A(LenA)
*     Output
      Real*8 M(DimM**2)
*     Explanation:
*     Calculate M=exp(-A)
*     A is stored as a 1 by LenA vector, and M is a DimM by DimM matrix
*     This subroutine assumes LenA=(DimM-1)*DimM/2
*     How it works:
*     First transform A into a skew-symmetric matrix form. Then call
*     ExpMat_Inner to calculate exp(A)
*     Auxiliary:
      Real*8 MatA(DimM**2)
      INTEGER I,J,N,iIJ1,iIJ2

      N=DimM**2
      CALL FZero(MatA,N)
      DO I=2,DimM
       Do J=1,I-1
        iIJ2=(I-2)*(I-1)/2+J
        iIJ1=(I-1)*DimM+J
        MatA(iIJ1)=A(iIJ2)
        iIJ1=(J-1)*DimM+I
        MatA(iIJ1)=-A(iIJ2)
       End Do
      END DO

      CALL ExpMat_Inner(M,MatA,DimM)
      RETURN
      End SUbroutine

      Subroutine ExpMat_Inner(R,X,nLen)
******Purpose: calculate R=exp(X)
******Explanation:
*     The subroutine uses the algorithm in Section 3.1.5 on Page 83 of
*     'Molecular Electronic-Structure Theory' by T. Helgaker, P.
*     J/orgensen and J. Olsen.
*     Jie Bao acknowledges Dr. Chen Zhou from Xiamen University, China, for
*     showing the resource of this algorithm and his code for reference.
      use stdalloc, only : mma_allocate, mma_deallocate
*     Input
      INTEGER nLen,nLen2
      Real*8 x(nLen**2)
*     Output
      Real*8 R(nLen**2)
*     Auxilliary
      Real*8 tau2(nLen),cospart(nLen**2),sinpart(nLen**2),
     &       scr(nLen**2),X2(nLen**2),tau(nLen)
      INTEGER nScrDiag,INFO,I
      Real*8,DIMENSION(:),Allocatable::ScrDiag
      Real*8 Coeff

      nLen2=nLen**2

*Step 1     calculate X2
      CALL DGEMM_('n','n',nLen,nLen,nLen,1.0d0,X,nLen,X,nLen,
     &                                   0.0d0,X2,nLen)

*Step 2     diagonalize X2
      CALL GetDiagScr(nScrDiag,X2,Scr,nLen)
      CALL mma_allocate(ScrDiag,nScrDiag)
      CALL DSYEV_('V','U',nLen,X2,nLen,tau2,ScrDiag,nScrDiag,INFO)
      CALL mma_deallocate(ScrDiag)

      DO I=1,nLen
       tau(I)=dsqrt(dabs(tau2(I)))
      END DO

*Step 3     build cos part of R matrix
      CALL DCopy_(nLen2,X2,1,CosPart,1)
      DO I=1,nLen
       CALL DScal_(nLen,cos(tau(I)),CosPart((I-1)*nLen+1),1)
      END DO

      CALL DGEMM_('n','t',nLen,nLen,nLen,1.0d0,CosPart,nLen,X2,nLen,
     &                                   0.0d0,Scr,nLen)
*     R=W * cos(tau) * W^T
      CALL DCopy_(nLen2,Scr,1,R,1)
*Step 4     build sin part of R matrix
      CALL DCopy_(nLen2,X2,1,SinPart,1)
      DO I=1,nLen
       IF(tau(I).lt.1.0d-8) THEN
        Coeff=1.0d0
       ELSE
        Coeff=sin(tau(I))/tau(I)
       END IF
       CALL DScal_(nLen,Coeff,SinPart((I-1)*nLen+1),1)
      END DO

      CALL DGEMM_('n','t',nLen,nLen,nLen,1.0d0,SinPart,nLen,X2,nLen,
     &                                   0.0d0,Scr,nLen)
*     R=R+W * tau^(-1) * sin(tau) * W^T * X
      CALL DGEMM_('n','n',nLen,nLen,nLen,1.0d0,Scr,nLen,X,nLen,
     &                                   1.0d0,R,nLen)
      RETURN
      End Subroutine
