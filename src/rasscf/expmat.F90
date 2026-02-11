!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Feb. 03, 2022, created this file.               *
! ****************************************************************
! Purpose: Return to the exponent of a skew-symmetric matrix
      Subroutine ExpMat(M,A,DimM,LenA)
!     Input
      INTEGER DimM, LenA
      Real*8 A(LenA)
!     Output
      Real*8 M(DimM**2)
!     Explanation:
!     Calculate M=exp(-A)
!     A is stored as a 1 by LenA vector, and M is a DimM by DimM matrix
!     This subroutine assumes LenA=(DimM-1)*DimM/2
!     How it works:
!     First transform A into a skew-symmetric matrix form. Then call
!     ExpMat_Inner to calculate exp(A)
!     Auxiliary:
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
