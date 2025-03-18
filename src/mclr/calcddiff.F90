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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************
      SUBROUTINE CalcDdiff(Ddiff,GDMat,M,K,nnA,nRoots)
      Implicit None
      INTEGER nnA,nRoots,M,K
      REAL*8,DIMENSION((nRoots+1)*nRoots/2,nnA,nnA)::GDMat
      REAL*8,DIMENSION(nnA**2)::Ddiff

      INTEGER it,iu,iMM,iKK

      iMM=(M+1)*M/2
      iKK=(K+1)*K/2

      DO it=1,nnA
       Do iu=1,nnA
        Ddiff((it-1)*nnA+iu)=GDMat(iMM,it,iu)-GDMat(iKK,it,iu)
       End Do
      END DO

      END SUBROUTINE CalcDdiff
