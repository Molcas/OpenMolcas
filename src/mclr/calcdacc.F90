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
      Subroutine CalcDacc(Dacc,GDMat,M,nnA,nRoots,zx)
      use Constants, only: Zero
      Implicit None
      INTEGER nnA,nRoots,M
      REAL*8,DIMENSION((nRoots+1)*nRoots/2,nnA,nnA)::GDMat
      REAL*8,DIMENSION(nnA**2)::Dacc
      REAL*8,DIMENSION((nRoots-1)*nRoots/2)::zx

      INTEGER it,iu,K,IKM,IKM2,iLoc1,iLoc2
      REAL*8 Fact

      Dacc(:)=Zero

      DO K=1,nRoots
       IF(K.eq.M) Cycle
       IF(M.gt.K) THEN
        IKM =(M-1)*M/2+K
        IKM2=(M-1)*(M-2)/2+K
       ELSE
        IKM =(K-1)*K/2+M
        iKM2=(K-1)*(K-2)/2+M
       END IF
       Fact=4.0d0*zx(IKM2)
       IF(K.gt.M) Fact=-Fact
       Do it=1,nnA
        iLoc1=(it-1)*nnA
        do iu=1,nnA
         iLoc2=iLoc1+iu
         Dacc(iLoc2)=Dacc(iLoc2)+GDMat(IKM,it,iu)*Fact
        end do
       End Do
      END DO
      END SUBROUTINE CalcDacc
