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
      Subroutine CalcbX(bX,LOK,R,H)
      use Constants, only: Zero
      use MCLR_Data, only: IRLXROOT
      use input_mclr, only: nRoots
      Implicit None
!***** Output
      Real*8,DIMENSION((nRoots-1)*nRoots/2)::bX
!***** Input
      Real*8,DIMENSION(nRoots**2)::R,H
      Real*8,DIMENSION(nRoots**2)::LOK
!**** Auxiliaries
      INTEGER I,K,L,M,N,IKL,IIM,IIN,IKOL,IIK,IIL,ILOK
      Real*8 TempD

      bX(:)=Zero
      I=irlxroot
      DO K=2,nRoots
       IIK=(I-1)*nRoots+K
      DO L=1,K-1
       IIL=IIK-K+L
       IKL=(K-2)*(K-1)/2+L
       IKOL=(L-1)*nRoots+K
       ILOK=(K-1)*nRoots+L
       bX(IKL)=R(IIK)**2*LOK(ILOK)-R(IIL)**2*LOK(IKOL)
       Do M=2,nRoots
        IIM=IIK-K+M
       Do N=1,M-1
        TempD=0.0d0
        IIN=IIK-K+N
        IF(M.eq.K) TempD=TempD+H((L-1)*nRoots+N)
        IF(N.eq.K) TempD=TempD+H((M-1)*nRoots+L)
        IF(M.eq.L) TempD=TempD-H((K-1)*nRoots+N)
        IF(N.eq.L) TempD=TempD-H((M-1)*nRoots+K)
        bX(IKL)=bX(IKL)+TempD*R(IIM)*R(IIN)
       End Do
       End Do
       bX(IKL)=bX(IKL)*2.0d0
      END DO
      END DO
      END SUBROUTINE CalcbX
