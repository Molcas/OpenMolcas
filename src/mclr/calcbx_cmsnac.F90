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
! Copyright (C) 2021, Paul B Calio                                     *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Based on cmsbxbp.f from Jie J. Bao                             *
! ****************************************************************
      Subroutine CalcbX_CMSNAC(bX,LOK,R,H,E_Final)
      use Constants, only: Zero
      use MCLR_Data, only: ISMECIMSPD,NACSTATES
      use input_mclr, only: nRoots
      Implicit None
!***** Output
      Real*8,DIMENSION((nRoots-1)*nRoots/2)::bX
!***** Input
      Real*8,DIMENSION(nRoots**2)::R,H
      Real*8,DIMENSION(nRoots)::E_Final
      Real*8,DIMENSION(nRoots**2)::LOK
!**** Auxiliaries
      INTEGER I,J,K,L,M,N,IKL,IIN,IJM,IKOL,IIK,IJK,IIL,IJL,ILOK
      Real*8 TempD, dE_IJ

      bX(:)=Zero
      I=NACstates(1)
      J=NACstates(2)
      dE_IJ = E_Final(I)-E_Final(J)

!***** R_JK*HKL*RIL
      DO K=2,nRoots
       IIK=(I-1)*nRoots+K
       IJK=(J-1)*nRoots+K
      DO L=1,K-1
       IIL=IIK-K+L
       IJL=IJK-K+L
       IKL=(K-2)*(K-1)/2+L
       IKOL=(L-1)*nRoots+K
       ILOK=(K-1)*nRoots+L
!******Diagonal elements of R_JK * H_KL * R_IK
       bX(IKL)=2.0d0*(R(IJK)*R(IIK)*LOK(ILOK)-R(IJL)*R(IIL)*LOK(IKOL))
!****** Additional NAC term (Requires only one line)
!****** R_JK * <K|L> * R_IL =
!****** R_JL * R_IK - R_JK * R_IL
       if(.not.isMECIMSPD) bX(IKL) = bX(IKL) +                          &
     & dE_IJ * ( R(IJL)*R(IIK)-R(IJK)*R(IIL) )

!******Off-Diagonal elements of R_JK * H_KL * R_IK
       Do M=1,nRoots
        IJM=IJK-K+M
       Do N=1,nRoots
        if(M.eq.N) cycle
        TempD=0.0d0
        IIN=IIK-K+N
        IF(M.eq.K) TempD=TempD+H((L-1)*nRoots+N)
        IF(N.eq.K) TempD=TempD+H((M-1)*nRoots+L)
        IF(M.eq.L) TempD=TempD-H((K-1)*nRoots+N)
        IF(N.eq.L) TempD=TempD-H((M-1)*nRoots+K)
        bX(IKL)=bX(IKL)+TempD*R(IJM)*R(IIN)
       End Do
       End Do
      END DO
      END DO
      END SUBROUTINE CalcbX_CMSNAC
