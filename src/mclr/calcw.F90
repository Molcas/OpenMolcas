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
      subroutine CalcW(W,GDMat,PUVX,NPUVX,IndTUVX)
      use Constants, only: Zero
      use MCLR_Data, only: nNA
      use input_mclr, only: nRoots,ntAsh
      Implicit None

!*****Output
      Real*8,DIMENSION((nRoots+1)*nRoots/2,(nRoots+1)*nRoots/2)::W
!*****Input
      Integer NPUVX
      Real*8,DIMENSION((nRoots+1)*nRoots/2,nnA,nnA)::GDMat
      Real*8,DIMENSION(NPUVX)::PUVX
      INTEGER,DIMENSION(ntAsh,ntAsh,ntAsh,ntAsh)::IndTUVX
!*****Auxiliary Quantities
      INTEGER K,L,M,N,IKL,IMN,it,iu,iv,ix

      DO K=1,nRoots
       DO L=1,K
       IKL=(K-1)*K/2+L
       Do M=1,nRoots
        Do N=1,M
         IMN=(M-1)*M/2+N
         W(IKL,IMN)=Zero
         do it=1,nnA
          do iu=1,nnA
           do iv=1,nnA
            do ix=1,nnA
             IF(IndTUVX(it,iu,iv,ix).ne.0) THEN
            W(IKL,IMN)=W(IKL,IMN)+GDMat(IKL,it,iu)*GDMat(IMN,iv,ix)*    &
     &       PUVX(IndTUVX(it,iu,iv,ix))
             END IF
            end do
           end do
          end do
         end do
        End Do
       End Do
       END DO
      END DO

      End Subroutine CalcW
