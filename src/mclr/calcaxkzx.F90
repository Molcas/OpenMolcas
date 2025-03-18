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
      subroutine CalcAXkzx(AXkzx,GDMat,PUVX,NPUVX,IndPUVX,zx)
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: Zero
      use MCLR_Data, only: nNA, nDens2
      use input_mclr, only: nRoots,ntBas,ntAsh,nSym,nAsh,nOrb
      Implicit None

      Integer NPUVX
      Real*8,DIMENSION((nRoots+1)*nRoots/2,nnA,nnA),Intent(In)::GDMat
      Real*8,DIMENSION(NPUVX),Intent(In)::PUVX
      INTEGER,DIMENSION(ntBas,ntAsh,ntAsh,ntAsh),Intent(In)::IndPUVX
      Real*8,DIMENSION((nRoots-1)*nRoots/2),Intent(In)::zx
      Real*8,DIMENSION(nDens2), Intent(Out)::AXkzx

!*****Auxiliary Quantities
      INTEGER,DIMENSION(nSym)::Off_Act,Off_Orb
      Real*8,DIMENSION(:),Allocatable::DKL1,DKL2,AXktmp
      INTEGER K,L,iKL,iKL2,iKK,iLL
      INTEGER p,q,nTOrb, iSym, i, j, iTri

      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

      Off_Act(1)=0
      Off_Orb(1)=0
      ntOrb=nOrb(1)
      DO ISym=2,nSym
       Off_Act(ISym)=Off_Act(ISym-1)+nAsh(iSym-1)
       Off_Orb(ISym)=Off_Orb(ISym-1)+nOrb(iSym-1)
       ntOrb=ntOrb+nOrb(ISym)
      END DO

      AXkzx(:)=Zero

      CALL mma_allocate(DKL1,ntAsh**2)
      CALL mma_allocate(DKL2,ntAsh**2)
      CALL mma_allocate(AXktmp,nDens2)

      DO K=2,nRoots
       Do L=1,K-1
       iKL=itri(K,L)
       iKK=itri(K,K)
       iLL=itri(L,L)
       iKL2=(K-1)*(K-2)/2+L
       do p=1,ntash
        do q=1,ntash
         DKL1((p-1)*ntash+q)=GDMat(IKL,p,q)+GDMat(IKL,q,p)
         DKL2((p-1)*ntash+q)=GDMat(IKK,p,q)-GDMat(ILL,p,q)
        end do
       end do
       AXktmp(:)=Zero
       CALL CalcAXk2(AXktmp,DKL1,DKL2,PUVX,                             &
     & NPUVX,IndPUVX,Off_Act,Off_Orb)
       CALL CalcAXk2(AXktmp,DKL2,DKL1,PUVX,                             &
     & NPUVX,IndPUVX,Off_Act,Off_Orb)
       CALL Daxpy_(nDens2,zx(IKL2),AXktmp,1,Axkzx,1)
       End Do
      END DO

      CALL mma_deallocate(DKL1)
      CALL mma_deallocate(DKL2)
      CALL mma_deallocate(Axktmp)
      END SUBROUTINE CalcAXkzx
