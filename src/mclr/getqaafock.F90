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
      Subroutine GetQaaFock(FOccMO,P2MOt,GDMat,zX,nP2)
      use stdalloc, only : mma_allocate, mma_deallocate
      use MCLR_Data, only: nNA, nDens2
      use input_mclr, only: nRoots,ntAsh,ntBas
      Implicit None
!*****Input
      Integer nP2
      Real*8,DIMENSION((nRoots-1)*nRoots/2)::zX
      Real*8,DIMENSION(nRoots*(nRoots+1)/2,nnA,nnA)::GDMat
      Real*8,DIMENSION(nP2)::P2MOt
!*****Output
      Real*8,DIMENSION(nDens2)::FOccMO
!*****For Debugging
      INTEGER NPUVX
      Real*8,DIMENSION(:),Allocatable::PUVX
      INTEGER,DIMENSION(ntAsh,ntAsh,ntAsh,ntAsh)::IndTUVX
      INTEGER,DIMENSION(ntBas,ntAsh,ntAsh,ntAsh)::IndPUVX
      Logical debug2
!*****Auxiliaries
      Real*8,DIMENSION(:),Allocatable::G1r,G2r,G2q,Fock,T,PQaa
      INTEGER K,L,nG2r,IKL,IKL2,IKK,ILL, nG1, nG2, nG1r
!***********************************************************************
!                                                                      *
       Integer i,j,itri
       itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
!                                                                      *
!***********************************************************************
      ng1=itri(ntash,ntash)
      ng2=itri(ng1,ng1)

      nG1r=ntash**2
      nG2r=(nG1r+1)*nG1r/2
      CALL mma_allocate(Fock,nDens2)
      CALL mma_allocate(T,nDens2)
      CALL mma_allocate(G1r,nG1r)
      CALL mma_allocate(G2r,nG2r)
      CALL mma_allocate(G2q,ng2)
      CALL mma_allocate(PQaa,ng2)
      Debug2=.false.
      CALL FZero(PQaa,nP2)
      CALL FZero(G1r,nG1r)
      CALL FZero(G2r,nG2r)

!*****************

      IF(Debug2) THEN
      CALL Get_PUVXLen(NPUVX)
      CALL mma_allocate(PUVX,NPUVX)
      CALL Get_Two_Ind(IndPUVX,IndTUVX)
      END IF

      DO K=1,nRoots
       IKK=(K+1)*K/2
       Do L=1,K-1
        ILL=(L+1)*L/2
        IKL=(K-1)*K/2+L
        IKL2=(K-1)*(K-2)/2+L
        CALL QaaP2MO(G2q,ng2,GDMat,IKL,IKK,ILL)
        IF(Debug2) CALL QaaVerif(G2q,ng2,PUVX,NPUVX,IndTUVX)
        CALL G2qtoG2r(G2r,G2q,nG2,nG2r)
        Call daxpy_(ng2,zX(IKL2),G2q,1,PQaa,1)
       End Do
      END DO


      CALL Daxpy_(nG2,1.0d0,PQaa,1,P2MOt,1)
      CALL Put_dArray('P2MOt',P2MOt,nG2)

      CALL G2qtoG2r(G2r,PQaa,nG2,nG2r)
      CALL FockGen(0.0d0,G1r,G2r,T,Fock,1)
      Call DAxPy_(nDens2,1.0d0,T,1,FOccMO,1)
      IF(Debug2) Call mma_deallocate(PUVX)
      CALL mma_deallocate(Fock)
      CALL mma_deallocate(T)
      CALL mma_deallocate(G1r)
      CALL mma_deallocate(G2r)
      CALL mma_deallocate(G2q)
      CALL mma_deallocate(PQaa)
      End Subroutine GetQaaFock
