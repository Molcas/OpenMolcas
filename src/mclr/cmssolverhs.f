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
* Copyright (C) 2021, Jie J. Bao                                       *
************************************************************************
* ****************************************************************
* history:                                                       *
* Jie J. Bao, on Aug. 06, 2020, created this file.               *
* ****************************************************************
      Subroutine SolveforRHS(Fock,CICSF,AXkzx,AXPzx,bk,bP)
#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "real.fh"
#include "sa.fh"
#include "crun_mclr.fh"
****** Output
      Real*8,DIMENSION(nDens2+6)::Fock
      Real*8,DIMENSION(nconf1*nroots)::CICSF
****** Input
      Real*8,DIMENSION(nDens2)::AXkzx
      Real*8,DIMENSION(NConf1*nRoots)::AXPzx
      Real*8,DIMENSION(nDens2)::bk
      Real*8,DIMENSION(nConf1*nRoots)::bP
****** Assistants
      INTEGER nRow

*****  Orbital Rotation Part
      nRow=nDens2
      CALL FZero(Fock,nDens2)
      CALL DCopy_(nRow,Axkzx,1,Fock,1)
      CALL DAXPY_(nRow,1.0d0,bk,1,Fock,1)

****** State-CSF Rotation Part
      nRow=nRoots*nConf1
      CALL FZero(CICSF,nRow)
      CALL DCopy_(nRow,AXPzx,1,CICSF,1)
      CALL DAXPY_(nRow,-1.0d0,bP,1,CICSF,1)

      RETURN
      END SUBROUTINE
******************************************************

******************************************************
      subroutine SolveforzX(zX,AXX,bX)
      use stdalloc, only : mma_allocate, mma_deallocate
      use cmslag,   only : ResQaaLag2
#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "real.fh"
#include "sa.fh"
#include "crun_mclr.fh"
#include "warnings.h"
****** Output
      Real*8,DIMENSION((nRoots-1)*nRoots/2)::zX
****** Input
      Real*8,DIMENSION((nRoots-1)*nRoots/2)::bX
      Real*8,DIMENSION(((nRoots-1)*nRoots/2)**2)::AXX
****** Assistants
      Real*8,DIMENSION(:),Allocatable::EigVal,bxscr,zXscr,Scr
      Real*8 TwoPi
      INTEGER NDim,nSPair,iPair,nScr,INFO


      NDim=((nRoots-1)*nRoots/2)
      nSPair=nDim
      TwoPi=2.0d0*Pi
      ResQaaLag2=0.0d0
      CALL mma_allocate(EigVal,nDim)
      CALL mma_allocate(bxScr ,nDim)
      CALL mma_allocate(zXScr ,nDim)

      CALL GetDiagScr(nScr,AXX,EigVal,nDim)
      CALL mma_allocate(Scr   ,nScr)

      CALL DSYEV_('V','U',nDim,AXX,nDim,EigVal,Scr,nScr,INFO)

      CALL DGEMM_('n','n',1,nDim,nDim,1.0d0,bx,1,AXX,nDim,
     &                                0.0d0,bxScr,1)


      DO iPair=1,nDim
       zxScr(iPair)=-bxScr(iPair)/EigVal(iPair)
       IF(Abs(zxScr(iPair)).gt.TwoPi) THEN
        zxScr(iPair)=0.0d0
        ResQaaLag2=ResQaaLag2+bxScr(iPair)**2
       END IF
      END DO

      write(6,'(6X,A37,2X,ES17.9)')
     & 'Residual in Qaa Lagrange Multipliers:',SQRT(ResQaaLag2)
      IF(ResQaaLag2.gt.epsilon**2) THEN
        write(6,*)
        write(6,'(6X,A)')
     &    'ERROR: RESIDUAL(S) FOR INTERMEDIATE STATE TOO BIG!'
        write(6,*)
        write(6,'(6X,A)')
     &    'This may come from a linear molecular or a linear'
        write(6,'(6X,A)')
     &    'fragment.'
        write(6,'(6X,A)')
     &    'CMS-PDFT Lagrange multipliers are not solved.'
        CALL WarningMessage(2,
     &    'Residual in Lagrange Multipliers for Qaa Too Big')
        CALL Quit(_RC_EXIT_EXPECTED_)
      END IF

      CALL DGEMM_('n','t',    1,nSPair,nSPair,
     &            1.0d0,zXScr,1,AXX,nSPair,
     &            0.0d0,zx   ,1)

      CALL mma_deallocate(EigVal)
      CALL mma_deallocate(bxScr )
      CALL mma_deallocate(zXScr )
      CALL mma_deallocate(Scr   )
      RETURN
      END SUBROUTINE
******************************************************
******************************************************
      Subroutine GetQaaFock(FOccMO,P2MOt,GDMat,zX,nP2)
      use stdalloc, only : mma_allocate, mma_deallocate
#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "real.fh"
#include "sa.fh"
#include "crun_mclr.fh"
******Input
      Real*8,DIMENSION((nRoots-1)*nRoots/2)::zX
      Real*8,DIMENSION(nRoots*(nRoots+1)/2,nnA,nnA)::GDMat
      Real*8,DIMENSION(nP2)::P2MOt
******Output
      Real*8,DIMENSION(nDens2)::FOccMO
******For Debugging
      INTEGER NPUVX
      Real*8,DIMENSION(:),Allocatable::PUVX
      INTEGER,DIMENSION(ntAsh,ntAsh,ntAsh,ntAsh)::IndTUVX
      INTEGER,DIMENSION(ntBas,ntAsh,ntAsh,ntAsh)::IndPUVX
      Logical debug2
******Auxiliaries
      Real*8,DIMENSION(:),Allocatable::G1r,G2r,G2q,Fock,T,PQaa
      INTEGER K,L,nG2r,IKL,IKL2,IKK,ILL
************************************************************************
*                                                                      *
       itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
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

******************

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
      RETURN
      End Subroutine
******************************************************

******************************************************
      Subroutine G2qtoG2r(G2r,G2q,nG2,nG2r)
#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "real.fh"
#include "sa.fh"
#include "crun_mclr.fh"
      INTEGER nG2,nG2r
      Real*8,DIMENSION(nG2 )::G2q
      Real*8,DIMENSION(nG2r)::G2r
      INTEGER iB,jB,iDij,iRij,iDkl,iRkl,iijkl,iRijkl
      Real*8 Fact
       itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
      Do iB=1,ntash
       Do jB=1,ntash
        iDij=iTri(ib,jB)
        iRij=jb+(ib-1)*ntash
        Do kB=1,ntash
         Do lB=1,ntash
          iDkl=iTri(kB,lB)
          iRkl=lb+(kb-1)*ntash
          fact=One
          if(iDij.ge.iDkl .and. kB.eq.lB) fact=Two
          if(iDij.lt.iDkl .and. iB.eq.jB) fact=Two
          iijkl=itri(iDij,iDkl)
          iRijkl=itri(iRij,iRkl)
          G2r(iRijkl)=Fact*G2q(iijkl)
         End Do
        End Do
       End Do
      End Do
      RETURN
      End Subroutine
******************************************************

******************************************************
      Subroutine QaaVerif(G2q,ng2,PUVX,NPUVX,IndTUVX)
#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "real.fh"
#include "sa.fh"
#include "crun_mclr.fh"
      INTEGER nG2,nPUVX
      Real*8,DIMENSION(nG2)::G2q
      Real*8,DIMENSION(NPUVX)::PUVX
      INTEGER,DIMENSION(ntAsh,ntAsh,ntAsh,ntAsh)::IndTUVX
      INTEGER I,J,K,L,IJKL,lMax
      Real*8 dQdX

      ijkl=0
      dQdX=0.0d0
      do i=1,nna
        do j=1,i
          do k=1,i
            if(i.eq.k) then
              lmax = j
            else
              lmax = k
            end if
            do l=1,lmax
              ijkl = ijkl + 1
              dQdX=dQdX+G2q(ijkl)*PUVX(IndTUVX(I,J,K,L))
            end do
          end do
        end do
      end do

      write(6,*) 'dQdX in QaaVerif=',dQdX

      RETURN
      End Subroutine
******************************************************

******************************************************
      Subroutine QaaP2MO(G2q,ng2,GDMat,IKL,IKK,ILL)
      use stdalloc, only : mma_allocate, mma_deallocate
#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "real.fh"
#include "sa.fh"
#include "crun_mclr.fh"
******  Input
      INTEGER nG2,IKL,IKK,ILL
      Real*8,DIMENSION(nRoots*(nRoots+1)/2,nnA,nnA)::GDMat
******  Output
      Real*8,DIMENSION(nG2)::G2q
******  Auxiliaries
      INTEGER i,j,k,l,ij,kl,ijkl,nD
      Real*8 Fact
      Real*8,DIMENSION(:),Allocatable::Dsum,Ddif
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)

******* Calculating dQ_aa/dX_KL, original purpose of this subroutine
      nD=nnA*(nnA+1)/2
      CALL mma_allocate(Dsum,nD)
      CALL mma_allocate(Ddif,nD)
      DO i=1,nnA
       Do j=1,i
        ij=iTri(i,j)
        Dsum(ij)=GDMat(IKL,i,j)+GDMat(IKL,j,i)
        Ddif(ij)=GDMat(IKK,i,j)-GDMat(ILL,i,j)
       End Do
      END DO
       ijkl=0
       do i=1,nna
         do j=1,i
           ij = iTri(i,j)
           do k=1,i
             if(i.eq.k) then
               lmax = j
             else
               lmax = k
             end if
             do l=1,lmax
               kl = iTri(k,l)
               ijkl = ijkl + 1
               fact=0.5d0
               if(k.eq.l) fact=0.25d0
               G2q(ijkl)=
     & fact*(Dsum(ij)*Ddif(kl)+Dsum(kl)*Ddif(ij))*2.0d0
             end do
           end do
         end do
       end do
      CALL mma_deallocate(Dsum)
      CALL mma_deallocate(Ddif)
      RETURN
      End Subroutine





