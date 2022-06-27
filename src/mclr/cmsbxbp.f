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
      subroutine CalcbXbP(bX,bP,FMO1t,FMO2t,R,H,nTri)
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


****** Output
       Real*8,DIMENSION((nRoots-1)*nRoots/2)::bX
       Real*8,DIMENSION(nConf1*nRoots)::bP
****** Input
       INTEGER nTri
       Real*8,DIMENSION(nRoots*nTri)::FMO1t
       Real*8,DIMENSION(nRoots*nacpr2)::FMO2t
       Real*8,DIMENSION(nRoots**2)::R,H
****** Auxiliaries
       Real*8,DIMENSION(:),Allocatable::LOK,CSFOK

       CALL mma_allocate(CSFOK,nRoots*nConf1)
       CALL mma_allocate(LOK,nRoots**2)
       CALL CalcOMat(CSFOK,LOK,FMO1t,FMO2t,nTri)
       CALL CalcbP(bP,CSFOK,LOK,R)
       CALL CalcbX(bX,LOK,R,H)
       CALL mma_deallocate(CSFOK)
       CALL mma_deallocate(LOK)
       return
       end subroutine
******************************************************

      Subroutine CalcbX(bX,LOK,R,H)
#include "stdalloc.fh"
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
****** Output
      Real*8,DIMENSION((nRoots-1)*nRoots/2)::bX
****** Input
      Real*8,DIMENSION(nRoots**2)::R,H
      Real*8,DIMENSION(nRoots**2)::LOK
***** Auxiliaries
      INTEGER I,K,L,M,N,IKL,IIM,IIN,IKOL,IIK,IIL
      Real*8 TempD

      CALL FZero(bX,(nRoots-1)*nRoots/2)
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
      RETURN
      END SUBROUTINE
******************************************************


******************************************************
      subroutine CalcbP(bP,CSFOK,LOK,R)
      use ipPage, only: W
#include "stdalloc.fh"
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
***** Output
      Real*8,DIMENSION(nConf1*nRoots)::bP
***** Input
      Real*8,DIMENSION(nRoots*nConf1)::CSFOK
      Real*8,DIMENSION(nRoots**2)::LOK
      Real*8,DIMENSION(nRoots**2)::R
***** Kind quantities that help
      INTEGER I,L,K,iLoc1,iLoc2
      Real*8 tempd
      I=irlxroot
      DO K=1,nRoots
       iLoc1=(K-1)*nConf1+1
       CALL DCopy_(nConf1,CSFOK(iLoc1),1,bP(iLoc1),1)
       Do L=K,K
        tempd=-LOK((K-1)*nRoots+L)
        iLoc2=(L-1)*nConf1+1
        CALL dAXpY_(nConf1,tempd,W(ipci)%Vec(iLoc2),1,bP(iLoc1),1)
       End Do

       Do L=1, nRoots
        IF (L.eq.K) CyCle
        tempd=-LOK((K-1)*nRoots+L)
        iLoc2=(L-1)*nConf1+1
        CALL dAXpY_(nConf1,tempd,W(ipci)%Vec(iLoc2),1,bP(iLoc1),1)
       End Do
       CALL DScal_(nConf1,2.0d0*R((I-1)*nRoots+K)**2,
     & bP(iLoc1),1)
      END DO
      RETURN
      End Subroutine
******************************************************

******************************************************
      subroutine CalcOMat(CSFOK,LOK,FMO1t,FMO2t,nTri)
      use ipPage, only: W
      Implicit Real*8 (a-h,o-z)
#include "stdalloc.fh"
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

*******Output
      Real*8,DIMENSION(nRoots*nConf1)::CSFOK
      Real*8,DIMENSION(nRoots**2)::LOK
*******Input
      INTEGER nTri
      Real*8,DIMENSION(nRoots*nTri)::FMO1t
      Real*8,DIMENSION(nRoots*nacpr2)::FMO2t
*******A little help from
      Real*8,DIMENSION(1)::rdum
      Real*8,DIMENSION(:),Allocatable::FMO1
      INTEGER ILoc1,ILoc2,ILoc3,iOff,iS,jS,iB,jB,ji,ij,I,iK
      INTEGER ILoc4,iptmp,nConf3
*                                                                      *
************************************************************************
*                                                                      *
      Interface
       SubRoutine CISigma_sa(iispin,iCsym,iSSym,Int1,nInt1,Int2s,nInt2s,
     &                       Int2a,nInt2a,ipCI1,ipCI2, Have_2_el)
       Integer iispin, iCsym, iSSym
       Integer nInt1, nInt2s, nInt2a
       Real*8, Target:: Int1(nInt1), Int2s(nInt2s), Int2a(nInt2a)
       Integer ipCI1, ipCI2
       Logical Have_2_el
       End SubRoutine CISigma_sa
      End Interface
*                                                                      *
************************************************************************
*                                                                      *


      nConf3=nint(Max(xispsm(State_SYM,1),xispsm(State_SYM,1)))
      iptmp=ipGet(nconf3*nRoots)
      CALL FZero(CSFOK,nRoots*nConf1)
      CALL FZero(W(iptmp)%Vec,nRoots*nConf3)
      CALL mma_allocate(FMO1,nDens2)

      DO i=1,nRoots
       iK=I
       ILoc1=  (IK-1)*nTri
       ILoc2=1+(IK-1)*NACPR2
       ILoc3=1+(IK-1)*nConf1
       iOff=0
       Do iS=1,nSym
         jS=iS
         IF (NBAS(IS)*NBAS(JS).NE.0) THEN
          do iB=1,nBas(iS)
            do jB=1,iB
             ioff = ioff+1
             ji= ipMat(is,js)-1 +(iB-1)*nbas(iS)+jB
             FMO1(ji) = FMO1t(ioff+iLoc1)
             if (iB.ne.jB) then
              ij= ipMat(is,js)-1 +(jB-1)*nbas(iS)+iB
              FMO1(ij) = FMO1t(ioff+iLoc1)
             end if
            end do
          end do
         END IF
       End do
*       irc=ipin(ipCI)
       CALL CISigma_SA(0,State_Sym,State_Sym,FMO1,nDens2,
     & FMO2t(iLoc2),NACPR2,rdum,1,ipci,iptmp,.True.)
       CALL Daxpy_(nConf1,Real(nRoots,8),W(iptmp)%Vec(iLoc3),1,
     &                                          CSFOK(iLoc3),1)

       Do L=1,nRoots
        ILoc4=(L-1)*NConf1+1
        LOK((I-1)*nRoots+L)=ddot_(nConf1,CSFOK(iLoc3),1,
     &                             W(ipCI)%Vec(iLoc4),1)
       End Do
      END DO
      CALL mma_deallocate(FMO1)
      RETURN
      End Subroutine
******************************************************
