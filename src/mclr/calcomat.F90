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
      subroutine CalcOMat(CSFOK,LOK,FMO1t,FMO2t,nTri)
      use ipPage, only: W
      use stdalloc, only: mma_allocate, mma_deallocate
      use MCLR_Data, only: nConf1, nAcPr2, ipCI, ipMat, nDens2
      use MCLR_Data, only: XISPSM
      use input_mclr, only: nRoots,State_Sym,nSym,nBas
      Implicit None

!******Output
      Real*8,DIMENSION(nRoots*nConf1)::CSFOK
      Real*8,DIMENSION(nRoots**2)::LOK
!******Input
      INTEGER nTri
      Real*8,DIMENSION(nRoots*nTri)::FMO1t
      Real*8,DIMENSION(nRoots*nacpr2)::FMO2t
!******A little help from
      Real*8,DIMENSION(1)::rdum
      Real*8,DIMENSION(:),Allocatable::FMO1
      INTEGER ILoc1,ILoc2,ILoc3,iOff,iS,jS,iB,jB,ji,ij,I,iK
      INTEGER ILoc4,iptmp,nConf3, L
      INTEGER, EXTERNAL:: ipGet
      Real*8, External:: DDot_
!                                                                      *
!***********************************************************************
!                                                                      *
      Interface
       SubRoutine CISigma_sa(iispin,iCsym,iSSym,Int1,nInt1,Int2s,nInt2s,&
     &                       Int2a,nInt2a,ipCI1,ipCI2, Have_2_el)
       Integer iispin, iCsym, iSSym
       Integer nInt1, nInt2s, nInt2a
       Real*8, Target:: Int1(nInt1), Int2s(nInt2s), Int2a(nInt2a)
       Integer ipCI1, ipCI2
       Logical Have_2_el
       End SubRoutine CISigma_sa
      End Interface
!                                                                      *
!***********************************************************************
!                                                                      *


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
!       irc=ipin(ipCI)
       CALL CISigma_SA(0,State_Sym,State_Sym,FMO1,nDens2,               &
     & FMO2t(iLoc2),NACPR2,rdum,1,ipci,iptmp,.True.)
       CALL Daxpy_(nConf1,Real(nRoots,8),W(iptmp)%Vec(iLoc3),1,         &
     &                                          CSFOK(iLoc3),1)

       Do L=1,nRoots
        ILoc4=(L-1)*NConf1+1
        LOK((I-1)*nRoots+L)=ddot_(nConf1,CSFOK(iLoc3),1,                &
     &                             W(ipCI)%Vec(iLoc4),1)
       End Do
      END DO
      CALL mma_deallocate(FMO1)
      End Subroutine CalcOMat
