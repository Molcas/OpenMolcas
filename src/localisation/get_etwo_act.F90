!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine Get_Etwo_act(Dma,Dmb,nBDT,nBas,nSym,Etwo)

      Implicit Real*8 (a-h,o-z)
      Integer nBDT, nBas(8)
      Real*8  Dma(nBDT), Dmb(nBDT)
#include "WrkSpc.fh"
#include "choscf.fh"
#include "choscreen.fh"
#include "chotime.fh"
      Integer ipFLT(2), ipKLT(2), nIorb(8,2), ipPorb(2)
      Integer ipDm(2)
!      Real*8   Get_ExFac
!      External Get_ExFac
!      Character*16 KSDFT
!
      timings=.false.
      Estimate=.false.
      REORD=.false.
!
      Update=.true.
      DECO=.true.
      dmpk=1.0d0
      dFKmat=0.0d0
      ALGO=4
      NSCREEN=10
!
!      nDMat=2
      nBB=0
      Do i=1,nSym
!         nForb(i,1)=0
!         nForb(i,2)=0
         nBB=nBB+nBas(i)**2
      End Do
!      Call Get_cArray('DFT functional',KSDFT,16)
!      ExFac=Get_ExFac(KSDFT)
!      FactXI=1.0d0*ExFac
!      FactXI=1.0d0  ! always HF energy
      Call GetMem('PLTc','Allo','Real',ipPLT,nBDT)
      call dcopy_(nBDT,Dma,1,Work(ipPLT),1)
      Call daxpy_(nBDT,1.0d0,Dmb,1,Work(ipPLT),1)
!
      Call GetMem('ChMc','Allo','Real',ipPorb(1),2*nBB)
      ipPorb(2)=ipPorb(1)+nBB
      Call GetMem('DSQc','Allo','Real',ipDm(1),2*nBB)
      ipDm(2)=ipDm(1)+nBB
      Call UnFold(Dma,nBDT,Work(ipDm(1)),nBB,nSym,nBas)
      Call UnFold(Dmb,nBDT,Work(ipDm(2)),nBB,nSym,nBas)
      iOff=0
      Do i=1,nSym
         ipV=ipPorb(1)+iOff
         ipDai=ipDm(1)+iOff
         Call CD_InCore(Work(ipDai),nBas(i),Work(ipV),nBas(i),          &
     &                  nIorb(i,1),1.0d-12,irc)
         If (irc.ne.0) Then
            write(6,*) ' Alpha density. Sym= ',i,'   rc= ',irc
            Call Abend()
         EndIf
         ipV=ipPorb(2)+iOff
         ipDbi=ipDm(2)+iOff
         Call CD_InCore(Work(ipDbi),nBas(i),Work(ipV),nBas(i),          &
     &                  nIorb(i,2),1.0d-12,irc)
         If (irc.ne.0) Then
            write(6,*) ' Beta density. Sym= ',i,'   rc= ',irc
            Call Abend()
         EndIf
         iOff=iOff+nBas(i)**2
      End Do
!
      Call GetMem('FCNO','Allo','Real',ipFCNO,2*nBDT)
      Call FZero(Work(ipFCNO),2*nBDT)
      ipFLT(1)=ipFCNO
      ipFLT(2)=ipFLT(1)+nBDT
      Call GetMem('KLTc','Allo','Real',ipKLT(1),2*nBDT)
      Call FZero(Work(ipKLT(1)),2*nBDT)
      ipKLT(2)=ipKLT(1)+nBDT
!
      Call Cho_X_init(irc,ChFracMem)
      if (irc.ne.0) then
         Call WarningMessage(2,'Get_CNOs. Non-zero rc in Cho_X_init.')
         Call Abend
      endif

! BIGOT FIXME
      Call WarningMessage(2,                                            &
     &     'There is probably a bug here, ipPLT should have two '//     &
     &     'elements.')
      Call Abend()
!     Call CHO_LK_SCF(irc,nDMat,ipFLT,ipKLT,nForb,nIorb,
!    &                    ipPorb,ipPLT,FactXI,nSCReen,dmpk,dFmat)
      if (irc.ne.0) then
         Call WarningMessage(2,'Get_CNOs. Non-zero rc in Cho_LK_scf.')
         CALL Abend
      endif

      Call Cho_X_Final(irc)
      if (irc.ne.0) then
         Call WarningMessage(2,'Get_CNOs. Non-zero rc in Cho_X_Final.')
         CALL Abend
      endif
!
      Etwo = 0.5d0*(ddot_(nBDT,Dma,1,Work(ipFLT(1)),1)                  &
     &     +        ddot_(nBDT,Dmb,1,Work(ipFLT(2)),1))
!
      Call GetMem('KLTc','Free','Real',ipKLT(1),2*nBDT)
      Call GetMem('FCNO','Free','Real',ipFCNO,2*nBDT)
      Call GetMem('DSQc','Free','Real',ipDm(1),2*nBB)
      Call GetMem('ChMc','Free','Real',ipPorb(1),2*nBB)
      Call GetMem('PLTc','Free','Real',ipPLT,nBDT)
!
      Return
      End
