************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine Get_CNOs(irc,nIF,nRASO,xNrm)
************************************************************************
*                                                                      *
*     purpose: Generate constrained orbitals from INPORB               *
*              Analysis of the spin configurations                     *
*                                                                      *
************************************************************************
*
      Implicit Real*8 (a-h,o-z)
*
#include "real.fh"
#include "Molcas.fh"
#include "inflocal.fh"
#include "WrkSpc.fh"
      Integer irc
      Integer nRASO(nSym), nIF(nSym)
      Real*8  xNrm
      Character Line*62
      Integer  iOffS(0:8), indxC(16,2,8), nOcc(8)
      Integer  Cho_irange
      External Cho_irange
      Logical  DoneCholesky
************************************************************************
      Match(k,i) = iWork(ipMatch-1+2*(i-1)+k)
************************************************************************
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*
      Call DecideonCholesky(DoneCholesky)
      If(.not.DoneCholesky) then
       write(6,*) '*** Constrained NOs implemented only with CD or RI.'
       write(6,*) '*** Use Cholesky or RICD in Seward and rerun! *****'
       Call Abend
      Endif
*
      irc=0
      nSconf=1
      iOffS(0)=0
      nBB=0
      nBT=0
      nBLT=0
      MaxBas=0
      lConstr=0
      Do iSym=1,nSym
         iOffS(iSym)=iOffS(iSym-1)+nConstr(iSym)
         Do j=1,nConstr(iSym)
            nSconf=2*nSconf
         End Do
         lConstr=lConstr+nConstr(iSym)
         nOcc(iSym)=nIF(iSym)+nConstr(iSym)
         nBB=nBB+nBas(iSym)**2
         nBT=nBT+nBas(iSym)
         nBLT=nBLT+nBas(iSym)*(nBas(iSym)+1)/2
         MaxBas=Max(MaxBas,nBas(iSym))
      End Do
      write(6,'(A,I6)') ' Total number of spin configurations: ',nSconf
      write(6,*)
*
      Call GetMem('Occb','Allo','Real',mAdOcc_ab,nBT)
      Call GetMem('CMOb','Allo','Real',mAdCMO,2*nBB)
      mAdCMO_ab=mAdCMO+nBB
      CALL GetMem('DLT','ALLO','Real',ip_Da,2*nBLT)
      ip_Db=ip_Da+nBLT
*
      Call GetMem('Match','Allo','Inte',ipMatch,2*MxConstr)
      Call GetMem('Corb','Allo','Real',ipCorb,MaxBas)

      Do iCount=0,nSconf-1
         Do jCount=0,lConstr-1
            jSym=Cho_Irange(jCount+1,iOffS,nSym,.true.)
            lCount=jCount-iOffS(jSym-1)+1
            kbit=ibits(iCount,jCount,1)
            If (kbit.eq.0) Then
               indxC(lCount,1,jSym)=1
               indxC(lCount,2,jSym)=2
            Else
               indxC(lCount,1,jSym)=2
               indxC(lCount,2,jSym)=1
            EndIf
         End Do
         write(6,*) ' -------------------------------------------------'
         write(6,*) ' Configuration of the constrained spins (up/down) '
         write(6,*) ' -------------------------------------------------'
         write(6,'(1X,A,I3)') ' nr ',iCount+1
         Do iSym=1,nSym
            write(6,'(1X,A,I1)') ' sym: ',iSym
            Line(1:14)='         (+) '
            k=15
            Do j=1,nConstr(iSym)
               If (indxC(j,1,iSym).eq.1) Then
                  Line(k:k+2)=' u '
               ElseIf (indxC(j,1,iSym).eq.2) Then
                  Line(k:k+2)=' d '
               Else
                  Line(k:k+2)='   '
               EndIf
               k=k+3
            End Do
            write(6,*) Line(1:k-1)
            Line(1:14)='         (-) '
            k=15
            Do j=1,nConstr(iSym)
               If (indxC(j,2,iSym).eq.1) Then
                  Line(k:k+2)=' u '
               ElseIf (indxC(j,2,iSym).eq.2) Then
                  Line(k:k+2)=' d '
               Else
                  Line(k:k+2)='   '
               EndIf
               k=k+3
            End Do
            write(6,*) Line(1:k-1)
         End Do
         write(6,*) ' -------------------------------------------------'
*
         xNrm = 0.0d0
         iOff=0
         jOff=0
         Do iSym=1,nSym
            call dcopy_(nBas(iSym)**2,Work(ipCMO+iOff),1,
     &                 Work(mAdCMO+iOff),1)
            call dcopy_(nBas(iSym)**2,Work(mAdCMO+iOff),1,
     &                 Work(mAdCMO_ab+iOff),1)
            lOcc_=ipOcc+jOff+nIF(iSym)
            call dcopy_(nRASO(iSym),Work(lOcc_),1,Work(mAdOcc_ab),1)
            Call BestMatch(nConstr(iSym),nRASO(iSym),Work(mAdOcc_ab),
     &                     iWork(ipMatch),MxConstr)
            Do i=1,nConstr(iSym)
               k=Match(1,i)
               jOcc=ipOcc-1+jOff+nIF(iSym)+k
               xOkk=Work(jOcc)/2.0d0
               kc=mAdCMO+iOff+nBas(iSym)*(nIF(iSym)+k-1)
               xNrm=xNrm+ddot_(nBas(iSym),Work(kc),1,Work(kc),1)
               l=Match(2,i)
               iOcc=ipOcc-1+jOff+nIF(iSym)+l
               yOkk=Work(iOcc)/2.0d0
               xnorm=sqrt(abs(xOkk)+abs(yOkk)) !ensures correct normaliz
               lc=mAdCMO+iOff+nBas(iSym)*(nIF(iSym)+l-1)
               xOkk=sqrt(abs(xOkk))/xnorm
               yOkk=sqrt(abs(yOkk))/xnorm
               call dscal_(nBas(iSym),xOkk,Work(kc),1)
               call dscal_(nBas(iSym),yOkk,Work(lc),1)
               call dcopy_(nBas(iSym),Work(lc),1,Work(ipCorb),1)
               Call daxpy_(nBas(iSym),1.0d0,Work(kc),1,Work(ipCorb),1)
               Call daxpy_(nBas(iSym),-1.0d0,Work(kc),1,Work(lc),1)
               call dscal_(nBas(iSym),-1.0d0,Work(lc),1)
               call dcopy_(nBas(iSym),Work(ipCorb),1,Work(kc),1)
            End Do
            jc=1
            kc=nConstr(iSym)+1
            Do i=1,nConstr(iSym)
               l=Match(indxC(i,2,iSym),i)
               lc1=mAdCMO+iOff+nBas(iSym)*(nIF(iSym)+l-1)
               lc2=mAdCMO_ab+iOff+nBas(iSym)*(nIF(iSym)+jc-1)
               call dcopy_(nBas(iSym),Work(lc1),1,Work(lc2),1)
               k=Match(indxC(i,1,iSym),i)
               kc1=mAdCMO+iOff+nBas(iSym)*(nIF(iSym)+k-1)
               kc2=mAdCMO_ab+iOff+nBas(iSym)*(nIF(iSym)+kc-1)
               call dcopy_(nBas(iSym),Work(kc1),1,Work(kc2),1)
               jc=jc+1
               kc=kc+1
            End Do
            kc=nConstr(iSym)+1
            Do i=1,nConstr(iSym)
               ic1=mAdCMO_ab+iOff+nBas(iSym)*(nIF(iSym)+i-1)
               ic2=mAdCMO+iOff+nBas(iSym)*(nIF(iSym)+kc-1)
               call dcopy_(nBas(iSym),Work(ic1),1,Work(ic2),1)
               kc1=mAdCMO_ab+iOff+nBas(iSym)*(nIF(iSym)+kc-1)
               kc2=mAdCMO+iOff+nBas(iSym)*(nIF(iSym)+i-1)
               call dcopy_(nBas(iSym),Work(kc1),1,Work(kc2),1)
               kc=kc+1
            End Do
            iOff=iOff+nBas(iSym)**2
            jOff=jOff+nBas(iSym)
         End Do
*
         iOff=0
         kOff=0
         Do iSym=1,nSym
            ipDaa=ip_Da+kOff
            mAdCMOO=mAdCMO+iOff+nBas(iSym)*nIF(iSym)
            Call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nConstr(iSym),
     &                       1.0d0,Work(mAdCMOO),nBas(iSym),
     &                             Work(mAdCMOO),nBas(iSym),
     &                       0.0d0,Work(ipDaa),nBas(iSym))
            ipDbb=ip_Db+kOff
            mAdCMOO=mAdCMO_ab+iOff+nBas(iSym)*nIF(iSym)
            Call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nConstr(iSym),
     &                       1.0d0,Work(mAdCMOO),nBas(iSym),
     &                             Work(mAdCMOO),nBas(iSym),
     &                       0.0d0,Work(ipDbb),nBas(iSym))
            Do j=1,nBas(iSym)
               Do i=1,j-1
                  ji=j*(j-1)/2+i
                  iDaa=ipDaa-1+ji
                  Work(iDaa)=2.0d0*Work(iDaa)
                  iDbb=ipDbb-1+ji
                  Work(iDbb)=2.0d0*Work(iDbb)
               End Do
            End Do
            iOff=iOff+nBas(iSym)**2
            kOff=kOff+nBas(iSym)*(nBas(iSym)+1)/2
         End Do
*
         Call Get_Etwo_act(Work(ip_Da),Work(ip_Db),nBLT,nBas,nSym,Etwo)
*
         write(6,'(1X,A,F12.7,A)') ' Active-Active repulsion : ',
     &                             Etwo,'  a.u.'
         write(6,*) ' -------------------------------------------------'
         write(6,*)
         xNrm=sqrt(xNrm)

      End Do
*
      Call GetMem('Corb','Free','Real',ipCorb,MaxBas)
      Call GetMem('Match','Free','Inte',ipMatch,2*MxConstr)
      Call GetMem('DLT','Free','Real',ip_Da,2*nBLT)
      Call GetMem('CMOb','Free','Real',mAdCMO,2*nBB)
      Call GetMem('Occb','Free','Real',mAdOcc_ab,nBT)

      Return
      End

************************************************************************
*                                                                      *
************************************************************************
      Subroutine Get_Etwo_act(Dma,Dmb,nBDT,nBas,nSym,Etwo)

      Implicit Real*8 (a-h,o-z)
      Integer nBDT, nBas(8)
      Real*8  Dma(nBDT), Dmb(nBDT)
#include "WrkSpc.fh"
      Integer ALGO,NSCREEN
      Logical REORD,DECO
      Common /CHOSCF / REORD,DECO,dmpk,dFKmat,ALGO,NSCREEN
      Logical Estimate, Update, timings
      COMMON /CHOSCREEN/ Estimate,Update
      COMMON  /CHOTIME /timings
      Integer ipFLT(2), ipKLT(2), nForb(8,2), nIorb(8,2), ipPorb(2)
      Integer ipDm(2)
c      Real*8   Get_ExFac
c      External Get_ExFac
c      Character*16 KSDFT
*
      timings=.false.
      Estimate=.false.
      REORD=.false.
*
      Update=.true.
      DECO=.true.
      dmpk=1.0d0
      dFKmat=0.0d0
      ALGO=4
      NSCREEN=10
*
      nDMat=2
      nBB=0
      Do i=1,nSym
         nForb(i,1)=0
         nForb(i,2)=0
         nBB=nBB+nBas(i)**2
      End Do
c      Call Get_cArray('DFT functional',KSDFT,16)
c      ExFac=Get_ExFac(KSDFT)
c      FactXI=1.0d0*ExFac
      FactXI=1.0d0  ! always HF energy
      Call GetMem('PLTc','Allo','Real',ipPLT,nBDT)
      call dcopy_(nBDT,Dma,1,Work(ipPLT),1)
      Call daxpy_(nBDT,1.0d0,Dmb,1,Work(ipPLT),1)
*
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
         Call CD_InCore(Work(ipDai),nBas(i),Work(ipV),nBas(i),
     &                  nIorb(i,1),1.0d-12,irc)
         If (irc.ne.0) Then
            write(6,*) ' Alpha density. Sym= ',i,'   rc= ',irc
            Call Abend()
         EndIf
         ipV=ipPorb(2)+iOff
         ipDbi=ipDm(2)+iOff
         Call CD_InCore(Work(ipDbi),nBas(i),Work(ipV),nBas(i),
     &                  nIorb(i,2),1.0d-12,irc)
         If (irc.ne.0) Then
            write(6,*) ' Beta density. Sym= ',i,'   rc= ',irc
            Call Abend()
         EndIf
         iOff=iOff+nBas(i)**2
      End Do
*
      Call GetMem('FCNO','Allo','Real',ipFCNO,2*nBDT)
      Call FZero(Work(ipFCNO),2*nBDT)
      ipFLT(1)=ipFCNO
      ipFLT(2)=ipFLT(1)+nBDT
      Call GetMem('KLTc','Allo','Real',ipKLT(1),2*nBDT)
      Call FZero(Work(ipKLT(1)),2*nBDT)
      ipKLT(2)=ipKLT(1)+nBDT
*
      Call Cho_X_init(irc,ChFracMem)
      if (irc.ne.0) then
         Call WarningMessage(2,'Get_CNOs. Non-zero rc in Cho_X_init.')
         Call Abend
      endif

! BIGOT FIXME
      Call WarningMessage(2,
     &     'There is probably a bug here, ipPLT should have two '//
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
*
      Etwo = 0.5d0*(ddot_(nBDT,Dma,1,Work(ipFLT(1)),1)
     &     +        ddot_(nBDT,Dmb,1,Work(ipFLT(2)),1))
*
      Call GetMem('KLTc','Free','Real',ipKLT(1),2*nBDT)
      Call GetMem('FCNO','Free','Real',ipFCNO,2*nBDT)
      Call GetMem('DSQc','Free','Real',ipDm(1),2*nBB)
      Call GetMem('ChMc','Free','Real',ipPorb(1),2*nBB)
      Call GetMem('PLTc','Free','Real',ipPLT,nBDT)
*
      Return
      End
