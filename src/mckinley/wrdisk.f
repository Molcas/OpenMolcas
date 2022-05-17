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
      SubRoutine WrDisk(rIn,nrIn,jDisp,iIrrep)
*
*      Sorry about this litle subroutine
*      The reason is just that right now I want to use AO for SCF
*      and MO for RASSCF, this will hopefully be changed, but if you
*      see this mess before that I apologize
*
      use Basis_Info, only: nBas
      use pso_stuff
      use Symmetry_Info, only: nIrrep, iOper
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "buffer.fh"
#include "etwas.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "disp.fh"
#include "disp2.fh"
#include "print.fh"
      Real*8 rIn(nrIn)
      Character*8 Label
      Integer ip(0:7),ip2(0:7),nA(0:7),ipCM(0:7)
      Real*8, Allocatable:: Act(:), InAct(:), Out(:), TempX(:), TempY(:)
*                                                                      *
************************************************************************
*                                                                      *
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      If (Show) Then
         Write (6,*)
         Write (6,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
         Write (6,*)
         Write (6,*) 'jDisp=',jDisp
      End If
      nin=0
      nIn2=0
      nna=0
      ipCC=1
      Do jIrrep=0,nIrrep-1
         kIrrep=NrOpr(iEOr(iOper(jIrrep),iOper(iIrrep)))
         If (kIrrep.lt.jIrrep) Then
            ip(jIrrep)=nIn
            nIn=nIN+nBas(kIrrep)*nBas(jIrrep)
          Else If (kIrrep.eq.jIrrep) Then
            ip(jIrrep)=nIn
            nIn=nIN+nBas(kIrrep)*(1+nBas(jIrrep))/2
          End If
          ip2(kIrrep)=nIn2
          nIn2=nIn2+nBas(kIrrep)*nBas(jIrrep)
          nA(jIrrep)=nnA
          nnA=nnA+nAsh(jIrrep)
          ipCM(jIrrep)=ipCC
          ipCC=ipCC+nBas(jIrrep)**2
#ifdef __INTEL_COMPILER
*
*       To avoid error in intel optimization -O3
*
          If (.False.) Write (6,*) ip(jIrrep)
#endif
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(Act,nIn2,Label='Act')
      Act(:)=Zero
      Call mma_allocate(InAct,nIn2,Label='InAct')
      InAct(:)=Zero
      Call mma_allocate(Out,nIn2,Label='Out')
      Out(:)=Zero
      Call mma_allocate(TempX,nIn2,Label='TempX')
      Call mma_allocate(TempY,nIn2,Label='TempY')
*                                                                      *
************************************************************************
*                                                                      *
*     Fock1
*
      If (Show) Then
         Write (6,*)
         Write (6,*) 'Fock1'
         Write (6,*)
      End If
      Do jIrrep=0,nIrrep-1
         kIrrep=NrOpr(iEOr(iOper(jIrrep),iOper(iIrrep)))
         If (nBAs(jIrrep).gt.0.and.nbas(kIrrep).gt.0) Then
            If (kIrrep.lt.jIrrep) Then
               Call DGEMM_('N','N',
     &                     nBas(jIrrep),nBas(kIrrep),nBas(kIrrep),
     &                One,rIn(ipDisp(jDisp)+ip(jIrrep)),nBas(jIrrep),
     &                     CMO(ipCM(kIrrep),1),nBas(kIrrep),
     &                Zero,TempY,nBas(jIrrep))
               Call DGEMM_('T','N',
     &                     nBas(jIrrep),nBas(kirrep),nBas(jIrrep),
     &                     One,CMO(ipCM(jIrrep),1),nBas(jIrrep),
     &                         TempY,nBas(jIrrep),
     &                     Zero,Act(1+ip2(kIrrep)),nBas(jIrrep))
            If (Show) Then
               Write (6,*)
               Write (6,*) 'ipDisp(jDisp),ip(jIrrep)=',
     &                      ipDisp(jDisp),ip(jIrrep)
               Call RecPrt('ipDisp',' ',
     &             rIn(ipDisp(jDisp)+ip(jIrrep)),
     &             nBas(jIrrep),nBas(kIrrep))
               Write (6,'(A,G20.10)') 'ipDisp:',
     &               DDot_(nBas(jIrrep)*nBas(kIrrep),
     &                    rIn(ipDisp(jDisp)+ip(jIrrep)),1,
     &                    rIn(ipDisp(jDisp)+ip(jIrrep)),1)
               Write (6,'(A,G20.10)') 'ipCM(kIrrep):',
     &               DDot_(nBas(kIrrep)*nBas(kIrrep),
     &                    CMO(ipCM(kIrrep),1),1,
     &                    CMO(ipCM(kIrrep),1),1)
               Write (6,'(A,G20.10)') 'ipCM(jIrrep):',
     &               DDot_(nBas(jIrrep)*nBas(jIrrep),
     &                    CMO(ipCM(jIrrep),1),1,
     &                    CMO(ipCM(jIrrep),1),1)
            End If
               Call DGetMO(Act(1+ip2(kIrrep)),Nbas(jIrrep),
     &                  nbas(jIrrep),nBas(kIrrep),
     &                  Act(1+ip2(jIrrep)),nBas(kIrrep))
            Else If (kIrrep.eq.jIrrep) Then
               Call Square(rIn(ipDisp(jDisp)+ip(jIrrep)),TempX,
     &                  1,nBas(kirrep),nBas(kirrep))
               Call DGEMM_('N','N',
     &                     nBas(jIrrep),nBas(kIrrep),nBas(kIrrep),
     &                     One,TempX,nBas(jIrrep),
     &                         CMO(ipCM(kIrrep),1),nBas(kIrrep),
     &                     Zero,TempY,nBas(jIrrep))
               Call DGEMM_('T','N',
     &                     nBas(jIrrep),nBas(kirrep),nBas(jIrrep),
     &                     One,CMO(ipCM(jIrrep),1),nBas(jIrrep),
     &                         TempY,nBas(jIrrep),
     &                     Zero,Act(1+ip2(jIrrep)),nBas(jIrrep))
            End If
            If (Show) Then
               Write (6,*) 'jIrrep,kIrrep=',jIrrep,kIrrep
               Write (6,'(A,G20.10)') 'Act:',
     &               DDot_(nIn2,Act,1,Act,1)
            End If
         End If
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*
      If (nMethod.eq.RASSCF) Then
*
*        Fock2
*
         If (Show) Then
            Write (6,*)
            Write (6,*) 'Fock2'
            Write (6,*)
         End If
         Do jIrrep=0,nIrrep-1
            kIrrep=NrOpr(iEOr(iOper(jIrrep),iOper(iIrrep)))
            If (nBas(jIrrep).gt.0.and.nBas(kIrrep).gt.0) Then
               If (kIrrep.eq.jIrrep) Then
                  Call Square(rIn(ipDisp2(jDisp)+ip(jIrrep)),TempX,
     &                        1,nBas(kirrep),nBas(kirrep))
                  Call DGEMM_('N','N',
     &                        nBas(jIrrep),nBas(kIrrep),nBas(kIrrep),
     &                        One,TempX,nBas(jIrrep),
     &                            CMO(ipCM(kIrrep),1),nBas(kIrrep),
     &                        Zero,TempY,nBas(jIrrep))
                  Call DGEMM_('T','N',
     &                     nBas(jIrrep),nBas(kirrep),nBas(jIrrep),
     &                     One,CMO(ipCM(jIrrep),1),nBas(jIrrep),
     &                         TempY,nBas(jIrrep),
     &                     Zero,InAct(1+ip2(jIrrep)),nBas(jIrrep))
                Else If (kirrep.lt.jirrep) Then
                   Call DGEMM_('N','N',
     &                nBas(jIrrep),nBas(kIrrep),nBas(kIrrep),
     &                One,rin(ipDisp2(jDisp)+ip(jIrrep)),nBas(jIrrep),
     &                    CMO(ipCM(kIrrep),1),nBas(kIrrep),
     &                Zero,TempY,nBas(jIrrep))
                   Call DGEMM_('T','N',
     &                    nBas(jIrrep),nBas(kirrep),nBas(jIrrep),
     &                    One,CMO(ipCM(jIrrep),1),nBas(jIrrep),
     &                        TempY,nBas(jIrrep),
     &                    Zero,InAct(1+ip2(kIrrep)),nBas(jIrrep))
                   Call DGetMO(InAct(1+ip2(kIrrep)),Nbas(jIrrep),
     &                         nBas(jIrrep),nBas(kIrrep),
     &                         InAct(1+ip2(jIrrep)),nBas(kIrrep))
                End If
                If (Show) Then
                   Write (6,*) 'jIrrep,kIrrep=',jIrrep,kIrrep
                   Write (6,'(A,G20.10)') 'InAct:',
     &                   DDot_(nIn2,InAct,1,InAct,1)
                End If
             End If
          End Do
*                                                                      *
************************************************************************
*                                                                      *
*        Fock Tot
*
         If (Show) Then
            Write (6,*)
            Write (6,*) 'Fock Tot'
            Write (6,*)
         End If
         iii=0
         Do jIrrep=0,nIrrep-1
            kIrrep=NrOpr(iEOr(iOper(jIrrep),iOper(iIrrep)))
*
            If (nBas(jIrrep)*nIsh(kIrrep).gt.0) Then
               Call DaXpY_(nIsh(kIrrep)*nBas(jIrrep),2.0d0,
     &                    Act(1+ip2(kIrrep)),1,
     &                    Out(1+ip2(kIrrep)),1)
               Call DaXpY_(nIsh(kIrrep)*nBas(jIrrep),2.0d0,
     &                    InAct(1+ip2(kIrrep)),1,
     &                    Out(1+ip2(kIrrep)),1)
            End If
*
            If (nBas(jIrrep).gt.0) Then
               Do jAsh=1,nAsh(kIrrep)
                  Do kAsh=1,nAsh(kIrrep)
                     rDe=  G1(iTri(nA(kIrrep)+jAsh,nA(kIrrep)+kAsh),1)
                     ipOut=1+ip2(kIrrep)+nIsh(kIrrep)*nBas(jIrrep)+
     &                     nBas(jIrrep)*(kAsh-1)
                     ipIn1=1+ip2(kIrrep)
     &                    +nBas(jIrrep)*(jAsh-1+nIsh(kIrrep))
                     Call DaXpY_(nBas(jIrrep),rde,InAct(ipIn1),1,
     &                                              Out(ipOut),1)
                  End Do
               End Do
            End If
*
            If (nBas(jIrrep)*nAsh(kIrrep).gt.0) Then
               ipOut=1+ip2(kIrrep)+nIsh(kIrrep)*nBas(jIrrep)
               If (Show) Then
                  Write (6,*) 'jIrrep,kIrrep=',jIrrep,kIrrep
                  Write (6,'(A,G20.10)') 'ipDisp3:',
     &                  DDot_(nBas(jIrrep)*nAsh(kIrrep),
     &                            rin(ipDisp3(jDisp)+iii),1,
     &                            rin(ipDisp3(jDisp)+iii),1)
               End If
               Call DGEMM_('T','N',
     &                     nBas(jIrrep),nAsh(kIrrep),nBas(jIrrep),
     &                     One,CMO(ipCM(jIrrep),1),nBas(jIrrep),
     &                         rin(ipDisp3(jDisp)+iii),nBas(jIrrep),
     &                     Zero,TempY,nBas(jIrrep))
               Call DaXpY_(nAsh(kIrrep)*nBas(jIrrep),One,
     &                    TempY,1,
     &                    Out(ipOut),1)
               iii=iii+nBas(jIrrep)*nAsh(kIrrep)
            End If
#ifdef __INTEL_COMPILER
            If (.False.) Write (6,*) kIrrep, iii
#endif
            If (Show) Then
               Write (6,*) 'jIrrep,kIrrep=',jIrrep,kIrrep
               Write (6,'(A,G20.10)') 'Out:',
     &               DDot_(nIn2,Out,1,Out,1)
            End If
         End Do
         If (Show) Write (6,*)
*                                                                      *
************************************************************************
*                                                                      *
         irc=-1
         iopt=0
         Label='TOTAL'
         Call dWrMck(irc,iopt,Label,jdisp,Out,2**iIrrep)
         If (iRc.ne.0) Then
            Write (6,*) 'WrDisk: Error writing to MCKINT'
            Write (6,'(A,A)') 'Label=',Label
            Call Abend()
         End If
         If (Show) Then
            Write (6,'(A,G20.10)') 'TOTAL:',
     &            DDot_(nIn2,Out,1,Out,1)
         End If
*
         irc=-1
         iopt=0
         Label='INACTIVE'
         Call dWrMck(irc,iopt,Label,jdisp,InAct,2**iIrrep)
         If (iRc.ne.0) Then
            Write (6,*) 'WrDisk: Error writing to MCKINT'
            Write (6,'(A,A)') 'Label=',Label
            Call Abend()
         End If
         If (Show) Then
            Write (6,'(A,G20.10)') 'INACTIVE:',
     &            DDot_(nIn2,InAct,1,InAct,1)
            Write (6,*)
         End If
*
         irc=-1
         iopt=0
         Label='MOPERT'
         nt=nna*(nna+1)/2
         nt=nt*(nt+1)/2
         Call dWrMck(irc,iopt,Label,jdisp,rIn(ipMO(jdisp,1)),2**iIrrep)
         If (iRc.ne.0) Then
            Write (6,*) 'WrDisk: Error writing to MCKINT'
            Write (6,'(A,A)') 'Label=',Label
            Call Abend()
         End If
*                                                                      *
************************************************************************
*                                                                      *
*         SCF case
*
      Else
*                                                                      *
************************************************************************
*                                                                      *
         irc=-1
         iopt=0
         Label='TOTAL'
         Call dWrMck(irc,iopt,Label,jdisp,Act,2**iIrrep)
         If (iRc.ne.0) Then
            Write (6,*) 'WrDisk: Error writing to MCKINT'
            Write (6,'(A,A)') 'Label=',Label
            Call Abend()
         End If
         If (Show) Then
            Write (6,'(A,G20.10)') 'TOTAL:',
     &            DDot_(nIn2,Act,1,Act,1)
         End If
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(TempY)
      Call mma_deallocate(TempX)
      Call mma_deallocate(Out)
      Call mma_deallocate(InAct)
      Call mma_deallocate(Act)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
