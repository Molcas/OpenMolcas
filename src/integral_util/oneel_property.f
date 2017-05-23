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
      SubRoutine OneEl_Property(Kernel,KrnlMm,Label,ip,lOper,nComp,
     &                          CCoor,nOrdOp,rNuc,rHrmt,iChO,
     &                          D_tot,nDens,Property,Sig)
      Implicit Real*8 (A-H,O-Z)
      External Kernel, KrnlMm
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "real.fh"
      Character Label*8
      Real*8 CCoor(3,nComp), rNuc(nComp), Property(nComp), D_tot(nDens)
      Integer ip(nComp), lOper(nComp), iChO(nComp)
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 112
      iPrint = nPrint(iRout)
      Call qEnter('OneEl')
      If (rHrmt.ne.One) Then
         Call WarningMessage(2,'OneEl_Property: rHrmt.ne.One')
         Call Abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute the one-electron integrals
*
      Call OneEl_Integrals(Kernel,KrnlMm,Label,ip,lOper,nComp,
     &                     CCoor,nOrdOp,rHrmt,iChO)
*                                                                      *
************************************************************************
*                                                                      *
*                    P O S T P R O C E S S I N G                       *
*                                                                      *
************************************************************************
*                                                                      *
      If (iPrint.ge.10)    Call PrMtrx(Label,lOper,nComp,ip,Work)
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute properties
*
      LenTot=0
      Do iComp = 1, nComp
         iSmLbl = lOper(iComp)
*                                                                      *
************************************************************************
*                                                                      *
*--------Compute properties directly from integrals
*
         nInt=n2Tri(iSmLbl)
         LenTot = LenTot + nInt + 4
         If (nInt.ne.0) Then
            Call CmpInt(Work(ip(iComp)),nInt,nBas,nIrrep,iSmLbl)
            If (nInt.ne.nDens) Then
               Call WarningMessage(2,'OneEl_Property: nInt.ne.nDens')
               Write (6,*) 'nInt=',nInt
               Write (6,*) 'nDens',nDens
               Call Abend()
            End If
            Property(iComp)=rNuc(iComp)
     &                     -Sig*DDot_(nDens,D_Tot,1,Work(ip(iComp)),1)
         Else
            Property(iComp)=rNuc(iComp)
         End If
*
      End Do  ! iComp
*                                                                      *
************************************************************************
*                                                                      *
*---- Deallocate memory for integral
*
      Call GetMem(' ','FREE','REAL',ip(1),LenTot)
*
*                                                                      *
************************************************************************
*                                                                      *
      Call qExit('OneEl')
      Return
      End
      SubRoutine OneEl_Integrals(Kernel,KrnlMm,Label,ip,lOper,nComp,
     &                           CCoor,nOrdOp,rHrmt,iChO)
      Implicit Real*8 (A-H,O-Z)
      External Kernel, KrnlMm
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "real.fh"
      Character Label*8
      Real*8 CCoor(3,nComp)
      Integer ip(nComp), lOper(nComp), iChO(nComp), iStabO(0:7)
      Integer iTwoj(0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 112
      iPrint = nPrint(iRout)
*     Call qEnter('OneEl_I')
      If (iPrint.ge.19) Then
         Write (6,*) ' In OneEl: Label', Label
         Write (6,*) ' In OneEl: nComp'
         Write (6,'(1X,8I5)') nComp
         Write (6,*) ' In OneEl: lOper'
         Write (6,'(1X,8I5)') lOper
         Write (6,*) ' In OneEl: n2Tri'
         Do iComp = 1, nComp
            ip(iComp) = n2Tri(lOper(iComp))
         End Do
         Write (6,'(1X,8I5)') (ip(iComp),iComp=1,nComp)
         Call RecPrt(' CCoor',' ',CCoor,3,nComp)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute the number of blocks from each component of the operator
*     and the irreps it will span.
*
      nIC = 0
      llOper = 0
      Do iComp = 1, nComp
         llOper = iOr(llOper,lOper(iComp))
         Do iIrrep = 0, nIrrep-1
            If (iAnd(lOper(iComp),iTwoj(iIrrep)).ne.0) nIC = nIC + 1
         End Do
      End Do
      If (iPrint.ge.20) Write (6,*) ' nIC =',nIC
      If (nIC.eq.0) Then
         Call WarningMessage(2,'OneEl_Integrals: nIC.eq.0')
         Call Abend()
      End If
      Call SOS(iStabO,nStabO,llOper)
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate memory for symmetry adapted one electron integrals.
*     Will just store the unique elements, i.e. low triangular blocks
*     and lower triangular elements in the diagonal blocks.
*
      Call ICopy(nComp,-1,0,ip,1)
      LenTot=0
      Do iComp = 1, nComp
         ip(1)=1+LenTot
         LenInt=n2Tri(lOper(iComp))
         LenTot=LenTot+LenInt+4
      End Do
      Call GetMem(Label,'ALLO','REAL',ibase,LenTot)
      call dcopy_(LenTot,Zero,0,Work(ibase),1)
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute all SO integrals for all components of the operator.
*
      Call OneEl_(Kernel,KrnlMm,Label,ip,lOper,nComp,CCoor,
     &            nOrdOp,rHrmt,iChO,
     &            dum,1,dum,idum,0,0,
     &            iStabO,nStabO,nIC,
     &            Dum,1,0,Work(ibase),LenTot)
*
*
*     Fix pointer to be relative to Work(1)
*
      Do iComp = 1, nComp
         ip(iComp)=ip(iComp)-ip(1)+ibase
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Call qExit('OneEl_I')
      Return
      End
