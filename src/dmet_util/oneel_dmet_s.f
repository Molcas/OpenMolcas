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
* Copyright (C) 1990,1991,1993,1999, Roland Lindh                      *
*               1990, IBM                                              *
************************************************************************
      SubRoutine OneEl_DMET_s(Kernel,KrnlMm,Label,ip,lOper,nComp,CCoor,
     &                      nOrdOp,rNuc,rHrmt,iChO,DMET_s,nBfn)
      use PrpPnt
      Implicit Real*8 (A-H,O-Z)
      External Kernel, KrnlMm
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "real.fh"
      Real*8 DMET_s(nBfn,nBfn)
      Real*8, Dimension(:), Allocatable :: Array
      Character Label*8
      Character L_Temp*8
      Real*8 CCoor(3,nComp), rNuc(nComp)
      Integer ip(nComp), lOper(nComp), iChO(nComp), iStabO(0:7)
      Integer iTwoj(0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/
*                                                                      *
************************************************************************
*                                                                      *
      Dum=rHrmt
      iDum=nOrdOp
      iDum=iCho(1)
      iRout = 112
      iPrint = nPrint(iRout)
      Call qEnter('OneEl')
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
      If (nIC.eq.0) Go To 999
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
         LenInt=n2Tri(lOper(iComp))
         LenTot=LenTot+LenInt+4
      End Do
      Call mma_allocate(Array,LenTot,label='Array')
      ip(1)=1
      call dcopy_(LenTot,Zero,0,Array(ip(1)),1)
      iadr=ip(1)
      do iComp = 1, nComp
         LenInt=n2Tri(lOper(iComp))
         ip(icomp)=iadr
         iadr=iadr+LenInt+4
*        Copy center of operator to work area.
         call dcopy_(3,Ccoor(1,iComp),1,Array(ip(iComp)+LenInt),1)
*        Copy nuclear contribution to work area.
         Array(ip(iComp)+LenInt+3) = rNuc(iComp)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute all SO integrals for all components of the operator.
*
      Do iBfn = 1, nBfn
         Do jBfn = 1, iBfn
            ijBfn = iBfn*(iBfn-1)/2 + jBfn - 1 + ip(1)
            Array(ijBfn)=DMET_s(iBfn,jBfn)
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*                    P O S T P R O C E S S I N G                       *
*                                                                      *
************************************************************************
*                                                                      *
      write(6,*) 'Overlap'
      Call PrMtrx(Label,lOper,nComp,ip,Array)
      If (iPrint.ge.10) Call PrMtrx(Label,lOper,nComp,ip,Array)
*                                                                      *
************************************************************************
*                                                                      *
*     Make a square sum on all the integrals for verification
      Call VrfMtrx(Label,lOper,nComp,ip,Array)
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute properties or write integrals to disc.
*
      Do iComp = 1, nComp
         iSmLbl = lOper(iComp)
*                                                                      *
************************************************************************
*                                                                      *
*---------- Write integrals to disc
*
            iOpt = 0
            iRC = -1
            L_Temp= 'Mltpl  0'
            iComp_=iComp
            Call WrOne(iRC,iOpt,L_Temp,iComp_,Array(ip(iComp)),iSmLbl)

            If (iRC.ne.0) then
               Call WarningMessage(2,
     &               ' *** Error in subroutine ONEEL ***,'//
     &               '     Abend in subroutine WrOne')
               Call Abend()
            End If
      End Do  ! iComp

*                                                                      *
************************************************************************
*                                                                      *
*---- Deallocate memory for integral
*
      Call mma_deallocate(Array)
*                                                                      *
************************************************************************
*VB                                                                      *
      write(6,*) "Array after deallocate"
 999  Continue
      Call qExit('OneEl')
      Return
      End
