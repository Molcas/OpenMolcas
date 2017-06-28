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
* Copyright (C) 2008, Roland Lindh                                     *
*               2017, Varinia Bernales, Roland Lindh                   *
************************************************************************
      Subroutine Mk_DMET_Shell(Info,nInfo,nBfn)
      Implicit Real*8 (A-H,O-Z)
      External Integral_RICD, Integral_RI_2
#include "itmax.fh"
#include "info.fh"
#include "SysDef.fh"
#include "real.fh"
#include "WrkSpc.fh"
*                                                                      *
************************************************************************
*                                                                      *
      iShll = Mx_Shll - 1
      If (nCnttp.eq.0) Then
         mdc = 0
      Else
         mdc = mdciCnttp(nCnttp) + nCntr(nCnttp)
      End If
      nCnttp = nCnttp + 1
      If (nCnttp.gt.Mxdbsc) Then
         Call WarningMessage(2,'Mk_DMET_Shell: Increase Mxdbsc')
         Call Abend()
      End If
      ipVal(nCnttp) = iShll + 1
      ipPrj(nCnttp) = -1
      ipSRO(nCnttp) = -1
      ipSOC(nCnttp) = -1
*
      Bsl(nCnttp)='.....DMET'
      Charge(nCnttp)=Zero
      iAtmNr(nCnttp)=1
      AuxCnttp(nCnttp)=.False.
      aCD_Thr(nCnttp)=One
      NoPairL(nCnttp)=.False.
      CrRep(nCnttp)=Zero
      pChrg(nCnttp)=.False.
      Fixed(nCnttp)=.False.
      nOpt(nCnttp) = 0
      nPrj_Shells(nCnttp) = 0
      nSRO_Shells(nCnttp) = 0
      nSOC_Shells(nCnttp) = 0
*
      nPrim=1
      nCntrc=nBfn
      nTot_Shells(nCnttp) = 1
      nVal_Shells(nCnttp) = 1
*
      iShll = iShll + 1
      AuxShell(iShll) = .False.
      iStrt = ipExp(iShll)
      nExp(iShll) = nPrim
      nBasis(iShll) = nCntrc
      nBasis_Cntrct(iShll) = nCntrc
      ipBk(iShll) = ip_Dummy
      ip_Occ(iShll) = ip_Dummy
      ipAkl(iShll) = ip_Dummy
      iEnd = iStrt + nPrim - 1
*     Exponent
      Work(iStrt)=Zero
*     Coefficients
      iStrt = iEnd + 1
      ipCff(iShll) = iStrt
      ipCff_Cntrct(iShll) = iStrt
      ipCff_Prim(iShll) = iStrt
      iEnd = iStrt + nPrim*nCntrc -1
      Call DCopy_(nPrim,One,0,Work(iStrt),1)
      call dcopy_(nPrim*nCntrc,Work(iStrt),1,Work(iEnd+1),1)
      iEnd = iEnd + nPrim*nCntrc
      If (iShll.lt.MxShll) ipExp(iShll+1) = iEnd + 1
*
      Transf(iShll)=.False.
      Prjct(iShll)=.False.
      iAngMx=Max(iAngMx,0)
*
*-----The coordinates
*
      ipCntr(nCnttp) = ipExp(iShll+1)
      nCnt = 1
      If (mdc+nCnt.gt.Mxdc) Then
         Call WarningMessage(2,'Mk_DMET_Shell: Increase Mxdb')
         Write (6,*) mdc,nCnt,Mxdc
         Call Abend()
      End If
      mdciCnttp(nCnttp)=mdc
      LblCnt(mdc+nCnt) = 'Origin'
      If (mdc+nCnt.gt.1) Call ChkLbl(LblCnt(mdc+nCnt),LblCnt,mdc+nCnt-1)
      iOff=ipCntr(nCnttp)+(nCnt-1)*3
      Work(iOff  )=Zero
      Work(iOff+1)=Zero
      Work(iOff+2)=Zero
      nCntr(nCnttp) = nCnt
      mdc = mdc + nCnt
      If (iShll.lt.MxShll) ipExp(iShll+1) = ipExp(iShll+1) + nCnt*3
*
*     Compute the number of elements stored in the dynamic memory so
*     far.
*
      nInfo = ipExp(iShll+1) - Info
*                                                                      *
************************************************************************
*                                                                      *
      Mx_Shll=iShll+1
      Mx_mdc=mdc
*
      iCnttp_Dummy=nCnttp
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
