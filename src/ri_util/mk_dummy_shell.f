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
************************************************************************
      Subroutine Mk_Dummy_Shell(Info,nInfo,DInf,nDInf)
************************************************************************
*                                                                      *
*     Add the final DUMMY SHELL!                                       *
*                                                                      *
* 2008 R. Lindh, Dept. of Theor. Chem., Univ. of Lund, Sweden          *
************************************************************************
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
      External Integral_RICD, Integral_RI_2
#include "itmax.fh"
#include "info.fh"
#include "SysDef.fh"
#include "real.fh"
#include "stdalloc.fh"
      Real*8 DInf(nDInf)
*                                                                      *
************************************************************************
*                                                                      *
      iShll = Mx_Shll - 1
      mdc = mdciCnttp(nCnttp) + dbsc(nCnttp)%nCntr
      nCnttp = nCnttp + 1
      If (nCnttp.gt.Mxdbsc) Then
         Call WarningMessage(2,'Mk_Dummy_Shell: Increase Mxdbsc')
         Call Abend()
      End If
      ipVal(nCnttp) = iShll + 1
      ipPrj(nCnttp) = -1
      ipSRO(nCnttp) = -1
      ipSOC(nCnttp) = -1
*
      Bsl(nCnttp)='.....RI_Dummy'
      Charge(nCnttp)=Zero
      iAtmNr(nCnttp)=1
      AuxCnttp(nCnttp)=.True.
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
      nCntrc=1
      nTot_Shells(nCnttp) = 1
      nVal_Shells(nCnttp) = 1
*
      iShll = iShll + 1
      AuxShell(iShll) = .True.
      iStrt = ipExp(iShll)
      nExp(iShll) = nPrim
      nBasis(iShll) = nCntrc
      nBasis_Cntrct(iShll) = nCntrc
      ipBk(iShll) = ip_Dummy
      ip_Occ(iShll) = ip_Dummy
      ipAkl(iShll) = ip_Dummy
      iEnd = iStrt + nPrim - 1
*     Exponent
      DInf(iStrt)=Zero
*     Coefficients
      iStrt = iEnd + 1
      ipCff(iShll) = iStrt
      ipCff_Cntrct(iShll) = iStrt
      ipCff_Prim(iShll) = iStrt
      iEnd = iStrt + nPrim*nCntrc -1
      DInf(iStrt) =One
      DInf(iEnd+1)=999999.0D0
      call dcopy_(nPrim*nCntrc,DInf(iStrt),1,DInf(iEnd+1),1)
      iEnd = iEnd + nPrim*nCntrc
      If (iShll.lt.MxShll) ipExp(iShll+1) = iEnd + 1
*
      Transf(iShll)=.False.
      Prjct(iShll)=.False.
*
*-----The coordinates
*
      nCnt = 1
      If (mdc+nCnt.gt.Mxdc) Then
         Call WarningMessage(2,'Mk_Dummy_Shell: Increase Mxdbsc')
         Call Abend()
      End If
      mdciCnttp(nCnttp)=mdc
      LblCnt(mdc+nCnt) = 'Origin'
      If (mdc+nCnt.gt.1) Call ChkLbl(LblCnt(mdc+nCnt),LblCnt,mdc+nCnt-1)
!     Call allocate(dbsc(nCnttp)%Coor(1:3,1:1))
      Call mma_allocate(dbsc(nCnttp)%Coor,3,1,Label='dbsc:C')
      dbsc(nCnttp)%Coor(1:3,1:1)=Zero
      dbsc(nCnttp)%nCntr = nCnt
      mdc = mdc + nCnt
*
*     Compute the number of elements stored in the dynamic memory so
*     far.
*
      nInfo = ipExp(iShll+1) - 1
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
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(Info)
*
      End
