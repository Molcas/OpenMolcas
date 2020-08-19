
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
      Subroutine Mk_Dummy_Shell()
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
      dbsc(nCnttp)%Aux=.True.
      aCD_Thr(nCnttp)=One
      NoPairL(nCnttp)=.False.
      CrRep(nCnttp)=Zero
      pChrg(nCnttp)=.False.
      Fixed(nCnttp)=.False.
      nOpt(nCnttp) = 0
      nPrj_Shells(nCnttp) = 0
      nSRO_Shells(nCnttp) = 0
      nSOC_Shells(nCnttp) = 0
      dbsc(nCnttp)%Parent_iCnttp = 0
*
      nPrim=1
      nCntrc=1
      nTot_Shells(nCnttp) = 1
      nVal_Shells(nCnttp) = 1
*
      iShll = iShll + 1
      Shells(iShll)%Aux = .True.
      Call mma_allocate(Shells(iShll)%Exp,nPrim,Label='ExpDummy')
      Shells(iShll)%nExp=nPrim
      Shells(iShll)%nBasis=nCntrc
      Shells(iShll)%nBasis_c = nCntrc
*     Exponent
      Shells(iShll)%Exp(1)=Zero
*     Coefficients
      Call mma_allocate(Shells(iShll)%Cff_c,nPrim,nCntrc,2,
     &                  Label='Cff_c')
      Call mma_allocate(Shells(iShll)%pCff,nPrim,nCntrc,
     &                  Label='pCff')
      Call mma_allocate(Shells(iShll)%Cff_p,nPrim,nPrim ,2,
     &                  Label='Cff_p')
      Shells(iShll)%Cff_c(1,1,1)=One
      Shells(iShll)%Cff_c(1,1,2)=One
      Shells(iShll)%pCff(:,:) = Shells(iShll)%Cff_c(:,:,1)
*
      Shells(iShll)%Transf=.False.
      Shells(iShll)%Prjct =.False.
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
      Call mma_allocate(dbsc(nCnttp)%Coor,3,1,Label='dbsc:C')
      dbsc(nCnttp)%Coor(1:3,1:1)=Zero
      dbsc(nCnttp)%nCntr = nCnt
      mdc = mdc + nCnt
*                                                                      *
************************************************************************
*                                                                      *
      Mx_Shll=iShll+1
      Max_Shells=Mx_Shll
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
