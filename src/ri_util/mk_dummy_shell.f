
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
      use Center_Info
      use Sizes_of_Seward, only: S
      Implicit Real*8 (A-H,O-Z)
      External Integral_RICD, Integral_RI_2
#include "Molcas.fh"
#include "itmax.fh"
#include "info.fh"
#include "SysDef.fh"
#include "real.fh"
#include "stdalloc.fh"
*                                                                      *
************************************************************************
*                                                                      *
      iShll = S%Mx_Shll - 1
      mdc = dbsc(nCnttp)%mdci + dbsc(nCnttp)%nCntr
      nCnttp = nCnttp + 1
      If (nCnttp.gt.Mxdbsc) Then
         Call WarningMessage(2,'Mk_Dummy_Shell: Increase Mxdbsc')
         Call Abend()
      End If
      dbsc(nCnttp)%iVal = iShll + 1
      dbsc(nCnttp)%nVal = 1
      dbsc(nCnttp)%nShells = dbsc(nCnttp)%nVal
*
      dbsc(nCnttp)%Bsl='.....RI_Dummy'
      dbsc(nCnttp)%AtmNr=1
      dbsc(nCnttp)%Aux=.True.
*
      nPrim=1
      nCntrc=1
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
      n_dc=max(mdc+nCnt,n_dc)
      If (mdc+nCnt.gt.MxAtom) Then
         Call WarningMessage(2,'Mk_Dummy_Shell: Increase MxAtom')
         Call Abend()
      End If
      dbsc(nCnttp)%mdci=mdc
      dc(mdc+nCnt)%LblCnt = 'Origin'
      If (mdc+nCnt.gt.1) Call Chk_LblCnt(dc(mdc+nCnt)%LblCnt,mdc+nCnt-1)
      Call mma_allocate(dbsc(nCnttp)%Coor_Hidden,3,1,Label='dbsc:C')
      dbsc(nCnttp)%Coor => dbsc(nCnttp)%Coor_Hidden(:,:)
      dbsc(nCnttp)%Coor(1:3,1:1)=Zero
      dbsc(nCnttp)%nCntr = nCnt
      mdc = mdc + nCnt
*                                                                      *
************************************************************************
*                                                                      *
      S%Mx_Shll=iShll+1
      Max_Shells=S%Mx_Shll
      S%Mx_mdc=mdc
*
      If (iCnttp_Dummy.ne.0) Then
         Write (6,*) 'Mk_dummy_shell: iCnttp_Dummy'
         Call Abend()
      End If
      iCnttp_Dummy=nCnttp
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
