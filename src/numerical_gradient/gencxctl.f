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
* Copyright (C) 2013, Roland Lindh                                     *
************************************************************************
      Subroutine genCxCTL(iStop,Cartesian)
      Implicit Real*8 (a-h,o-z)
************************************************************************
*                                                                      *
*     subroutine for automatic generation of coordinates for numerical *
*     differentiation based on the rlxctl.f routine.                   *
*                                                                      *
*     Author: R. Lindh, Uppsala University                             *
*             (c) 2013, November                                       *
************************************************************************
#include "info_slapaf.fh"
      Parameter(nLbl=10*MxAtom)
#include "real.fh"
#include "WrkSpc.fh"
#include "nadc.fh"
#include "weighting.fh"
#include "db.fh"
#include "print.fh"
      Logical Numerical, PrQ, Cartesian, Found, TSC
      Character*8 Lbl(nLbl)
*
      Lu=6
*
      Call QEnter('GenCxCTL')
      iRout = 32
      iPrint=nPrint(iRout)
*                                                                      *
************************************************************************
*                                                                      *
*-----Process the input
*
      LuSpool=21
      Call SpoolInp(LuSpool)
*
      Call RdCtl_Slapaf(iRow,iInt,nFix,LuSpool,.True.)
      Curvilinear=.FALSE.
      Cartesian  =.NOT.Curvilinear
      Numerical = .False. ! Just to define it, value is irrelevant here!
*
      Call Close_LuSpool(LuSpool)
*                                                                      *
************************************************************************
*                                                                      *
      jPrint=nPrint(iRout)
*
      If (nLbl.lt.nBVec) Then
         Call WarningMessage(2,'Error in GenCxCTL')
         Write (Lu,*)
         Write (Lu,*) '**********************'
         Write (Lu,*) ' ERROR: nLbl.lt.nBVec '
         Write (Lu,*) ' nLbl=',nLbl
         Write (Lu,*) ' nBVec=',nBVec
         Write (Lu,*) '**********************'
         Call Quit_OnUserError()
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute the Wilson B-matrix, these describe the transformations
*     between internal and Cartesian coordinates. Values of the
*     Internal coordinates are computed too.
*
      BSet=.True.
      HSet=.False.
      PrQ=.False.
*
      nFix=0
      nWndw=iter
      Call BMtrx(iRow,nBVec,ipB,nsAtom,mInt,ipqInt,Lbl,
     &           Work(ipCoor),nDimBC,Work(ipCM),AtomLbl,nSym,iOper,
     &           Smmtrc,Degen,BSet,HSet,iter,ipdqInt,ipShf,
     &           Work(ipGx),Work(ipCx),mTtAtm,iWork(ipANr),iOptH,
     &           User_Def,nStab,jStab,Curvilinear,Numerical,
     &           DDV_Schlegel,HWRS,Analytic_Hessian,iOptC,PrQ,mxdc,
     &           iCoSet,lOld,rHidden,nFix,nQQ,iRef,Redundant,nqInt,
     &           MaxItr,nWndw)
*
      nPrint(30) = nPrint(30)-1
*                                                                      *
************************************************************************
*                                                                      *
      Call Put_dArray('BMtrx',Work(ipB),3*nsAtom*mInt)
      Call Put_iScalar('No of Internal coordinates',mInt)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Generate list of coordinates for numerical differentiation
*                                                                      *
************************************************************************
*                                                                      *
*     Make lists for all Cartesian coordinates and the
*     corresponding displacement in internal coordinates
*
      Call Allocate_Work(ipCList,2*mInt*3*nsAtom)
      Call FZero(Work(ipCList),2*mInt*3*nsAtom)
      Call Allocate_Work(ipDList,mInt)
      Call FZero(Work(ipDList),mInt)
      Call Allocate_Work(ip_RefCoor,3*nsAtom)
*                                                                      *
************************************************************************
*                                                                      *
*     If in TS-search regime and without TSConstraints, remove
*     constraints (in slapaf this is done differently)
*
      Call qpg_iScalar('TS Search',Found)
      If (Found) Call Get_lScalar('TS Search',Found)
      Call f_Inquire('TSC',TSC)
      If (Found.and..not.TSC)
     &   Call Merge_Constraints('','','UDC',nLambda,iRow_c)
*                                                                      *
************************************************************************
*                                                                      *
*     Get the T-matrix
*
      Call Allocate_Work(ip_T,mInt**2)
      Call TMatrix(Work(ip_T))
      Call Put_iScalar('nLambda',nLambda)
      Call Put_dArray('T-matrix',Work(ip_T),mInt**2)
*     Call RecPrt('T-matrix',' ',Work(ip_T),mInt,mInt)
*                                                                      *
************************************************************************
*                                                                      *
*     Now start generating the displaced structure to be used in the
*     numerical differentiation.
*
*     Take a copy of the current structure - the reference
*     coordinates.
*
      call dcopy_(3*nsAtom,Work(ipCoor),1,Work(ip_RefCoor),1)
*
*     Loop over all displacements which are in the subspace in
*     which we like to minimize the energy. Hence, this will
*     eliminate naturally translational and rotational degrees
*     (3N-6) but also eliminate constrained degrees (3N-6-m)
*
      nDisp = mInt-nLambda
      Call GetMem('du','Allo','Real',ipdu,mInt)
*
*     Loop only over displacement which do not change the constraint.
*     Note that the constraints are placed first.
*
      Call GetMem(' DCF  ', 'Allo','Real',ipDCF, 3*nsAtom)
      Call GetMem(' dss  ', 'Allo','Real',ipdss, mInt)
      Call GetMem(' qTemp', 'Allo','Real',ipTmp, mInt)
      Do iDisp = 1+2*nLambda, 2*mInt
*
*        Get a virgin copy of the reference structure
*
         call dcopy_(3*nsAtom,Work(ip_RefCoor),1,Work(ipCoor),1)
*
*        Compute the effective index where to find the data
*
         Jter=Iter+1
*
*        Update the shift vector, in the space in which we split
*        the constraints and the reduced subspace.
*
         jpShf = ipShf + (Iter-1)*mInt
         Call FZero(Work(ipdu),mInt)
         Call FZero(Work(jpShf),mInt)
         jInter = (iDisp+1)/2
         If (Mod(iDisp,2).eq.0) Then
            Work(ipdu-1+jInter) = -Delta
         Else
            Work(ipdu-1+jInter) =  Delta
         End If
*        Call RecPrt('du',' ',Work(ipdu),mInt,1)
*
*        Transform displacement to the internal coordinate
*        space. This is a simple unitary transformation.
*
         Call DGEMM_('N','N',mInt,1,mInt,
     &               One,Work(ip_T),mInt,
     &                   Work(ipdu),mInt,
     &               Zero,Work(jpShf),mInt)
*        Call RecPrt('shf',' ',Work(jpShf),mInt,1)
*
*        Save the value of the displacement in the list.
*
         Work(ipDList-1+jInter) = Delta
*
*        Take a copy of the current values of the internal
*        coordinates.
*
         ip_From=ipqInt + (Iter-1)*mInt
         ip_To  =ipqInt + (Jter-1)*mInt
         call dcopy_(mInt,Work(ip_From),1,Work(ip_To),1)
*        Call RecPrt('Int_Ref',' ',Work(ip_To),1,mInt)
*
*        To the second set of coordinates add the shift.
*        This set of internal coordinates corresponds to
*        the set for which we like to get the Cartesian
*        coordinates.
*
         Call DaXpY_(mInt,One,Work(jpShf),1,Work(ip_To),1)
*        Call RecPrt('Int    ',' ',Work(ip_To),1,mInt)
*
*--------Transform the new internal coordinates to Cartesians
*
         PrQ=.False.
         BSet=.False.
         nWndw=Iter
         Call NewCar(Iter,nBVec,iRow,nsAtom,nDimBC,mInt,
     &               Work(ipCoor),ipB,Work(ipCM),
     &               Lbl,Work(ipShf),ipqInt,ipdqInt,
     &               Work(ipDCF),Work(ipdss),Work(ipTmp),Stop,
     &               AtomLbl,iOper,nSym,iSym,Smmtrc,
     &               Degen,Work(ipGx),Work(ipCx),mTtAtm,
     &               iWork(ipANr),iOptH,User_Def,
     &               nStab,jStab,Curvilinear,Numerical,
     &               DDV_Schlegel,HWRS, Analytic_Hessian,
     &               iOptC,PrQ,mxdc,iCoSet,rHidden,ipRef,
     &               Redundant,nqInt,MaxItr,nWndw)
*
         ip_To   = ipCList + (iDisp-1)*3*nsAtom
*
*        Move the new Cartesian coordinate to the list.
*
         call dcopy_(3*nsAtom,Work(ipCoor),1,Work(ip_To),1)
      End Do
      Call GetMem(' qTemp', 'Free','Real',ipTmp, mInt)
      Call GetMem(' dss  ', 'Free','Real',ipdss, mInt)
      Call GetMem(' DCF  ', 'Free','Real',ipDCF, 3*nsAtom)
*
      Call Free_Work(ipShf)
      Call Free_Work(ipdu)
      Call Free_Work(ip_T)
*
*     Call RecPrt('DList',' ',Work(ipDList),1,mInt)
*     Call RecPrt('CList',' ',Work(ipCList),3*nsAtom,2*mInt)
*
*     Save the lists on the runfile. To be used in the
*     numerical gradient module.
*
      Call Put_dArray('DList',Work(ipDList),mInt)
      Call Put_dArray('CList',Work(ipCList),2*mInt*3*nsAtom)
*
*     Deallocate temporary memory.
*
      Call Free_Work(ip_RefCoor)
      Call Free_Work(ipDList)
      Call Free_Work(ipCList)

*     Alaska only
      iStop=3
*
*     Done!
*
      If (ip_B.ne.ip_Dummy) Call Free_Work(ip_B)
      If (ip_dB.ne.ip_Dummy) Call Free_Work(ip_dB)
      If (ip_iB.ne.ip_iDummy) Call Free_iWork(ip_iB)
      If (ip_idB.ne.ip_iDummy) Call Free_iWork(ip_idB)
      If (ip_nqB.ne.ip_iDummy) Call Free_iWork(ip_nqB)
*
      If (ipNADC.ne.ip_Dummy) Call Free_Work(ipNADC)
      If (Ref_Geom) Call Free_Work(ipRef)
      If (Ref_Grad) Call Free_Work(ipGradRef)
      If (lRP)      Call Free_Work(ipR12)
      Call GetMem(' B ',    'Free','Real',ipB,   (nsAtom*3)**2)
*
      If (ipqInt.ne.ip_Dummy) Then
         Call GetMem('dqInt', 'Free','Real',ipdqInt, nqInt)
         Call GetMem('qInt', 'Free','Real',ipqInt, nqInt)
      End If
      Call GetMem('Relax', 'Free','Real',ipRlx, Lngth)
      Call GetMem('Grad',  'Free','Real',ipGrd, 3*nsAtom)
      Call GetMem('Coord', 'Free','Real',ipCoor,3*nsAtom)
      Call GetMem('Anr',   'Free','Inte',ipANr, nsAtom)
      Call GetMem('Charge','Free','Real',ipCM,  nsAtom)
*     The weights array length is actually the total number of atoms,
*     not just symmetry-unique, but that doesn't matter for deallocation
      Call GetMem('Weights','Free','Real',ipWeights,nsAtom)
*
*-----Terminate the calculations.
*
      Call QExit('GenCxCTL')
*
      Return
      End
