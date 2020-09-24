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
* Copyright (C) 2002, Laura Gagliardi                                  *
*               2002, Roland Lindh                                     *
************************************************************************
*
      Subroutine Polar(ireturn)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*      Author:                                                         *
*                                                                      *
*             Laura Gagliardi, Dipartimento di Chimica Fisica,         *
*             University of Palermo, ITALY. December 2002              *
*             Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN.                                         *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "Molcas.fh"
#include "timtra.fh"
      Parameter (nElem=(iTabMx*(iTabMx**2+6*iTabMx+11)+6)/6)
#include "WrkSpc.fh"
#include "real.fh"
      Real*8 Origin(3,0:iTabMx), CoC(3)
      Integer nBas(8), ip_mu(0:nElem-1), nOrb(8),ip_sq_mu(0:nElem-1),
     &        ip_D(0:6)
      Logical NoField, Standard, Utility, UserDen,PrintDen,SubtractDen
      Logical Restart, TDensity, XHole, Diffuse(3), Exist,  LIonize
      Character*(LENIN) LblCnt(MxAtom)
      Character*(LENIN4) LblCnt4(MxAtom)
      Character*12  Opt_Method
      Dimension dLimmo(2)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
      Utility = .False.
#ifdef _DEBUG_
#endif
      ireturn=99
*                                                                      *
************************************************************************
*                                                                      *
*     Prelude
*
      Call Init_LoProp(nSym,nBas,nOrb,CoC,nAtoms,ipC,ipQ_Nuc,ip_ANr,
     &                 ip_Type,ip_Center,nSize,nBas1,nBas2,nBasMax,ipP,
     &                 ipPInv)
*                                                                      *
************************************************************************
*                                                                      *
*     Read the input
*
*     NoField is defaulted to true if symmetry is implied.
*
      NoField = nSym.ne.1
      Call Readin_polar(NoField,Delta,MpProp_Level,Bond_Threshold,
     &                  iPlot,iPrint,Standard,Opt_Method,UserDen,
     &                  PrintDen,SubtractDen,SubScale,Restart,TDensity,
     &                  nStateI,nStateF,XHole,Diffuse,dLimmo,Thrs1,
     &                  Thrs2,nThrs,ThrsMul,Alpha,LIonize)
      Call InfoToMp(nSym,nBas,Energy_Without_FFPT,ip_Ene_Occ,nOcOb,
     &                  UserDen,Restart)
*                                                                      *
************************************************************************
*                                                                      *
*     Do the LoProp localization.
*
      Call GetMem('Ttot','Allo','Real',ip_Ttot,nBas1**2)
      Call GetMem('TtotInv','Allo','Real',ip_Ttot_Inv,nBas1**2)
*
      Call  Localize_LoProp_Drv(Work(ip_Ttot),Work(ip_Ttot_Inv),nBas,
     &                          iWork(ip_Center),iWork(ip_Type),
     &                          nBas1,nBas2,nSym,nBasMax,ipPInv,Restart)
*
      Call Free_iWork(ip_type)
*                                                                      *
************************************************************************
*                                                                      *
*     Read in the multipole moment integrals once and for all.
*
      Call Get_iScalar('Highest Mltpl',lMax)
      Write (6,'(A,I2)') ' Multipole moments will be processed up to'//
     &                   ' order ',lMax
      Write(6,*)
      mElem=(lMax*(lMax**2+6*lMax+11)+6)/6
*
      nTemp = nBas1**2
      Call Allocate_Work(ip_tmp,nTemp)
*
      Call Allocate_Work(ipMPq,mElem)
      Call Read_Multipole_Int(lMax,ip_sq_mu,nBas,ip_mu,
     &                        Work(ip_Ttot),Work(ip_tmp),Origin,
     &                        Work(ipMPq),mElem,nBas1,nBas2,nBasMax,
     &                        nTemp,nSym,ipPInv,Restart,Utility)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the 1-particle density matrix
*
      Dlt=-Delta
      iPert = 0
      Call Get_Density_Matrix(ip_D(0),nBas1,nBas2,nBasMax,nBas,nSym,ipP
     &                      ,UserDen,PrintDen,SubtractDen,SubScale,
     &                       Work(ipQ_Nuc),nAtoms,iPert,Restart,Utility
     &                      ,TDensity,nStateI,nStateF)
*
*   If computing local xhole-dipole moments. Should come after
*   get_density_matrix so modified densities are used.
      If(XHole) Call Compute_XHole_Int(nBas,nSym,ipXHole2,dMolExpec)

      If (.Not.NoField) Then
*        Read the one-electron hamiltonian.
         Call Read_h0(nSize,nBas(1),ip_h0,Restart)
         Do iPert = 1, 6
            iF = (iPert+1)/2
            Dlt = - Dlt
            If ((.NOT. Restart).and.(.not.UserDen)) Then
               Call Comp_F(Work(ip_h0),Work(ip_mu(iF)),nBas(1),Dlt, Ep,
     &                     Work(ip_mu(0)),CoC(iF),Origin(iF,1))
            End If
            Call Get_Density_Matrix(ip_D(iPert),nBas1,nBas2,nBasMax,
     &                  nBas,nSym,ipP,UserDen,PrintDen,SubtractDen,
     &              SubScale,Work(ipQ_Nuc),nAtoms,iPert,Restart,Utility,
     &              TDensity,nStateI,nStateF)
         End Do
         Call Free_Work(ip_h0)
      End If
      Do i = mElem, 1, -1
         Call Free_Work(ip_mu(i-1))
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     S T A T I C   P R O P E R T I E S                                *
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute localized multipole moments
*
      nPert = 2*3+1
      If (NoField) nPert=1
      nij= (nAtoms*(nAtoms+1)/2)
      nmu=nij*mElem*nPert
      Call Allocate_Work(ipMP,nmu)
      Call Allocate_Work(ip_sq_temp,nTemp)
      Call Allocate_Work(ip_EC,3*nij)
*
      Call Local_Properties(Work(ipC),nAtoms,ip_sq_mu,mElem,
     &                      Work(ip_sq_temp),Origin,
     &                      iWork(ip_center),Work(ip_Ttot_Inv),
     &                      Work(ip_tmp),nij,nPert,ip_D,
     &                      Work(ipMP),lMax,Work(ipMPq),CoC,Work(ip_EC),
     &                      iWork(ip_ANr),Standard,nBas1,nTemp,
     &                      Work(ipQ_Nuc),Bond_Threshold,
     &                      Utility,Opt_Method,iPlot,iPrint,nSym)
*
*---- If XHole integrals are available, localize them. Most unfortunate,
*     the local_properties routine is focused on multipole moments,
*     hence we rather write a new routine for Xhole, than significantly
*     edit the local_properties routine.
      If(XHole) then
        Call Allocate_Work(ipXHLoc2,nij)
        Call Local_Xhole(ipXHole2,dMolExpec,nAtoms,nBas1,nTemp
     &                  ,iWork(ip_center),Work(ip_Ttot)
     &                  ,Work(ip_Ttot_Inv),Work(ipC)
     &                  ,nij,Work(ip_EC),iWork(ip_ANr)
     &                  ,Bond_Threshold,iPrint,ipXHLoc2)
        Call Free_Work(ipXHole2)
      Else
        ipXHLoc2 = ip_Dummy
      Endif
*
*---- If the dear user has requested to get diffuse distributions
*     associated to the multipoles, go here.
      If(Diffuse(1)) then
        Call GetMem('ToPoint','Allo','Real',iTP,nAtoms)
        Call GetMem('NotToPoint','Allo','Real',ipMPp,nmu)
        call dcopy_(nmu,Work(ipMP),1,Work(ipMPp),1)
        Call CoreToPoint(nAtoms,ipMPp,iTP)
        LuYou=IsFreeUnit(81)
        Call OpnFl('DIFFPR',LuYou,Exist)
        Call Diff_MotherGoose(Diffuse,nAtoms,nBas1,ipMPp,ipC,nij,ip_EC
     &                       ,ip_ANr,ip_Ttot
     &                       ,ip_Ttot_Inv,lMax,iTP,dLimmo
     &                       ,Thrs1,Thrs2,nThrs,iPrint
     &                       ,ThrsMul,LuYou)
        Close(LuYou)
        Call GetMem('ToPoint','Free','Real',iTP,nAtoms)
        Call GetMem('NotToPoint','Free','Real',ipMPp,nmu)
      Endif
*
      Do i = mElem, 1, -1
         Call Free_Work(ip_sq_mu(i-1))
      End Do
      Call Free_Work(ip_Ttot)
      Call Free_Work(ip_Ttot_Inv)
      Call Free_Work(ip_sq_temp)
      Call Free_iWork(ip_center)
*                                                                      *
************************************************************************
*                                                                      *
      Call Allocate_Work(ipPol,6*nij)
      Call Allocate_Work(ipCpl,6*nij)
      Call Allocate_Work(ipCplT,6*nij)

      If (.Not.NoField) Then
*                                                                      *
************************************************************************
*                                                                      *
*        D Y N A M I C   P R O P E R T I E S                           *
*                                                                      *
************************************************************************
*                                                                      *
*        Compute the fluctuating charges
*
         Call Make_Fluctuating_Charges(nAtoms,iWork(ip_ANr),nij,nPert,
     &                                 Work(ipMP),mElem,Work(ip_EC),
     &                                 Alpha)
*                                                                      *
************************************************************************
*                                                                      *
*        Assemble the localized polarizabilities
*
         Call Dynamic_Properties(Work(ip_tmp),nAtoms,Work(ipMP),nij,
     &                           nPert,mElem,Delta,Work(ip_EC),
     &                           Work(ipPol),iWork(ip_ANr),
     &                           Bond_Threshold, Work(ipCpl),
     &                           Work(ipCplT))

*
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Print out the properties
*
      Call Get_cArray('LP_L',LblCnt4,(LENIN4)*nAtoms)
      Do i=1,nAtoms
       LblCnt(i)(1:LENIN)=LblCnt4(i)(1:LENIN)
      EndDo
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate arrays needed in Print_Local
*
      nDim = nij*mElem
      Call Allocate_Work(ipxMP,nDim)
      Call Allocate_Work(ipxxMP,nDim)
      Call Allocate_Work(ipnxMP,nDim)
      Call Print_Local(Work(ipMP),nij,mElem,Work(ipC),nAtoms,CoC,
     &                 Work(ipQ_Nuc),lMax,LblCnt,Work(ipMPq),
     &                 Work(ip_EC),Work(ipPol),NoField,Work(ip_Tmp),
     &                 Work(ipxMP),Work(ipxxMP), Work(ipnxMP),
     &                 iWork(ip_ANr),nOcOb,Energy_Without_FFPT,
     &                 ip_Ene_Occ,MpProp_Level,Bond_Threshold,
     &                 XHole,Work(ipXHLoc2),dMolExpec,
     &                 Work(ipCpl), Work(ipCplT), LIonize)

      Call Free_Work(ipxMP)
      Call Free_Work(ipxxMP)
      Call Free_Work(ipnxMP)
*                                                                      *
************************************************************************
*                                                                      *
      Call Free_Work(ipPol)
      Call Free_Work(ipCpl)
      Call Free_Work(ipCplT)

      Call Free_Work(ipQ_Nuc)
      Call Free_Work(ipMPq)
      Call Free_Work(ip_EC)
      Call Free_Work(ipMP)
      Call Free_Work(ip_Tmp)
      Call Free_iWork(ip_ANr)
      Call Free_Work(ipC)
      If(Xhole)Call Free_Work(ipXHLoc2)
      If (nSym.ne.1) Then
         Call Free_Work(ipP)
         Call Free_Work(ipPInv)
      End If
*
* Set flag on runfile to identify that it can be used to restart a
* LoProp calculation in the future.
*
      If (.Not. Restart) Then
         LoProp_Mode = 1
         If (NoField) Then
            LoProp_Mode = 2
         End If
         Call Put_iScalar('LoProp Restart',LoProp_Mode)
      End If
*
*                                                                      *
************************************************************************
*                                                                      *
*     Cleanup so that finish will not scream.
*
      nfld_tim =0
      nfld_stat=0
*                                                                      *
************************************************************************
*                                                                      *
      ireturn=0
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
#endif
      Return
      End
