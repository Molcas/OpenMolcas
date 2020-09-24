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
      Subroutine LoProp(ireturn)
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
*             Roland Lindh, Department of Chemical Physics,            *
*             University of Lund, SWEDEN.                              *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "Molcas.fh"
      Parameter (nElem=(iTabMx*(iTabMx**2+6*iTabMx+11)+6)/6)
#include "WrkSpc.fh"
#include "real.fh"
      Real*8 Origin(3,0:iTabMx), CoC(3)
      Integer nBas(8), ip_mu(0:nElem-1), nOrb(8),ip_sq_mu(0:nElem-1),
     &        ip_D(0:6)
      Logical NoField, Standard, Utility, UserDen,PrintDen,SubtractDen
      Logical Restart, TDensity,lSave, Reduce_Prt
      External Reduce_Prt
      Character*(LENIN4) LblCnt(MxAtom)
      Character*12  Opt_Method
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
      Utility = .True.
      Utility = .False.
#ifdef _DEBUGPRINT_
         Call QEnter('LoProp')
#endif
      lSave = ireturn.eq.0
      ireturn=99

      iPL=iPrintLevel(-1)
      If (Reduce_Prt().and.iPL.lt.3) iPL=0

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
      NoField = .True.
      Standard=.True.
      UserDen=.False.
      PrintDen=.False.
      SubtractDen=.False.
      Restart=.False.
      TDensity=.false.
      nStateI=1
      nStateF=1
      Delta=0.001D0
      MpProp_Level=0
      Bond_Threshold=-1.0D0
      iPlot = 0
      iPrint = 0
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
      lMax=0   ! do only charges
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
      Call Free_Work(ip_Ttot)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the 1-particle density matrix
*
      Dlt=-Delta
      iPert = 0
      Call Get_Density_Matrix(ip_D(0),nBas1,nBas2,nBasMax,nBas,nSym,
     &                       ipP,UserDen,PrintDen,SubtractDen,SubScale,
     &                       Work(ipQ_Nuc),nAtoms,iPert,Restart,Utility,
     &                       TDensity,nStateI,nStateF)
*
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
      If (iPL.ge.2) Then
         Write(6,*)
         Call CollapseOutput(1,'   Static properties:')
         Write(6,'(3X,A)')     '   ------------------'
         Write(6,*)
      End If
      Call Local_Properties(Work(ipC),nAtoms,ip_sq_mu,mElem,
     &                      Work(ip_sq_temp),Origin,
     &                      iWork(ip_center),Work(ip_Ttot_Inv),
     &                      Work(ip_tmp),nij,nPert,ip_D,
     &                      Work(ipMP),lMax,Work(ipMPq),CoC,Work(ip_EC),
     &                      iWork(ip_ANr),Standard,nBas1,nTemp,
     &                      Work(ipQ_Nuc),Bond_Threshold,
     &                      Utility,Opt_Method,iPlot,iPrint,nSym)
*
      Do i = mElem, 1, -1
         Call Free_Work(ip_sq_mu(i-1))
      End Do
      Call Free_Work(ip_Ttot_Inv)
      Call Free_Work(ip_sq_temp)
      Call Free_iWork(ip_center)
*                                                                      *
************************************************************************
*                                                                      *
      Call Allocate_Work(ipPol,6*nij)
*                                                                      *
************************************************************************
*                                                                      *
*     Print out the properties
*
      Call Get_cArray('LP_L',LblCnt,(LENIN4)*nAtoms)
      Call LoProp_Print(Work(ipMP),nij,nElem,nAtoms,Work(ipQ_Nuc),
     &                       LblCnt,lSave)
      If (iPL.ge.2) Then
         Call CollapseOutput(0,'   Static properties:')
         Write(6,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call Free_Work(ipPol)
      Call Free_Work(ipQ_Nuc)
      Call Free_Work(ipMPq)
      Call Free_Work(ip_EC)
      Call Free_Work(ipMP)
      Call Free_Work(ip_Tmp)
      Call Free_iWork(ip_ANr)
      Call Free_Work(ipC)
      If (nSym.ne.1) Then
         Call Free_Work(ipP)
         Call Free_Work(ipPInv)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      ireturn=0
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
         Call QExit('LoProp')
#endif
      Return
      End
