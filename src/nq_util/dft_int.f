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
      Subroutine DFT_Int(Weights,mGrid,list_s,nlist_s,AOInt,nAOInt,
     &                   FckInt,nFckInt,SOTemp,nSOTemp,
     &                   TabAO,ipTabAO,nTabAO,dF_dRho,ndF_dRho,
     &                   nSym,iSpin,Flop,Rho,nRho,Scr,mScr,
     &                   Fact,ndc,mAO,TabAOMax,T_X,
     &                   list_bas,Functional_type,nAOMax)
************************************************************************
*                                                                      *
* Object: Front-end for compting DFT integrals                         *
*                                                                      *
*      Author:Roland Lindh, Dept. of Theor. Chem., Lund Unibersity,   ,*
*             SWEDEN. November 2008                                    *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      External Do_NInt1_d, Do_nInt1,
     &         Do_NInt2_d, Do_nInt2,
     &         Do_NInt3_d, Do_nInt3,
     &         Do_NInt4_d, Do_nInt4
#include "functional_types.fh"
      Integer Functional_type
      Real*8 Weights(mGrid), SOTemp(nSOTemp,iSpin), Fact(ndc**2),
     &       TabAO(nTabAO), Scr(mScr), Rho(nRho,mGrid),
     &       AOInt(nAOInt*nAOInt,iSpin), FckInt(nFckInt,iSpin),
     &       dF_dRho(ndF_dRho,mGrid), TabAOMax(nlist_s)
      Integer list_s(2,nlist_s), ipTabAO(nlist_s), list_bas(2,nlist_s)
*                                                                      *
************************************************************************
*                                                                      *
      If (Functional_type.eq.LDA_type) Then
         nFn=1
         nScr=iSpin*nFn*nAOMax
         Call DFT_IntX(Do_NInt1_d,Do_NInt1,
     &                 Weights,mGrid,list_s,nlist_s,AOInt,nAOInt,
     &                 FckInt,nFckInt,SOTemp,nSOTemp,
     &                 TabAO,ipTabAO,nTabAO,dF_dRho,ndF_dRho,
     &                 nSym,iSpin,Flop,Rho,nRho,Scr,nScr,
     &                 Fact,ndc,mAO,TabAOMax,T_X,
     &                 list_bas,nFn)
*                                                                      *
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.GGA_type) Then
         nFn=4
         nScr=iSpin*nFn*nAOMax
         Call DFT_IntX(Do_NInt2_d,Do_NInt2,
     &                 Weights,mGrid,list_s,nlist_s,AOInt,nAOInt,
     &                 FckInt,nFckInt,SOTemp,nSOTemp,
     &                 TabAO,ipTabAO,nTabAO,dF_dRho,ndF_dRho,
     &                 nSym,iSpin,Flop,Rho,nRho,Scr,nScr,
     &                 Fact,ndc,mAO,TabAOMax,T_X,
     &                 list_bas,nFn)
*                                                                      *
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.meta_GGA_type1) Then
         nFn=4
         nScr=iSpin*nFn*nAOMax
         Call DFT_IntX(Do_NInt4_d,Do_NInt4,
     &                 Weights,mGrid,list_s,nlist_s,AOInt,nAOInt,
     &                 FckInt,nFckInt,SOTemp,nSOTemp,
     &                 TabAO,ipTabAO,nTabAO,dF_dRho,ndF_dRho,
     &                 nSym,iSpin,Flop,Rho,nRho,Scr,nScr,
     &                 Fact,ndc,mAO,TabAOMax,T_X,
     &                 list_bas,nFn)
*                                                                      *
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.meta_GGA_type2) Then
         nFn=5
         nScr=iSpin*nFn*nAOMax
         Call DFT_IntX(Do_NInt3_d,Do_NInt3,
     &                 Weights,mGrid,list_s,nlist_s,AOInt,nAOInt,
     &                 FckInt,nFckInt,SOTemp,nSOTemp,
     &                 TabAO,ipTabAO,nTabAO,dF_dRho,ndF_dRho,
     &                 nSym,iSpin,Flop,Rho,nRho,Scr,nScr,
     &                 Fact,ndc,mAO,TabAOMax,T_X,
     &                 list_bas,nFn)
*                                                                      *
************************************************************************
*                                                                      *
      Else
         Write (6,*) 'DFT_Int: Illegal functional type!'
         Call Abend()
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
