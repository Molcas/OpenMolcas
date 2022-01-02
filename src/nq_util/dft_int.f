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
      Subroutine DFT_Int(Weights,mGrid,list_s,nlist_s,
     &                   FckInt,nFckInt,dF_dRho,ndF_dRho,
     &                   iSpin,Flop,Fact,ndc,mAO,
     &                   list_bas,Functional_type)
************************************************************************
*                                                                      *
* Object: Front-end for compting DFT integrals                         *
*                                                                      *
*      Author:Roland Lindh, Dept. of Theor. Chem., Lund Unibersity,   ,*
*             SWEDEN. November 2008                                    *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      External Do_NInt1_d, Do_nInt1X,
     &         Do_NInt2_d, Do_nInt2X,
     &         Do_NInt3_d, Do_nInt3X,
     &         Do_NInt4_d, Do_nInt4X
#include "functional_types.fh"
      Integer Functional_type
      Real*8 Weights(mGrid), Fact(ndc**2),
     &       FckInt(nFckInt,iSpin),dF_dRho(ndF_dRho,mGrid)
      Integer list_s(2,nlist_s), list_bas(2,nlist_s)
*                                                                      *
************************************************************************
*                                                                      *
      If (Functional_type.eq.LDA_type) Then
         nFn=1
         Call DFT_IntX(Do_NInt1_d,Do_nInt1X,
     &                 Weights,mGrid,list_s,nlist_s,
     &                 FckInt,nFckInt,
     &                 dF_dRho,ndF_dRho,
     &                 iSpin,Flop,
     &                 Fact,ndc,mAO,
     &                 list_bas,nFn)
*                                                                      *
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.GGA_type) Then
         nFn=4
         Call DFT_IntX(Do_NInt2_d,Do_nInt2X,
     &                 Weights,mGrid,list_s,nlist_s,
     &                 FckInt,nFckInt,
     &                 dF_dRho,ndF_dRho,
     &                 iSpin,Flop,
     &                 Fact,ndc,mAO,
     &                 list_bas,nFn)
*                                                                      *
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.meta_GGA_type1) Then
         nFn=4
         Call DFT_IntX(Do_NInt4_d,Do_nInt4X,
     &                 Weights,mGrid,list_s,nlist_s,
     &                 FckInt,nFckInt,
     &                 dF_dRho,ndF_dRho,
     &                 iSpin,Flop,
     &                 Fact,ndc,mAO,
     &                 list_bas,nFn)
*                                                                      *
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.meta_GGA_type2) Then
         nFn=5
         Call DFT_IntX(Do_NInt3_d,Do_nInt3X,
     &                 Weights,mGrid,list_s,nlist_s,
     &                 FckInt,nFckInt,
     &                 dF_dRho,ndF_dRho,
     &                 iSpin,Flop,
     &                 Fact,ndc,mAO,
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
