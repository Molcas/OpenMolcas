************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine DrvXV(h1,TwoHam,D,RepNuc,nh1,First,Dff,NonEq,
     &                 lRF,KSDFT,ExFac,iCharge,iSpin,
     &                 D1I,D1A,nD1,DFTFOCK,Do_DFT)
#ifdef _EFP_
      use EFP_Module
#endif
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "nq_info.fh"
#include "debug.fh"
      Real*8 h1(nh1), TwoHam(nh1), D(nh1,2)
      Real*8 D1I(nD1),D1A(nD1)
      Logical First, Dff, lRF, NonEq, Do_Grad, Do_DFT
      Dimension Grad(1),RN(1)
*
      Logical Do_OFemb,KEonly,OFE_first
      COMMON  / OFembed_L / Do_OFemb,KEonly,OFE_first
      Character*16  OFE_KSDFT
      COMMON  / OFembed_C / OFE_KSDFT
*
      Logical Do_ESPF
#ifdef _EFP_
      Logical EFP_On
#endif
      Character*(*) KSDFT
      Character*4 DFTFOCK
      Character*8 Label
*
      Debug=.False.
*                                                                      *
************************************************************************
*                                                                      *
*---- If first iterations save original h1 and RepNuc, the reaction
*     field and/or the dft section might later modify these. Any static
*     contributions to these are updated internally in the subroutine
*     rctfld.
*
*     We also save one set which are not to be modified at all.
*
      If (First) Then
         Label='PotNuc00'
         Call Put_Temp(Label,RepNuc,1)
         Label='h1_raw  '
         Call Put_Temp(Label,h1,nh1)
*
         Label='PotNucXX'
         Call Put_Temp(Label,RepNuc,1)
         Label='h1    XX'
         Call Put_Temp(Label,h1,nh1)
      End If

*
*---- Start always with the original h1 and RepNuc! Observe that they
*     now might contain static updates from the RF or the Langevin
*     codes through the action of the subroutine rctfld.
*
      Label='PotNuc00'
      Call Get_Temp(Label,RN,1)
      RepNuc=RN(1)
      Label='h1_raw  '
      Call Get_Temp(Label,h1,nh1)
cnf
*                                                                      *
************************************************************************
*                                                                      *
*                        ESPF section                                  *
*                                                                      *
************************************************************************
*                                                                      *
      Call DecideOnESPF(Do_ESPF)
      If (Do_ESPF) Call h1_espf(h1,RepNuc,nh1,First,Do_DFT)
cnf
*                                                                      *
************************************************************************
*                                                                      *
*              Reaction Field section                                  *
*                                                                      *
************************************************************************
*                                                                      *
      If (lRF) Call DrvRF(h1,TwoHam,D,RepNuc,nh1,First,Dff,NonEq,
     &                    iCharge)
*
************************************************************************
*                                                                      *
*                         DFT section                                  *
*                                                                      *
************************************************************************
*                                                                      *
      Do_Grad=.False.
      Grad=Zero
      nGrad=1
      If (KSDFT.ne.'SCF'.and.Do_DFT)
     &   Call DrvDFT(h1,TwoHam,D,RepNuc,nh1,First,Dff,lRF,KSDFT,ExFac,
     &               Do_Grad,Grad,nGrad,iSpin,D1I,D1A,nD1,DFTFOCK)
*
************************************************************************
*                                                                      *
*                         Orbital-Free Embedding section               *
*                                                                      *
************************************************************************
*                                                                      *
      If (Do_OFemb)
     &   Call DrvEMB(h1,D,RepNuc,nh1,OFE_KSDFT,ExFac,
     &               Do_Grad,Grad,nGrad,D1I,D1A,nD1,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _EFP_
      If (EFP_On())
     &   Call DrvEFP(First)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
