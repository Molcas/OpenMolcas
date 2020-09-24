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
      SubRoutine DrvRF(h1,TwoHam,D,RepNuc,nh1,First,Dff,NonEq,iCharge)
      use External_Centers, only: iXPolType
      Implicit Real*8 (A-H,O-Z)
      Real*8 h1(nh1), TwoHam(nh1), D(nh1)
#include "SysDef.fh"
#include "print.fh"
#include "real.fh"
#include "rctfld.fh"
#include "stdalloc.fh"
      Logical First, Dff, NonEq
      Character*8 Label
      Real*8 RepNuc_Temp
      Save RepNuc_Temp
      Dimension RepNucXX(1)
      Real*8, Allocatable:: RFld(:,:), h1_RF(:), h1_XX(:)
*
      iRout = 1
      iPrint = nPrint(iRout)
*                                                                      *
************************************************************************
*                                                                      *
      If (.Not.lRF) Return

*
*
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
*
      Call Init_RctFld(NonEq,iCharge)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(RFld,nh1,2,Label='RFld')
      RFld(:,2)=Zero
      If (First) RepNuc_Temp=RepNuc
*                                                                      *
************************************************************************
*                                                                      *
      If (lLangevin.or.iXPolType.gt.0) Then

*
*------- Reaction field a la polarizabilities and Langevin dipole
*        moments on a grid in a cavity in a dielectric medium.
*
         Call Langevin(h1,RFld(1,2),D,RepNuc,nh1,First,Dff)

      Else If (PCM) Then
*
*------- Reaction field a la PCM
*
*
         Call DrvPCM(h1,RFld(1,2),D,RepNuc,nh1,First,Dff,NonEq)

      Else If (lRFCav) Then
*
*------- Reaction field a la cavity in dielectric medium.
*
         Call RctFld(h1,RFld(1,2),D,RepNuc,nh1,First,Dff,NonEq)
*
*
      Else
         Call WarningMessage(2,
     &              'I do not know what reaction field type to use.')
         Call Abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Get the original one-electron hamiltonian and add the reaction
*     field contribution to it.
*
      Label='h1    XX'
      Call Get_Temp(Label,RFld(1,1),nh1)
      Call DaXpY_(nh1,-One,h1,1,RFld(1,1),1)
      Call DScal_(nh1,-One,RFld(1,1),1)
*     Add contribution to the TwoHam array
      Call DaXpY_(nh1,One,RFld(1,2),1,TwoHam,1)
*
*-----Store away information for perturbative calculations
*
      Call DaXpY_(nh1,One,RFld(1,2),1,RFld(1,1),1)
*
      ERFSelf=RepNuc-RepNuc_Temp
      EEE=DDot_(nh1,RFld(1,2),1,D,1)
      ERFSelf=ERFSelf-Half*EEE

*
      Call Put_dScalar('RF Self Energy',ERFSelf)
      Call Put_dArray('Reaction field',RFld(1,1),nh1)
*
      Call mma_deallocate(RFld)
*                                                                      *
************************************************************************
*                                                                      *
*---- Write out the matrix elements of the reaction field potential
*     to be exploited in subsequent codes like, CC, MP2, and CASPT2.
*
      Label='PotNucXX'
      Call Get_Temp(Label,RepNucXX,1)
      RepNuc_RF=RepNuc-RepNucXX(1)
*
      Call mma_allocate(h1_RF,nh1+4,Label='h1_RF')
      Call mma_allocate(h1_XX,nh1  ,Label='h1_XX')
*
      Label='h1    XX'
      Call Get_Temp(Label,h1_XX,nh1)
      call dcopy_(nh1,h1,1,h1_RF,1)
      Call DaXpY_(nh1,-One,h1_XX,1,h1_RF,1)
      Call mma_deallocate(h1_XX)
*
      h1_RF(nh1+3)=RepNuc_RF
*
      iSyLbl=1
      iRc=-1
      iOpt=0
      iComp=1
      Label='OneHamRF'
      Call WrOne(iRc,iOpt,Label,iComp,h1_RF,iSyLbl)
*
      Call mma_deallocate(h1_RF)
*                                                                      *
************************************************************************
*                                                                      *
      Call Free_iSD()
      Return
      End
