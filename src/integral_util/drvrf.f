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
      Implicit Real*8 (A-H,O-Z)
      Real*8 h1(nh1), TwoHam(nh1), D(nh1)
#include "SysDef.fh"
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
#include "rctfld.fh"
#include "WrkSpc.fh"
      Logical First, Dff, NonEq
      Character*8 Label
      Real*8 RepNuc_Temp
      Save RepNuc_Temp
*
      iRout = 1
      iPrint = nPrint(iRout)
*                                                                      *
************************************************************************
*                                                                      *
      If (.Not.lRF) Return

*
      Call QEnter('DrvRF')
*
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
*
      Call Init_RctFld(NonEq,iCharge)
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('RFld1','Allo','Real',ipRFld1,nh1)
      Call GetMem('RFld2','Allo','Real',ipRFld2,nh1)
      Call FZero(Work(ipRFld2),nh1)
      If (First) RepNuc_Temp=RepNuc
*                                                                      *
************************************************************************
*                                                                      *
      If (lLangevin.or.iXPolType.gt.0) Then

*
*------- Reaction field a la polarizabilities and Langevin dipole
*        moments on a grid in a cavity in a dielectric medium.
*
         Call Langevin(h1,Work(ipRFld2),D,RepNuc,nh1,First,Dff)

      Else If (PCM) Then
*
*------- Reaction field a la PCM
*
*
         Call DrvPCM(h1,Work(ipRFld2),D,RepNuc,nh1,First,Dff,NonEq)

      Else If (lRFCav) Then
*
*------- Reaction field a la cavity in dielectric medium.
*
         Call RctFld(h1,Work(ipRFld2),D,RepNuc,nh1,First,Dff,NonEq)
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
      Call Get_Temp(Label,Work(ipRFld1),nh1)
      Call DaXpY_(nh1,-One,h1,1,Work(ipRFld1),1)
      Call DScal_(nh1,-One,Work(ipRFld1),1)
*     Add contribution to the TwoHam array
      Call DaXpY_(nh1,One,Work(ipRFld2),1,TwoHam,1)
*
*-----Store away information for perturbative calculations
*
      Call DaXpY_(nh1,One,Work(ipRFld2),1,Work(ipRFld1),1)
*
      ERFSelf=RepNuc-RepNuc_Temp
      EEE=DDot_(nh1,Work(ipRFld2),1,D,1)
      ERFSelf=ERFSelf-Half*EEE

*
      Call Put_dScalar('RF Self Energy',ERFSelf)
      Call Put_dArray('Reaction field',Work(ipRFld1),nh1)
*
      Call GetMem('RFld2','Free','Real',ipRFld2,nh1)
      Call GetMem('RFld1','Free','Real',ipRFld1,nh1)
*                                                                      *
************************************************************************
*                                                                      *
*---- Write out the matrix elements of the reaction field potential
*     to be exploited in subsequent codes like, CC, MP2, and CASPT2.
*
      Label='PotNucXX'
      Call Get_Temp(Label,RepNucXX,1)
      RepNuc_RF=RepNuc-RepNucXX
*
      Call GetMem('h1_RF','Allo','Real',ip_h1_RF,nh1+4)
      Call GetMem('h1_XX','Allo','Real',ip_h1_XX,nh1)
*
      Label='h1    XX'
      Call Get_Temp(Label,Work(ip_h1_XX),nh1)
      call dcopy_(nh1,h1,1,Work(ip_h1_RF),1)
      Call DaXpY_(nh1,-One,Work(ip_h1_XX),1,Work(ip_h1_RF),1)
      Call GetMem('h1_XX','Free','Real',ip_h1_XX,nh1)
*
      Work(ip_h1_RF+nh1+3)=RepNuc_RF
*
      iSyLbl=1
      iRc=-1
      iOpt=0
      iComp=1
      Label='OneHamRF'
      Call WrOne(iRc,iOpt,Label,iComp,Work(ip_h1_RF),iSyLbl)
*
      Call GetMem('h1_RF','Free','Real',ip_h1_RF,nh1+4)
*                                                                      *
************************************************************************
*                                                                      *
      Call Free_iSD()
      Call QExit('DrvRF')
      Return
      End
