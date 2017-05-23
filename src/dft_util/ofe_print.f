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
      Subroutine OFE_print(Energy_A)

      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
      Real*8 ReCharge(MxAtom)
      COMMON  / OFembed_R / Rep_EN,Func_AB,Func_A,Func_B,Energy_NAD,
     &                      V_Nuc_AB,V_Nuc_BA,V_emb
      COMMON  / OFembed_R2/ dFMD
      Character*16 NamRfil
      Character*10 Fmt
      Integer  Cho_X_GetTol
      External Cho_X_GetTol
*
      Call Get_iScalar('nSym',nSym)
      Call Get_iScalar('Unique atoms',nAtoms)
      Call Get_dArray('Effective nuclear Charge',ReCharge,nAtoms)
*
      Call Get_NameRun(NamRfil)
      Call NameRun('AUXRFIL')
      Call PotNuc_nad(nSym,nAtoms,ReCharge,ZRE_nad)
*
      Call Get_dEnergy(Energy_B)
      If (dFMD.gt.0.0) Call Get_dScalar('KSDFT energy',Ec_A)
*
      Call NameRun(NamRfil)
*
      iTol = Cho_X_GetTol(8)
      Call Add_Info('V_OFE',V_emb,1,iTol)
      Call Add_Info('V_NUC',V_Nuc_AB,1,iTol)
      Call Add_Info('E_NAD',Energy_NAD,1,iTol)
      Call Add_Info('RP_EN',Rep_EN,1,iTol)
*
      Fmt='(A,F19.10)'
      write(6,*)
      write(6,*)  '     -----------------------------------------------'
      write(6,*)  '      Orbital-Free Embedding Calculation : Results  '
      write(6,*)  '     -----------------------------------------------'
      write(6,Fmt)'        DFT energy  (A)    : ', Func_A
      write(6,Fmt)'        DFT energy  (B)    : ', Func_B
      write(6,Fmt)'        DFT energy (A+B)   : ', Func_AB
      write(6,*)
      write(6,Fmt)'        Nonelectr. Vemb    : ', V_emb ! for <A|Vq|A>
      write(6,*)
      write(6,Fmt)'        Energy (A)         : ', Energy_A
      write(6,Fmt)'        Energy (B)         : ', Energy_B
      write(6,Fmt)'        DFT energy (NAD)   : ', Energy_NAD
      write(6,Fmt)'        Vnuc(B)*rhoA       : ', V_Nuc_AB
      write(6,Fmt)'        Vnuc(A)*rhoB       : ', V_Nuc_BA
      write(6,Fmt)'        Electr. repulsion  : ', Rep_EN
      write(6,*)  '     -----------------------------------------------'
      write(6,Fmt)'       Nuclear rep. (A--B) : ', ZRE_nad
      write(6,Fmt)'       Energy (A+B)        : ', Energy_B+Energy_A
     &                                                     +Energy_NAD
     &                                                     +V_Nuc_AB
     &                                                     +V_Nuc_BA
     &                                                     +Rep_EN
     &                                                     +ZRE_nad
      If(dFMD.gt.0.0) write(6,Fmt)'       SCF restoring Ec(A) : ', Ec_A
      write(6,*)  '     -----------------------------------------------'
      write(6,*)
      write(6,*)
*
      Call Put_dScalar('NAD dft energy',Energy_NAD)
      Return
      End
************************************************************************
*
************************************************************************
      Subroutine Get_dEnergy(Energy)
      Implicit Real*8 (a-h,o-z)
      Real*8  Energy
      Logical Found_EAV

      Found_EAV=.false.
      Call Qpg_dScalar('Average energy',Found_EAV)

      If (Found_EAV) Then
         Call Get_dScalar('Average energy',Energy)
      Else
         Call Get_dScalar('Last energy',Energy)
      EndIf

      Return
      End
