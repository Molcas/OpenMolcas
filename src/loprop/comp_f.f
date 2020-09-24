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
      Subroutine Comp_F(h0,Ei,nBas,Delta_i,Energy,S,Refx,Originx)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
      Real*8 h0(*),Ei(*),S(*)
*
      nInts = nBas*(nBas+1)/2
      Call Allocate_Work(ip_h0,nInts+4)
      Call Comp_F_(h0,Ei,nBas,Delta_i,Energy,Work(ip_h0),S,Refx,
     &             Originx,nInts)
      Call Free_Work(ip_h0)
*
      Return
      End
      Subroutine Comp_F_(h0,Ei,nBas,Delta_i,Energy,h0_temp,S,Refx,
     &                   Originx,nInts)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "real.fh"
      Real*8 h0(nInts+4), h0_temp(nInts+4),
     &       Ei(nInts+4), S(nInts+4)
      Character*8 Method, Label
      Integer mBas(8)
*
      iOpt1=1
      iOpt2=2
*
      Call Get_cArray('Relax Method',Method,8)
      Call Allocate_Work(ipC,1)
      Call Get_iScalar('nSym',i)
      Call Get_iArray('nBas',mBas,i)
      nsize=nInts + 4
*
*     Add perturbations to h0
*
      call dcopy_(nSize,h0,1,h0_temp,1)
      Call DaXpY_(nSize,Delta_i,Ei,1,h0_temp,1)
      Call DaXpY_(nSize,Delta_i*(Originx-Refx),S ,1,h0_temp,1)
      Call Get_dScalar('PotNuc',PotNuc_Save)
      Call Put_dScalar('PotNuc',h0_temp(nInts+4))
*
*---- Write the one-electron hamiltonian.
*
C     Call TriPrt('H0 before wrone',' ',h0,nBas)
      iComp=1
      iSyLbl=1
      Label='OneHam  '
      iRc=-1
      Call WrOne(iRc,0,Label,iComp,h0_temp,iSyLbl)
C     Call TriPrt('H0_temp after wrone',' ',h0_temp,nBas)
*
      If (Method.eq.'RHF-SCF' .or.
     &    Method.eq.'UHF-SCF' .or.
     &    Method.eq.'KS-DFT' ) Then
*
         Call StartLight('scf')
         Call Disable_Spool()
         Call xml_open('module',' ',' ',0,'scf')
         Call SCF(ireturn)
         Call xml_close('module')
         If (iReturn.ne.0) Go To 99
         Call GetMem('PT2','Flush','Real',ipC,iDum)
*
      Else If (Method(1:5).eq.'MBPT2') Then
*
         Call StartLight('scf')
         Call Disable_Spool()
         Call xml_open('module',' ',' ',0,'scf')
         Call SCF(ireturn)
         Call xml_close('module')
         If (iReturn.ne.0) Go To 99
         Call GetMem('PT2','Flush','Real',ipC,iDum)
         Call StartLight('mbpt2')
         Call Disable_Spool()
         Call mp2_driver(ireturn)
         If (iReturn.ne.0) Go To 99
         Call GetMem('PT2','Flush','Real',ipC,iDum)
*
      Else If (Method.eq.'RASSCF' .or.
     &         Method.eq.'CASSCF' ) Then
*
         Call StartLight('rasscf')
         Call Disable_Spool()
         Call RASSCF(ireturn)
         If (iReturn.ne.0) Go To 99
         Call GetMem('PT2','Flush','Real',ipC,iDum)
*
      Else If (Method.eq.'CASPT2') Then
*
         Call StartLight('rasscf')
         Call Disable_Spool()
         Call RASSCF(ireturn)
         If (iReturn.ne.0) Go To 99
         Call GetMem('PT2','Flush','Real',ipC,iDum)
         Call StartLight('caspt2')
         Call Disable_Spool()
         Call CASPT2(ireturn)
         If (iReturn.ne.0) Go To 99
         Call GetMem('PT2','Flush','Real',ipC,iDum)
*
      Else
         Write (6,*) 'Method=',Method
         Write (6,*) ' Oups!'
         Call QTrace
         Call Abend()
      End If
*
C     Call Get_Energy(Energy)
      Call Get_dScalar('Last energy',Energy)
      Call Free_Work(ipC)
*
      Call WrOne(iRc,0,Label,iComp,h0,iSyLbl)
      Call Put_dScalar('PotNuc',PotNuc_Save)
      Return
*
 99   Continue
      Write (6,*)
      Write (6,*) 'Comp_f: Wave function calculation failed!'
      Write (6,*)
      Call QTrace()
      Call Abend()
* Avoid unused argument warnings
      If (.False.) Call Unused_integer(nBas)
*
      End
