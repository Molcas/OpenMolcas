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
      Subroutine espf_mltp (natom,MltOrd,nMult,nGrdPt,ipTTT,ipMltp,
     &                      ipGrid,ipIsMM,ipExt,iPL)
      Implicit Real*8 (A-H,O-Z)
*
*     Compute the expected values of the ESPF operators
*     i.e. charges or charges+dipoles
*
#include "espf.fh"
*
      Character*(LENIN) CName(MxAtom)
      Save Axis
      Character*3 Axis(3)
      Data Axis/' x ',' y ',' z '/
*
      Call QEnter('espf_mltp')
*
      If (iPL.ge.5) Then
         Write(6,*) ' In espf_mltp:',MltOrd,nMult,nGrdPt,ipTTT,ipMltp,
     &               ipGrid,ipIsMM
         Call RecPrt('TTT',' ',Work(ipTTT),nGrdPt,nMult)
      End If
      Call GetMem('Nuclear charge','Allo','Real',ip_Charge,natom)
      Call Get_Nuc_Charge_All(Work(ip_Charge),natom)
      iMult = 0
      Do iAtom = 0, natom-1
         If (iWork(ipIsMM+iAtom).eq.0) Then
            Work(ipMltp+iMult) = Work(ip_Charge+iAtom)
            If (MltOrd.ne.1) then
               Do iOff = 1, MltOrd-1
                  Work(ipMltp+iMult+iOff) = Zero
               EndDo
            End If
            iMult = iMult + MltOrd
         End If
      EndDo
      Call GetMem('Nuclear charge','Free','Real',ip_Charge,natom)
      opnuc = Dum
      ncmp = 1
      iAddPot = -2
      Call GetMem('dESPF2','Allo','Real',ipD2,nGrdPt)
      Call DrvPot(Work(ipGrid),opnuc,ncmp,Work(ipD2),nGrdPt,iAddPot)
      If (iPL.ge.5) Call RecPrt('PV',' ',Work(ipD2),nGrdPt,1)
      Do jMlt = 0, nMult-1
         Do kPnt = 0, nGrdPt-1
            Work(ipMltp+jMlt) = Work(ipMltp+jMlt)
     &        +Work(ipD2+kPnt)*Work(ipTTT+jMlt*nGrdPt+kPnt)
         EndDo
      EndDo
      Call GetMem('dESPF2','Free','Real',ipD2,nGrdPt)
*
      If (iPL.ge.3) Then
         Write(6,'(/,A,/)')
     &           '      Expectation values of the ESPF operators:'
         Call GetMem('ElecInt','Allo','Real',ipEI,natom)
         Call Get_CArray('Unique Atom Names',CName,LENIN*natom)
         SumOfChg = Zero
         jMlt = 0
         TotElecInt = Zero
         Do iAtom = 0, natom-1
            Work(ipEI+iAtom) = Zero
            iCur = ipExt+iAtom*MxExtPotComp
            If (iWork(ipIsMM+iAtom).eq.1) goto 10
            Do kOrd = 0, MltOrd-1
               If (kOrd.eq.0) Then
                  Write(6,1000) CName(iAtom+1),Work(ipMltp+jMlt)
                  SumOfChg = SumOfChg + Work(ipMltp+jMlt)
               Else
                  Write(6,1001) Axis(kOrd),Work(ipMltp+jMlt+kOrd)
               End If
               Work(ipEI+iAtom) = Work(ipEI+iAtom)
     &                          + Work(ipMltp+jMlt+kOrd)*Work(iCur+kOrd)
            EndDo
            jMlt = jMlt + MltOrd
            TotElecInt = TotElecInt + Work(ipEI+iAtom)
10          Continue
         EndDo
         Write (6,1002) SumOfChg
         Write (6,1003) TotElecInt
         Do iAtom = 0, natom-1
            If (iWork(ipIsMM+iAtom) .eq. 0)
     &                   Write(6,1004) CName(iAtom+1),Work(ipEI+iAtom)
         EndDo
         Write(6,'(A)')
         Call GetMem('ElecInt','Free','Real',ipEI,natom)
1000  Format('        Charge on ',A,'      = ',F10.4)
1001  Format('        + Dipole component ',A3,'= ',F10.4)
1002  Format(/,'      Total ESPF charge     = ',F10.4,/)
1003  Format(/,'      Total ESPF QM/MM interaction energy = ',F10.6,/)
1004  Format('        ',A,' individual contribution =',F10.6)
      End If
*
      Call QExit('espf_mltp')
      Return
      End
