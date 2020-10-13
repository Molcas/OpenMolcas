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
      Subroutine Localize_LoProp_Drv(Ttot,Ttot_Inv,nBas,iCenter,iType,
     &                            nBas1,nBas2,nSym,nBasMax,ipP,Restart)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 Ttot(nBas2), Ttot_Inv(nBas2)
      Integer iCenter(nBas1), iType(nBas1), nBas(nSym)
      Character*8 Label
      Logical Found,Restart
      Dimension idum(1)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*---- Get the overlap matrix
*
      iOpt1=1
      iOpt0=0
      Label='Mltpl  0'
      iRc=-1
      iSyLbl=1
      If (Restart) Then
         Call Qpg_iArray('LoProp nInts',Found,nElem)
         Call Allocate_iWork(ip_restart,nElem)
         Call Get_iArray('LoProp nInts',iWork(ip_restart),nElem)
         nInts = iWork(ip_restart+0) - 4
         Call GetMem('Tmp','Allo','Real',ip_SSym,nInts+4)
         Call Qpg_dArray('LoProp Integrals',Found,nInts_tot)
         If (.Not. Found) Then
            Write(6,*) 'LoProp Integrals not available on the RunFile.'
            Call Abend()
         End If
         Call Allocate_Work(ip_all_ints,nInts_Tot)
         Call Get_dArray('LoProp Integrals',Work(ip_all_ints),nInts_tot)
         Call dCopy_(iWork(ip_restart+0),Work(ip_all_ints),1,
     &                                Work(ip_SSym),1)
         Call Get_iArray('LoProp iSyLbl',iWork(ip_restart),nElem)
         iSyLbl = iWork(ip_restart+0)
         Call Free_Work(ip_all_ints)
         Call Free_iWork(ip_restart)
      Else
         Call iRdOne(iRc,iOpt1,Label,1,idum,iSyLbl)
         If ( iRc.ne.0 ) Then
            Write (6,*) 'Polar: error reading length of mu!'
            Write (6,*) 'Mu=',0
            Call Abend()
         End If
         nInts=idum(1)
         Call GetMem('Tmp','Allo','Real',ip_SSym,nInts+4)
         Call RdOne(iRc,iOpt0,Label,1,Work(ip_SSym),iSyLbl)
         If ( iRc.ne.0 ) Then
            Write (6,*) 'Polar: error reading mu!'
            Write (6,*) 'Mu=',0
            Call Abend()
         End If
      End If
#ifdef _DEBUGPRINT_
      Call RecPrt('SSym',' ',Work(ip_SSym),1,nInts)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('SMatrix','Allo','Real',ip_Tmp,nBas2)
      If (nSym.eq.1) Then
         ipS = ip_Tmp
      Else
         Call GetMem('SMatrix','Allo','Real',ipS,nBas1**2)
      End If
*
      iOfft = ip_SSym
      iOffs = ip_Tmp
      Do iSym = 1, nSym
         If (nBas(iSym).eq.0) Go To 99
*
*        Now I square the overlap matrix because it has been stored as a
*        lower triangle
*
         Call Square(Work(iOfft),Work(iOffs),1,nBas(iSym),nBas(iSym))
*
         iOfft = iOfft + nBas(iSym)*(nBas(iSym)+1)/2
         iOffs = iOffs + nBas(iSym)**2
 99      Continue
      End Do
      Call Free_Work(ip_SSym)
*
      If (nSym.ne.1) Then
*
*------- Desymmetrize
*
         nScr=nBasMax*nBas1
         Call Allocate_Work(ipScr,nScr)
         Call FZero(Work(ipS),nBas1**2)
         Call Desymmetrize(Work(ip_Tmp),nBas2,Work(ipScr),nScr,
     &                     Work(ipS),nBas,nBas1,Work(ipP),nSym,
     &                     iSyLbl)
         Call Free_Work(ipScr)
         Call Free_Work(ip_Tmp)
*
      End If
#ifdef _DEBUGPRINT_
      Call RecPrt('Overlap matrix',' ',Work(ipS),nBas1,nBas1)
#endif
*
*
*---- Localize
*
      Call Localize_LoProp(Ttot,Ttot_Inv,nBas1,Work(ipS),iCenter,iType)
#ifdef _DEBUGPRINT_
      Call RecPrt('Ttot',' ',Ttot,nBas1,nBas1)
      Call RecPrt('Ttot_Inv',' ',Ttot_Inv,nBas1,nBas1)
      Call xSpot('Exit Localize_LoProp_Drv')
#endif
*
*                                                                      *
************************************************************************
*                                                                      *
      Call Free_Work(ipS)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
