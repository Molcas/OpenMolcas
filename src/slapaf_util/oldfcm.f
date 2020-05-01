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
      Subroutine OLDFCM(ipH,nQQ,nsAtom,iPrint,RunOld)
      Implicit Real*8 (a-h,o-z)
************************************************************************
*                                                                      *
*     Object : To read in a force constant matrix from another         *
*              interphase.                                             *
*                                                                      *
************************************************************************
#include "WrkSpc.fh"
      Character*8 Method
      Character*(*) RunOld
      Logical Found
*
*...  Prologue
      Call qEnter('OldFcm')
*
*...  Set runfile to be RUNOLD
      Call NameRun(RunOld)
*
*...  Get the method used in the old calculation
         Call Get_cArray('Relax Method',Method,8)
*
*...  Get the final energy obtained by the last calculation
*        Call Get_Energy(Energy)
         Call Get_dScalar('Last energy',Energy)
*
*...  Get the number of internal coordinates
         Call Get_iScalar('No of Internal coordinates',iInter)
         If ( iInter.le.0 ) Then
            Call WarningMessage(2,
     &                  'OldFCM: iInter.le.0')
            Write (6,*) 'iInter=',iInter
            Call Abend()
         End If
*
*...  Get the force constant matrix
         Call qpg_dArray('Hess',Found,nHess)
         If (.not.Found .or. nHess.eq.0) Then
            Call SysAbendmsg('OldFcm','Did not find:','Hess')
         End If
         Call GetMem('Hess','Allo','Real',ipRlx,nHess)
         Call get_dArray('Hess',Work(ipRlx),nHess)

         lHess = iInter**2
         If ( nHess.ne.lHess ) Then
            Call WarningMessage(2,'OldFCM: nHess.ne.lHess')
            Write (6,*) 'length,lHess=',length,lHess
            Call Abend()
         End If
*
*...  Reset runfile to be RUNFILE
      Call NameRun('RUNFILE')
*
*...  Echo the input information
         If (iPrint.ge.6) Then
         Write(6,*)
         Write(6,'(6X,A)')
     &   'SLAPAF has been supplied with an old force constant matrix.'
         Write(6,'(6X,3A)')
     &   'It is based on ',Method,' calculations.'
         Write(6,'(6X,A,F18.10)')
     &   'The final energy was',Energy
         End If
         If (iPrint.ge.99)
     &   Call RecPrt(' OldFcm',' ',Work(ipRlx),iInter,iInter)
*
         ipH = ipRlx
         nQQ = iINter
*
*...  Epilogue, end
      Call qExit('OldFcm')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nsAtom)
      End
