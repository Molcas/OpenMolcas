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
      Subroutine Read_h0(nSize,nBas,ip_h0,Restart)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
      Character*8 Label
      Logical Restart
      Dimension nInts(1)
*                                                                      *
************************************************************************
*                                                                      *
      iOpt0=0
      iOpt1=1
*
      Call Allocate_Work(ip_h0,nsize)
*
      iComp=1
      iSyLbl=1
      Label='OneHam  '
      iRc=-1
      If (Restart) Then
         Call Get_dArray('LoProp H0',Work(ip_h0),nSize)
      Else
         Call iRdOne(iRc,iOpt1,Label,iComp,nInts,iSyLbl)
         If ( iRc.ne.0 ) Then
            Write (6,*) 'Read_h0: Error reading ONEINT'
            Write (6,'(A,A)') 'Label=',Label
            Call QTrace()
            Call Abend()
         End If
         If (nInts(1)+4.ne.nSize) Then
            Write (6,*) 'Local_Polar: nInts+4.ne.nSize',nInts(1)+4,nSize
            Call QTrace
            Call Abend()
         End If
         iRc=-1
         Call RdOne(iRc,iOpt0,Label,iComp,Work(ip_h0),iSyLbl)
         Call Put_dArray('LoProp H0',Work(ip_h0),nSize)
      End If
C     Call TriPrt('H0 ',' ',Work(ip_h0),nBas)
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nBas)
      End
