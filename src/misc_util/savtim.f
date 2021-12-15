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
      SubRoutine SavTim(iFld,TCPU,TWall)
      Implicit Real*8 (a-h,o-z)
#include "para_info.fh"
#include "timtra.fh"
#include "WrkSpc.fh"
*
      if(nfld_tim.eq.0) return
      If (iFld.gt.nfld_tim) Then
         Call WarningMessage(2,'SavTim: iFld.gt.nfld_tim')
         Write (6,*) 'iFld=',iFld
         Write (6,*) 'nFld_tim=',nFld_tim
         Call Abend()
      End If
      iad=iGATim+myrank*nFld_Tim*2+iFld-1
      Work(iad)=Work(iad)+TCPU
      Work(iad+nFld_Tim)=Work(iad+nFld_Tim)+TWall
      Return
      End
      SubRoutine SavStat(iFld,Value,op)
      Implicit Real*8 (a-h,o-z)
      character*(*) op
#include "para_info.fh"
#include "timtra.fh"
#include "WrkSpc.fh"
*
      if(nfld_stat.eq.0) return
      If (iFld.gt.nfld_Stat) Then
         Call WarningMessage(2,'SavStat: iFld.gt.nfld_stat')
         Write (6,*) 'iFld=',iFld
         Write (6,*) 'nFld_Stat=',nFld_Stat
         Call Abend()
      End If
      iad=iGAStat+myrank*nFld_stat+iFld-1
      if(op.eq.'+') then
        Work(iad)=Work(iad)+Value
      else if(op.eq.'-') then
        Work(iad)=Work(iad)-Value
      else if(op.eq.'=') then
        Work(iad)=Value
      end if
      Return
      End
