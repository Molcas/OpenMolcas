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
      Subroutine Put_QLast
      Implicit Real*8 (a-h,o-z)
#include "tlist.fh"
#include "WrkSpc.fh"
#include "real.fh"
*
      if(ipTskQ.eq.0) return
      iTCnSt_c=iTCnSt-1
c     Call XFlush(6)
c     Write (*,*)
c     Write (*,*) 'Put_QLast: '
c     Call RecPrt('TskQ',' ',Work(ipTskQ),2,nTasks)
c     Write (*,'(A,2F10.1)') 'Last indices of the task (QLast)=',QLast
c     Write (*,'(A,4I9)') 'ipTskQ,iTCnSt_c,nTasks,iTskCan=',
c    &                     ipTskQ,iTCnSt_c,nTasks,iTskCan
      Work(ipTskQ+(iTskCan-1)*2  )=QLast(1)
      Work(ipTskQ+(iTskCan-1)*2+1)=QLast(2)
*
      QLast(1)=Not_Used
      QLast(2)=Not_Used
c     Call RecPrt('TskQ',' ',Work(ipTskQ),2,nTasks)
c     Write (*,*)
c     Call XFlush(6)
*
      Return
      End
