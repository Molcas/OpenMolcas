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
      use TList_Mod
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
*
      if(.NOT.Allocated(TskQ)) return
c     Call RecPrt('TskQ',' ',TskQ,2,nTasks)
c     Write (*,'(A,2F10.1)') 'Last indices of the task (QLast)=',QLast
c     Write (*,'(A,4I9)') 'iTCnSt_c,nTasks,iTskCan=',
c    &                     iTCnSt_c,nTasks,iTskCan
      TskQ(1,iTskCan) = QLast(1)
      TskQ(2,iTskCan) = QLast(2)
*
      QLast(1)=Not_Used
      QLast(2)=Not_Used
c     Call RecPrt('TskQ',' ',TskQ,2,nTasks)
c     Write (*,*)
c     Call XFlush(6)
*
      Return
      End
