!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Put_lArray(Label,Log,nLog)

character*(*) Label
#include "WrkSpc.fh"
logical log(nLog)

call Allocate_iWork(ip_Tmp,nLog)
call Log2Int(Log,iWork(ip_Tmp),nLog)
call Put_iArray(Label,iWork(ip_Tmp),nLog)
call Free_iWork(ip_Tmp)

return

end subroutine Put_lArray
