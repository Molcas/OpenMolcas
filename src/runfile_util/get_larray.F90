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

subroutine Get_lArray(Label,Logc,nLog)

use Definitions, only: iwp

implicit none
character(len=*) :: Label
integer(kind=iwp) :: nLog
logical(kind=iwp) :: Logc(nLog)
#include "WrkSpc.fh"
integer(kind=iwp) :: ip_Tmp

call Allocate_iWork(ip_Tmp,nLog)
call Get_iArray(Label,iWork(ip_Tmp),nLog)
call Int2Log(iWork(ip_Tmp),Logc,nLog)
call Free_iWork(ip_Tmp)

return

end subroutine Get_lArray
