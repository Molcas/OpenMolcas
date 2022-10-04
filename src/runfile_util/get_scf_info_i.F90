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

subroutine Get_SCF_Info_I(iCase,ipI,nI)

implicit real*8(A-H,O-Z)
#include "WrkSpc.fh"
character*24 Label
logical Found

Label = 'SCFInfoI'
if (iCase == 1) Label = 'SCFInfoI_ab'
call qpg_iArray(Label,Found,nI)
if ((.not. Found) .or. (nI == 0)) call SysAbendMsg('get_scf_info_i','Did not find:',Label)
call GetMem('ipI','Allo','Inte',ipI,nI)
call Get_iArray(Label,iWork(ipI),nI)

return

end subroutine Get_SCF_Info_I
