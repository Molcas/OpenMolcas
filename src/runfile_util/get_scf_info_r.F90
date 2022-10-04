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

subroutine Get_SCF_Info_R(iCase,ipR,nR)

implicit real*8(A-H,O-Z)
#include "WrkSpc.fh"
character*24 Label
logical Found

Label = 'SCFInfoR'
if (iCase == 1) Label = 'SCFInfoR_ab'
call qpg_dArray(label,Found,nR)
if ((.not. Found) .or. (nR == 0)) call SysAbendmsg('get_scf_info_r','Did not find:',Label)
call GetMem('ipR','Allo','Real',ipR,nR)
call get_dArray(label,Work(ipR),nR)

return

end subroutine Get_SCF_Info_R
