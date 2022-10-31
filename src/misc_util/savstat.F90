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

subroutine SavStat(iFld,val,op)

use Para_Info, only: MyRank
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: iFld
real(kind=wp) :: val
character(len=*) :: op
#include "timtra.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iad

if (nfld_stat == 0) return
if (iFld > nfld_Stat) then
  call WarningMessage(2,'SavStat: iFld > nfld_stat')
  write(u6,*) 'iFld=',iFld
  write(u6,*) 'nFld_Stat=',nFld_Stat
  call Abend()
end if
iad = iGAStat+myrank*nFld_stat+iFld-1
if (op == '+') then
  Work(iad) = Work(iad)+val
else if (op == '-') then
  Work(iad) = Work(iad)-val
else if (op == '=') then
  Work(iad) = val
end if

return

end subroutine SavStat
