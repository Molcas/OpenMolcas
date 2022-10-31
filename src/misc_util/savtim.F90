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

subroutine SavTim(iFld,TCPU,TWall)

use Para_Info, only: MyRank
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: iFld
real(kind=wp) :: TCPU, TWall
#include "timtra.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iad

if (nfld_tim == 0) return
if (iFld > nfld_tim) then
  call WarningMessage(2,'SavTim: iFld > nfld_tim')
  write(u6,*) 'iFld=',iFld
  write(u6,*) 'nFld_tim=',nFld_tim
  call Abend()
end if
iad = iGATim+myrank*nFld_Tim*2+iFld-1
Work(iad) = Work(iad)+TCPU
Work(iad+nFld_Tim) = Work(iad+nFld_Tim)+TWall

return

end subroutine SavTim
