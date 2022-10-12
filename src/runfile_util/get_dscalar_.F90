!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2003, Per-Olof Widmark                                 *
!***********************************************************************

subroutine Get_dScalar_(Label,rData)

use RunFile_data, only: i_run_DS_used, lw, nTocDS, sSpecialField
use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: Label
real(kind=wp), intent(out) :: rData
integer(kind=iwp) :: i, item, RecIdx(nTocDS)
real(kind=wp) :: RecVal(nTocDS)
character(len=lw) :: CmpLab1, CmpLab2, RecLab(nTocDS)

!----------------------------------------------------------------------*
! Read info from runfile.                                              *
!----------------------------------------------------------------------*
call cRdRun('dScalar labels',RecLab,lw*nTocDS)
call dRdRun('dScalar values',RecVal,nTocDS)
call iRdRun('dScalar indices',RecIdx,nTocDS)
!----------------------------------------------------------------------*
! Locate item                                                          *
!----------------------------------------------------------------------*
item = -1
CmpLab1 = Label
call UpCase(CmpLab1)
do i=1,nTocDS
  CmpLab2 = RecLab(i)
  call UpCase(CmpLab2)
  if (CmpLab1 == CmpLab2) then
    item = i
    exit
  end if
end do

if (item /= -1) then
  if (Recidx(item) == sSpecialField) then
    write(u6,*) '***'
    write(u6,*) '*** Warning, reading temporary dScalar field'
    write(u6,*) '***   Field: ',Label
    write(u6,*) '***'
#   ifndef _DEVEL_
    call AbEnd()
#   endif
  end if
end if
!----------------------------------------------------------------------*
! Transfer data to arguments                                           *
!----------------------------------------------------------------------*
i_run_DS_used(item) = i_run_DS_used(item)+1
if (item == -1) call SysAbendMsg('get_dScalar','Could not locate: ',Label)
if (RecIdx(item) == 0) call SysAbendMsg('get_dScalar','Data not defined: ',Label)
rData = RecVal(item)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine Get_dScalar_
