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

subroutine Get_iScalar_(Label,iData)

use RunFile_data, only: i_run_IS_used, lw, nTocIS, sSpecialField
use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: Label
integer(kind=iwp), intent(out) :: iData
integer(kind=iwp) :: i, item, RecIdx(nTocIS), RecVal(nTocIS)
character(len=lw) :: CmpLab1, CmpLab2, RecLab(nTocIS)

!----------------------------------------------------------------------*
! Read info from runfile.                                              *
!----------------------------------------------------------------------*
call cRdRun('iScalar labels',RecLab,lw*nTocIS)
call iRdRun('iScalar values',RecVal,nTocIS)
call iRdRun('iScalar indices',RecIdx,nTocIS)
!----------------------------------------------------------------------*
! Locate item                                                          *
!----------------------------------------------------------------------*
item = -1
CmpLab1 = Label
call UpCase(CmpLab1)
do i=1,nTocIS
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
    write(u6,*) '*** Warning, reading temporary iScalar field'
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
i_run_IS_used(item) = i_run_IS_used(item)+1
if (item == -1) call SysAbendMsg('get_iScalar','Could not locate: ',Label)
if (RecIdx(item) == 0) call SysAbendMsg('get_iScalar','Data not defined: ',Label)
iData = RecVal(item)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine Get_iScalar_
