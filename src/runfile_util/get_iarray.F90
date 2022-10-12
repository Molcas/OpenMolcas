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
!***********************************************************************
!                                                                      *
! This routine gets array integer data from the runfile.               *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************
!  Get_iArray
!
!> @brief
!>   Read array data from runfile
!> @author Per-Olof Widmark
!>
!> @details
!> This routine gets array double data from the runfile.
!>
!> @param[in]  Label Name of field
!> @param[out] Data  Data to read from runfile
!> @param[in]  nData Length of array
!>
!> @see ::Put_iArray
!***********************************************************************

subroutine Get_iArray(Label,iData,nData)

use RunFile_data, only: i_run_IA_used, lw, nTocIA, sSpecialField
use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: Label
integer(kind=iwp), intent(in) :: nData
integer(kind=iwp), intent(out) :: iData(nData)
integer(kind=iwp) :: i, item, RecIdx(nTocIA), RecLen(nTocIA)
character(len=lw) :: CmpLab1, CmpLab2, RecLab(nTocIA)

!----------------------------------------------------------------------*
! Read info from runfile.                                              *
!----------------------------------------------------------------------*
call cRdRun('iArray labels',RecLab,lw*nTocIA)
call iRdRun('iArray indices',RecIdx,nTocIA)
call iRdRun('iArray lengths',RecLen,nTocIA)
!----------------------------------------------------------------------*
! Locate item                                                          *
!----------------------------------------------------------------------*
item = -1
CmpLab1 = Label
call UpCase(CmpLab1)
do i=1,nTocIA
  CmpLab2 = RecLab(i)
  call UpCase(CmpLab2)
  if (CmpLab1 == CmpLab2) item = i
end do

if (item /= -1) then
  if (Recidx(item) == sSpecialField) then
    write(u6,*) '***'
    write(u6,*) '*** Warning, reading temporary iArray field'
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
i_run_IA_used(item) = i_run_IA_used(item)+1
if (item == -1) call SysAbendMsg('get_iArray','Could not locate: ',Label)
if (RecIdx(item) == 0) call SysAbendMsg('get_iArray','Data not defined: ',Label)
if (Reclen(item) /= nData) call SysAbendMsg('get_iArray','Data of wrong length: ',Label)
call iRdRun(RecLab(item),iData,nData)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine Get_iArray
