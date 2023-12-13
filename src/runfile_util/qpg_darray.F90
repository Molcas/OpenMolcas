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
! This routine queries the existence of array data on runfile.         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************
!  Qpg_dArray
!
!> @brief
!>   Check if a field exist on the runfile
!> @author Per-Olof Widmark
!>
!> @details
!> This routine queries the existence of array data on runfile.
!>
!> @param[in]  Label Name of field
!> @param[out] Found Was the field found
!> @param[out] nData Length of field
!>
!> @see ::Put_dArray
!***********************************************************************

subroutine Qpg_dArray(Label,Found,nData)

use RunFile_data, only: lw, nTocDA, sSpecialField
use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: Label
logical(kind=iwp), intent(out) :: Found
integer(kind=iwp), intent(out) :: nData
integer(kind=iwp) :: i, item, iTmp, nTmp, RecIdx(nTocDA), RecLen(nTocDA)
character(len=lw) :: CmpLab1, CmpLab2, RecLab(nTocDA)

!----------------------------------------------------------------------*
! Read info from runfile.                                              *
!----------------------------------------------------------------------*
call ffRun('dArray labels',nTmp,iTmp)
if (nTmp == 0) then
  Found = .false.
  nData = 0
  return
end if
call cRdRun('dArray labels',RecLab,lw*nTocDA)
call iRdRun('dArray indices',RecIdx,nTocDA)
call iRdRun('dArray lengths',RecLen,nTocDA)
!----------------------------------------------------------------------*
! Locate item                                                          *
!----------------------------------------------------------------------*
item = -1
CmpLab1 = Label
call UpCase(CmpLab1)
do i=1,nTocDA
  CmpLab2 = RecLab(i)
  call UpCase(CmpLab2)
  if (CmpLab1 == CmpLab2) item = i
end do

! Do we read an old temporary field?

if (item /= -1) then
  if (RecIdx(item) == sSpecialField) then
    write(u6,*) '***'
    write(u6,*) '*** Warning, querying temporary dArray field'
    write(u6,*) '***   Field: ',Label
    write(u6,*) '***'
#   ifndef _DEVEL_
    call AbEnd()
#   endif
  end if
end if
!----------------------------------------------------------------------*
! Did we manage to find it?                                            *
!----------------------------------------------------------------------*
Found = .true.
if (item == -1) Found = .false.
if ((item /= -1) .and. (RecIdx(item) == 0)) Found = .false.
if (Found) then
  nData = RecLen(item)
else
  nData = 0
end if
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine Qpg_dArray
