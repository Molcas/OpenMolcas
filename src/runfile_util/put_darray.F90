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
! This routine puts array double data to the runfile.                  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************
!  Put_dArray
!
!> @brief
!>   Add/update array data in runfile
!> @author Per-Olof Widmark
!>
!> @details
!> This routine is used to put array data of type
!> ``Real*8`` into the runfile. The data items are
!> identified by the \p label. Below is a list of the
!> data items that are recognized. The labels are
!> case insensitive and significant to 16 characters.
!> (May change to 24 characters.)
!>
!> @warning
!> Naming convention is under development
!> and labels may change to the next version.
!>
!> For development purposes you can use an unsupported
!> label. Whenever such a field is accessed a warning
!> message is printed in the output, to remind the
!> developer to update this routine.
!>
!> @param[in] Label Name of field
!> @param[in] rData Data to put on runfile
!> @param[in] nData Length of array
!***********************************************************************

subroutine Put_dArray(Label,rData,nData)

use RunFile_data, only: LabelsDA, lw, nTocDA, sNotUsed, sRegularField, sSpecialField
use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: Label
integer(kind=iwp), intent(in) :: nData
real(kind=wp), intent(in) :: rData(nData)
integer(kind=iwp) :: i, item, iTmp, nTmp, RecIdx(nTocDA) = 0, RecLen(nTocDA) = 0
character(len=lw) :: CmpLab1, CmpLab2, RecLab(nTocDA) = ''

!----------------------------------------------------------------------*
! Do setup if this is the first call.                                  *
!----------------------------------------------------------------------*
call ffRun('dArray labels',nTmp,iTmp)
if (nTmp == 0) then
  RecLab(:) = LabelsDA
  RecIdx(:) = sNotUsed
  RecLen(:) = 0
  call cWrRun('dArray labels',RecLab,lw*nTocDA)
  call iWrRun('dArray indices',RecIdx,nTocDA)
  call iWrRun('dArray lengths',RecLen,nTocDA)
else
  call cRdRun('dArray labels',RecLab,lw*nTocDA)
  call iRdRun('dArray indices',RecIdx,nTocDA)
  call iRdRun('dArray lengths',RecLen,nTocDA)
end if
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

! Do we create a new temporary field?

if (item == -1) then
  do i=1,nTocDA
    if (RecLab(i) == ' ') item = i
  end do
  if (item /= -1) then
    RecLab(item) = Label
    RecIdx(item) = sSpecialField
    call cWrRun('dArray labels',RecLab,lw*nTocDA)
    call iWrRun('dArray indices',RecIdx,nTocDA)
  end if
end if

! Is this a temporary field?

if (item /= -1) then
  if (Recidx(item) == sSpecialField) then
    write(u6,*) '***'
    write(u6,*) '*** Warning, writing temporary dArray field'
    write(u6,*) '***   Field: ',Label
    write(u6,*) '***'
#   ifndef _DEVEL_
    call AbEnd()
#   endif
  end if
end if
!----------------------------------------------------------------------*
! Write data to disk.                                                  *
!----------------------------------------------------------------------*
if (item == -1) call SysAbendMsg('put_dArray','Could not locate',Label)
call dWrRun(RecLab(item),rData,nData)
!write(u6,*) rData
if (RecIdx(item) == 0) then
  RecIdx(item) = sRegularField
  call iWrRun('dArray indices',RecIdx,nTocDA)
end if
if (RecLen(item) /= nData) then
  RecLen(item) = nData
  call iWrRun('dArray lengths',RecLen,nTocDA)
end if
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine Put_dArray
