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
! This routine puts array integer data to the runfile.                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************
!  Put_iArray
!
!> @brief
!>   Add/update array data in runfile
!> @author Per-Olof Widmark
!>
!> @details
!> This routine is used to put array data of type
!> ``Integer`` into the runfile. The data items are
!> identified by the label. Below is a list of the
!> data items that are recognized. The labels are
!> case insensitive and significant to 16 characters.
!>
!> For development purposes you can use an unsupported
!> label. Whenever such a field is accessed a warning
!> message is printed in the output, to remind the
!> developer to update this routine.
!>
!> @param[in] Label Name of field
!> @param[in] iData Data to put on runfile
!> @param[in] nData Length of array
!***********************************************************************

subroutine Put_iArray(Label,iData,nData)

use RunFile_data, only: LabelsIA, lw, nTocIA, sNotUsed, sRegularField, sSpecialField
use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: Label
integer(kind=iwp), intent(in) :: nData, iData(nData)
integer(kind=iwp) :: i, item, iTmp, nTmp, RecIdx(nTocIA) = 0, RecLen(nTocIA) = 0
character(len=lw) :: CmpLab1, CmpLab2, RecLab(nTocIA) = ''

!----------------------------------------------------------------------*
! Do setup if this is the first call.                                  *
!----------------------------------------------------------------------*
call ffRun('iArray labels',nTmp,iTmp)
if (nTmp == 0) then
  RecLab(:) = LabelsIA
  RecIdx(:) = sNotUsed
  RecLen(:) = 0
  call cWrRun('iArray labels',RecLab,lw*nTocIA)
  call iWrRun('iArray indices',RecIdx,nTocIA)
  call iWrRun('iArray lengths',RecLen,nTocIA)
else
  call cRdRun('iArray labels',RecLab,lw*nTocIA)
  call iRdRun('iArray indices',RecIdx,nTocIA)
  call iRdRun('iArray lengths',RecLen,nTocIA)
end if
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

! Do we create a new temporary field?

if (item == -1) then
  do i=1,nTocIA
    if (RecLab(i) == ' ') item = i
  end do
  if (item /= -1) then
    RecLab(item) = Label
    RecIdx(item) = sSpecialField
    call cWrRun('iArray labels',RecLab,lw*nTocIA)
    call iWrRun('iArray indices',RecIdx,nTocIA)
  end if
end if

! Is this a temporary field?

if (item /= -1) then
  if (Recidx(item) == sSpecialField) then
    write(u6,*) '***'
    write(u6,*) '*** Warning, writing temporary iArray field'
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
if (item == -1) call SysAbendMsg('put_iArray','Could not locate',Label)
call iWrRun(RecLab(item),iData,nData)
if (RecIdx(item) == 0) then
  RecIdx(item) = sRegularField
  call iWrRun('iArray indices',RecIdx,nTocIA)
end if
if (RecLen(item) /= nData) then
  RecLen(item) = nData
  call iWrRun('iArray lengths',RecLen,nTocIA)
end if
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine Put_iArray
