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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine puts scalar integer data to the runfile.                *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************
!  Put_iScalar
!
!> @brief
!>   Add/update scalar data in runfile
!> @author Per-Olof Widmark
!>
!> @details
!> This routine is used to put scalar data of type
!> ``Integer`` into the runfile. The data items are
!> identified by the \p label. Below is a list of the
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
!***********************************************************************

subroutine Put_iScalar(Label,iData)

use RunFile_data, only: IS_cache, LabelsIS, lw, nTocIS, num_IS_init, sNotUsed, sRegularField, sSpecialField
use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: Label
integer(kind=iwp), intent(in) :: iData
integer(kind=iwp) :: i, item, iTmp, nData, RecIdx(nTocIS) = 0, RecVal(nTocIS) = 0
character(len=lw) :: CmpLab1, CmpLab2, RecLab(nTocIS) = ''

!----------------------------------------------------------------------*
! Do setup if this is the first call.                                  *
!----------------------------------------------------------------------*
call ffRun('iScalar labels',nData,iTmp)
if (nData == 0) then
  RecLab(:) = LabelsIS
  RecVal(:) = 0
  RecIdx(:) = sNotUsed
  call cWrRun('iScalar labels',RecLab,lw*nTocIS)
  call iWrRun('iScalar values',RecVal,nTocIS)
  call iWrRun('iScalar indices',RecIdx,nTocIS)
else
  call cRdRun('iScalar labels',RecLab,lw*nTocIS)
  call iRdRun('iScalar values',RecVal,nTocIS)
  call iRdRun('iScalar indices',RecIdx,nTocIS)
end if
!----------------------------------------------------------------------*
! Locate item                                                          *
!----------------------------------------------------------------------*
item = -1
CmpLab1 = Label
call UpCase(CmpLab1)
do i=1,nTocIS
  CmpLab2 = RecLab(i)
  call UpCase(CmpLab2)
  if (CmpLab1 == CmpLab2) item = i
end do

! Do we create a new temporary field?

if (item == -1) then
  do i=1,nTocIS
    if (RecLab(i) == ' ') item = i
  end do
  if (item /= -1) then
    RecLab(item) = Label
    RecIdx(item) = sSpecialField
    call cWrRun('iScalar labels',RecLab,lw*nTocIS)
    call iWrRun('iScalar indices',RecIdx,nTocIS)
  end if
end if

! Is this a temporary field?

if (item /= -1) then
  if (Recidx(item) == sSpecialField) then
    write(u6,*) '***'
    write(u6,*) '*** Warning, writing temporary iScalar field'
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
if (item == -1) call SysAbendMsg('put_iScalar','Could not locate',Label)
RecVal(item) = iData
call iWrRun('iScalar values',RecVal,nTocIS)
if (RecIdx(item) == 0) then
  RecIdx(item) = sRegularField
  call iWrRun('iScalar indices',RecIdx,nTocIS)
end if
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
do i=1,num_IS_init
  if (IS_cache(i)%lab == CmpLab1) then
    IS_cache(i)%val = iData
    exit
  end if
end do

return

end subroutine Put_iScalar
