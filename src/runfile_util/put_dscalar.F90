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
! This routine puts scalar double data to the runfile.                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************
!  Put_dScalar
!
!> @brief
!>   To add/update scalar data in runfile
!> @author Per-Olof Widmark
!>
!> @details
!> This routine is used to put scalar data of type
!> ``Real*8`` into the runfile. The data items are
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
!> @param[in] rData Data to put on runfile
!***********************************************************************

subroutine Put_dScalar(Label,rData)

use RunFile_data, only: DS_cache, LabelsDS, lw, nTocDS, num_DS_init, sNotUsed, sRegularField, sSpecialField
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: Label
real(kind=wp), intent(in) :: rData
integer(kind=iwp) :: i, item, iTmp, nData, RecIdx(nTocDS)
real(kind=wp) :: RecVal(nTocDS)
character(len=lw) :: CmpLab1, CmpLab2, RecLab(nTocDS)

!----------------------------------------------------------------------*
! Do setup if this is the first call.                                  *
!----------------------------------------------------------------------*
! start pow mod ---
!write(6,'(3a)') 'Runfile: put_dscalar field "',Label,'"'
! end pow mod ---
call ffRun('dScalar labels',nData,iTmp)
if (nData == 0) then
  RecLab(:) = LabelsDS
  RecVal(:) = Zero
  RecIdx(:) = sNotUsed
  call cWrRun('dScalar labels',RecLab,lw*nTocDS)
  call dWrRun('dScalar values',RecVal,nTocDS)
  call iWrRun('dScalar indices',RecIdx,nTocDS)
else
  call cRdRun('dScalar labels',RecLab,lw*nTocDS)
  call dRdRun('dScalar values',RecVal,nTocDS)
  call iRdRun('dScalar indices',RecIdx,nTocDS)
end if
!----------------------------------------------------------------------*
! Locate item                                                          *
!----------------------------------------------------------------------*
item = -1
CmpLab1 = Label
call UpCase(CmpLab1)
do i=1,nTocDS
  CmpLab2 = RecLab(i)
  call UpCase(CmpLab2)
  if (CmpLab1 == CmpLab2) item = i
end do

! Do we create a new temporary field?

if (item == -1) then
  do i=1,nTocDS
    if (RecLab(i) == ' ') item = i
  end do
  if (item /= -1) then
    RecLab(item) = Label
    RecIdx(item) = sSpecialField
    call cWrRun('dScalar labels',RecLab,lw*nTocDS)
    call iWrRun('dScalar indices',RecIdx,nTocDS)
  end if
end if

! Is this a temporary field?

if (item /= -1) then
  if (Recidx(item) == sSpecialField) then
    write(u6,*) '***'
    write(u6,*) '*** Warning, writing temporary dScalar field'
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
if (item == -1) call SysAbendMsg('put_dScalar','Could not locate',Label)
RecVal(item) = rData
call dWrRun('dScalar values',RecVal,nTocDS)
if (RecIdx(item) == 0) then
  RecIdx(item) = sRegularField
  call iWrRun('dScalar indices',RecIdx,nTocDS)
end if
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
do i=1,num_DS_init
  if (DS_cache(i)%lab == CmpLab1) then
    DS_cache(i)%val = rData
    exit
  end if
end do

return

end subroutine Put_dScalar
