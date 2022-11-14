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
! This routine queries the existence of scalar data on runfile.        *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************
!  Qpg_iScalar
!
!> @brief
!>   Check if a field is available on runfile
!> @author Per-Olof Widmark
!>
!> @details
!> This routine queries the existence of scalar data on runfile.
!>
!> @param[in]  Label Name of field
!> @param[out] Found Was the field found
!>
!> @see ::Put_iScalar
!***********************************************************************

subroutine Qpg_iScalar(Label,Found)

use RunFile_data, only: lw, nTocIS, sSpecialField
use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: Label
logical(kind=iwp), intent(out) :: Found
integer(kind=iwp) :: i, item, iTmp, nTmp, RecIdx(nTocIS), RecVal(nTocIS)
character(len=lw) :: CmpLab1, CmpLab2, RecLab(nTocIS)

!----------------------------------------------------------------------*
! Read info from runfile.                                              *
!----------------------------------------------------------------------*
call ffRun('iScalar labels',nTmp,iTmp)
if (nTmp == 0) then
  Found = .false.
  return
end if
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
  if (CmpLab1 == CmpLab2) item = i
end do

! Do we read an old temporary field?

if (item /= -1) then
  if (RecIdx(item) == sSpecialField) then
    write(u6,*) '***'
    write(u6,*) '*** Warning, querying temporary iScalar field'
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
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine Qpg_iScalar
