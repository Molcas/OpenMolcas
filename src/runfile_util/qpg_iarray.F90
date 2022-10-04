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
!  Qpg_iArray
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
!> @see ::Put_iArray
!***********************************************************************

subroutine Qpg_iArray(Label,Found,nData)

implicit none
#include "pg_ia_info.fh"
!----------------------------------------------------------------------*
! Arguments                                                            *
!----------------------------------------------------------------------*
character*(*) Label
logical Found
integer nData
!----------------------------------------------------------------------*
! Define local variables                                               *
!----------------------------------------------------------------------*
character*16 RecLab(nTocIA)
integer RecIdx(nTocIA)
integer RecLen(nTocIA)
character*16 CmpLab1
character*16 CmpLab2
integer item
integer nTmp, iTmp
integer i
!----------------------------------------------------------------------*
! Initialize local variables                                           *
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
! Read info from runfile.                                              *
!----------------------------------------------------------------------*
call ffRun('iArray labels',nTmp,iTmp)
if (nTmp == 0) then
  Found = .false.
  nData = 0
  return
end if
call cRdRun('iArray labels',RecLab,16*nTocIA)
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

! Do we read an old temporary field?

if (item /= -1) then
  if (RecIdx(item) == sSpecialField) then
    write(6,*) '***'
    write(6,*) '*** Warning, querying temporary iArray field'
    write(6,*) '***   Field: ',Label
    write(6,*) '***'
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

end subroutine Qpg_iArray
