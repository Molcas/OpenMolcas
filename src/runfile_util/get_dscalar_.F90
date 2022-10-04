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

subroutine Get_dScalar_(Label,data)

implicit none
#include "pg_ds_info.fh"
!----------------------------------------------------------------------*
! Arguments                                                            *
!----------------------------------------------------------------------*
character*(*) Label
real*8 data
!----------------------------------------------------------------------*
! Define local variables                                               *
!----------------------------------------------------------------------*
real*8 RecVal(nTocDS)
character*16 RecLab(nTocDS)
integer RecIdx(nTocDS)
!
character*16 CmpLab1
character*16 CmpLab2
integer item
integer i

!----------------------------------------------------------------------*
! Initialize local variables                                           *
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
! Read info from runfile.                                              *
!----------------------------------------------------------------------*
call cRdRun('dScalar labels',RecLab,16*nTocDS)
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
  if (CmpLab1 == CmpLab2) item = i
end do

if (item /= -1) then
  if (Recidx(item) == sSpecialField) then
    write(6,*) '***'
    write(6,*) '*** Warning, reading temporary dScalar field'
    write(6,*) '***   Field: ',Label
    write(6,*) '***'
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
data = RecVal(item)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine Get_dScalar_
