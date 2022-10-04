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
! This routine gets array character data from the runfile.             *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************
!  Get_cArray
!
!> @brief
!>   Read array data from runfile
!> @author Per-Olof Widmark
!>
!> @details
!> This routine gets array character data from the runfile.
!>
!> @param[in]  Label Name of field
!> @param[out] Data  Data to read from runfile
!> @param[in]  nData Length of array
!>
!> @see ::Put_cArray
!***********************************************************************

subroutine Get_cArray(Label,data,nData)

implicit none
#include "pg_ca_info.fh"
!----------------------------------------------------------------------*
! Arguments                                                            *
!----------------------------------------------------------------------*
character*(*) Label
integer nData
character*(*) data(nData)
!----------------------------------------------------------------------*
! Define local variables                                               *
!----------------------------------------------------------------------*
character*16 RecLab(nTocCA)
integer RecIdx(nTocCA)
integer RecLen(nTocCA)
character*16 CmpLab1
character*16 CmpLab2
integer item
integer i

!----------------------------------------------------------------------*
! Read info from runfile.                                              *
!----------------------------------------------------------------------*
call cRdRun('cArray labels',RecLab,16*nTocCA)
call iRdRun('cArray indices',RecIdx,nTocCA)
call iRdRun('cArray lengths',RecLen,nTocCA)
!----------------------------------------------------------------------*
! Locate item                                                          *
!----------------------------------------------------------------------*
item = -1
CmpLab1 = Label
call UpCase(CmpLab1)
do i=1,nTocCA
  CmpLab2 = RecLab(i)
  call UpCase(CmpLab2)
  if (CmpLab1 == CmpLab2) item = i
end do

if (item /= -1) then
  if (Recidx(item) == sSpecialField) then
    write(6,*) '***'
    write(6,*) '*** Warning, reading temporary cArray field'
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
i_run_CA_used(item) = i_run_CA_used(item)+1
if (item == -1) call SysAbendMsg('get_cArray','Could not locate: ',Label)
if (RecIdx(item) == 0) call SysAbendMsg('get_cArray','Data not defined: ',Label)
if (Reclen(item) /= nData) call SysAbendMsg('get_cArray','Data of wrong length: ',Label)
call cRdRun(RecLab(item),data,nData)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine Get_cArray
