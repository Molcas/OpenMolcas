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
! This routine gets scalar double data from the runfile.               *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************
!  Get_dScalar
!
!> @brief
!>   Get scalar data form runfile
!> @author Per-Olof Widmark
!>
!> @details
!> This routine gets scalar double data from the runfile.
!>
!> @param[in]  Label Name of field
!> @param[out] Data  Data to get from runfile
!>
!> @see ::Put_dScalar
!***********************************************************************

subroutine Get_dScalar(Label,data)

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
character*16 CmpLab
integer dfirst, i
data dfirst/0/
save dfirst

if (dfirst == 0) then
  dfirst = 1
  num_DS_init = 0
  do i=1,nTocDS
    iLbl_DS_inmem(i) = ' '
    DS_init(i) = 0
  end do
end if

CmpLab = Label
call UpCase(CmpLab)

do i=1,num_DS_init
  if (iLbl_DS_inmem(i) == CmpLab) then
    if (DS_init(i) /= 0) then
      data = i_DS_inmem(i)
      return
    end if
  end if
end do

call Get_dScalar_(Label,data)
num_DS_init = num_DS_init+1

if (num_DS_init > nTocDS) then
# ifdef _DEBUGPRINT_
  do i=1,num_DS_init
    write(6,*) iLbl_DS_inmem(i),DS_init(i),i_DS_inmem(i),CmpLab
  end do
# endif
  call Abend()
end if

iLbl_DS_inmem(num_DS_init) = CmpLab
DS_init(num_DS_init) = 1
i_DS_inmem(num_DS_init) = data

return

end subroutine Get_dScalar
