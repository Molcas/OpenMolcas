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
! This routine gets scalar integer data from the runfile.              *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************
!  Get_iScalar
!
!> @brief
!>   Get scalar data from runfile
!> @author Per-Olof Widmark
!>
!> @details
!> This routine gets scalar integer data from the runfile.
!>
!> @param[in]  Label Name of field
!> @param[out] Data  Data to get from runfile
!>
!> @see ::Put_iScalar
!***********************************************************************

subroutine Get_iScalar(Label,data)

#include "pg_is_info.fh"
!----------------------------------------------------------------------*
! Arguments                                                            *
!----------------------------------------------------------------------*
character*(*) Label
integer data
!----------------------------------------------------------------------*
! Define local variables                                               *
!----------------------------------------------------------------------*
integer ifirst, i
character*16 CmpLab
data ifirst/0/
save ifirst

if (ifirst == 0) then
  ifirst = 1
  num_IS_init = 0
  do i=1,nTocIS
    iLbl_IS_inmem(i) = ' '
    IS_init(i) = 0
  end do
end if

CmpLab = Label
call UpCase(CmpLab)

do i=1,num_IS_init
  if (iLbl_IS_inmem(i) == CmpLab) then
    if (IS_init(i) /= 0) then
      data = i_IS_inmem(i)
      return
    end if
  end if
end do

call Get_iScalar_(Label,data)
num_IS_init = num_IS_init+1

if (num_IS_init > nTocIS) then
# ifdef _DEBUGPRINT__
  do i=1,num_IS_init
    write(6,*) iLbl_IS_inmem(i),IS_init(i),i_IS_inmem(i)
  end do
# endif
  call Abend()
end if

iLbl_IS_inmem(num_IS_init) = CmpLab
IS_init(num_IS_init) = 1
i_IS_inmem(num_IS_init) = data

return

end subroutine Get_iScalar
