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
! Copyright (C) 2008, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine get scalar integer data from the peek/poke buffer.      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: May 2008                                                    *
!                                                                      *
!***********************************************************************
!  Peek_iScalar
!
!> @brief
!>   Get scalar integer data from the peek/poke buffer
!> @author Per-Olof Widmark
!>
!> @details
!> This routine is used to get scalar data of type ``Integer``
!> from the peek/poke buffer. The data items are identified
!> by a text label.
!>
!> @param[in]  Label Name of field
!> @param[out] Data  Data to get from runfile
!***********************************************************************

subroutine Peek_iScalar(Label,data)

implicit none
#include "pp_is_info.fh"
!----------------------------------------------------------------------*
! Arguments                                                            *
!----------------------------------------------------------------------*
character*(*) Label
integer data
!----------------------------------------------------------------------*
! Define local variables                                               *
!----------------------------------------------------------------------*
logical Found
integer indx
integer i

!----------------------------------------------------------------------*
! Initialize local variables                                           *
!----------------------------------------------------------------------*
Found = .false.
!write(6,'(2a)') 'peek_iscalar: Label is ',Label
!write(6,'(a,i8)') 'peek_iscalar: is_no is ',is_no
!----------------------------------------------------------------------*
! Locate item                                                          *
!----------------------------------------------------------------------*
indx = -1
do i=1,is_no
  if (is_label(i) == Label) indx = i
end do
!write(6,'(a,i8)') 'peek_iscalar: indx is ',indx
!----------------------------------------------------------------------*
! Get data from buffer.                                                *
!----------------------------------------------------------------------*
if (indx == -1) then
  if (is_no >= nTabIS) then
    call SysAbendMsg('Peek_iScalar','Too many fields','Increase nTabIS and recompile')
  end if
  is_no = is_no+1
  indx = is_no
  call Qpg_iScalar(Label,Found)
  if (Found) then
    call Get_iScalar(Label,data)
  else
    call SysAbendMsg('Peek_iScalar','Field not found',Label)
  end if
  is_label(indx) = Label
  is_value(indx) = data
else
  data = is_value(indx)
end if
!write(6,'(a,i8)') 'peek_iscalar: Data is ',Data
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*

return

end subroutine Peek_iScalar
