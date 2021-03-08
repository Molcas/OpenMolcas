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
!> @param[out] Val   Data to get from runfile
!***********************************************************************

subroutine Peek_iScalar(Label,val)

use peekpoke, only: is_label, is_value, is_no
use Definitions, only: iwp

implicit none
!----------------------------------------------------------------------*
! Arguments                                                            *
!----------------------------------------------------------------------*
character(len=*), intent(in) :: Label
integer(kind=iwp), intent(out) :: val
!----------------------------------------------------------------------*
! Define local variables                                               *
!----------------------------------------------------------------------*
logical(kind=iwp) :: Found
integer(kind=iwp) :: indx, i

!----------------------------------------------------------------------*
! Initialize local variables                                           *
!----------------------------------------------------------------------*
Found = .false.
!write(u6,'(2a)') 'peek_iscalar: Label is ',Label
!write(u6,'(a,i8)') 'peek_iscalar: is_no is ',is_no
!----------------------------------------------------------------------*
! Locate item                                                          *
!----------------------------------------------------------------------*
indx = -1
do i=1,is_no
  if (is_label(i) == Label) then
    indx = i
    exit
  end if
end do
!write(u6,'(a,i8)') 'peek_iscalar: indx is ',indx
!----------------------------------------------------------------------*
! Get data from buffer.                                                *
!----------------------------------------------------------------------*
if (indx == -1) then
  if (is_no >= size(is_value)) then
    call SysAbendMsg('Peek_iScalar','Too many fields','Increase nTabIS and recompile')
  end if
  is_no = is_no+1
  indx = is_no
  call Qpg_iScalar(Label,Found)
  if (Found) then
    call Get_iScalar(Label,val)
  else
    call SysAbendMsg('Peek_iScalar','Field not found',Label)
  end if
  is_label(indx) = Label
  is_value(indx) = val
else
  val = is_value(indx)
end if
!write(u6,'(a,i8)') 'peek_iscalar: Val is ',Val
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
return

end subroutine Peek_iScalar
