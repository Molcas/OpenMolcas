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
! This routine put scalar integer data to the peek/poke buffer.        *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: May 2008                                                    *
!                                                                      *
!***********************************************************************
!  Poke_iScalar
!
!> @brief
!>   Put scalar integer data to the peek/poke buffer
!> @author Per-Olof Widmark
!>
!> @details
!> This routine is used to put scalar data of type ``Integer``
!> in the peek/poke buffer. The data items are identified
!> by a text label.
!>
!> @param[in] Label Name of field
!> @param[in] Val   Data to put on runfile
!***********************************************************************

subroutine Poke_iScalar(Label,val)

use peekpoke, only: is_label, is_value, is_no
use Definitions, only: iwp

implicit none
!----------------------------------------------------------------------*
! Arguments                                                            *
!----------------------------------------------------------------------*
character(len=*), intent(in) :: Label
integer(kind=iwp), intent(in) :: val
!----------------------------------------------------------------------*
! Define local variables                                               *
!----------------------------------------------------------------------*
integer(kind=iwp) :: indx, i

!----------------------------------------------------------------------*
! Initialize local variables                                           *
!----------------------------------------------------------------------*
!write(u6,'(2a)') 'poke_iscalar: Label is ',Label
!write(u6,'(a,i8)') 'poke_iscalar: Val is ',Val
!----------------------------------------------------------------------*
! Locate item                                                          *
!----------------------------------------------------------------------*
indx = -1
do i=1,is_no
  if (is_label(i) == Label) indx = i
end do
!write(u6,'(a,i8)') 'poke_iscalar: indx is ',indx
!----------------------------------------------------------------------*
! Save data in buffer.                                                 *
!----------------------------------------------------------------------*
if (indx == -1) then
  if (is_no >= size(is_value)) then
    call SysAbendMsg('Poke_iScalar','Too many fields','Increase nTabIS and recompile')
  end if
  is_no = is_no+1
  indx = is_no
end if
is_label(indx) = Label
is_value(indx) = val
!write(u6,'(a,i8)') 'poke_iscalar: is_no is ',is_no
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
return

end subroutine Poke_iScalar
