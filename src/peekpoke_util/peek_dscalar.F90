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
! This routine get scalar double data from the peek/poke buffer.       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: May 2008                                                    *
!                                                                      *
!***********************************************************************
!  Peek_dScalar
!
!> @brief
!>   Get scalar double data from the peek/poke buffer
!> @author Per-Olof Widmark
!>
!> @details
!> This routine is used to get scalar data of type ``Real*8``
!> from the peek/poke buffer. The data items are identified
!> by a text label.
!>
!> @param[in]  Label Name of field
!> @param[out] Val   Data to get from runfile
!***********************************************************************

subroutine Peek_dScalar(Label,val)

use peekpoke, only: ds_label, ds_value, ds_no
use Definitions, only: wp, iwp

implicit none
!----------------------------------------------------------------------*
! Arguments                                                            *
!----------------------------------------------------------------------*
character(len=*), intent(in) :: Label
real(kind=wp), intent(out) :: val
!----------------------------------------------------------------------*
! Define local variables                                               *
!----------------------------------------------------------------------*
logical(kind=iwp) :: Found
integer(kind=iwp) :: indx, i

!----------------------------------------------------------------------*
! Initialize local variables                                           *
!----------------------------------------------------------------------*
Found = .false.
!write(u6,'(2a)') 'peek_dscalar: Label is ',Label
!write(u6,'(a,i8)') 'peek_dscalar: ds_no is ',ds_no
!----------------------------------------------------------------------*
! Locate item                                                          *
!----------------------------------------------------------------------*
indx = -1
do i=1,ds_no
  if (ds_label(i) == Label) then
    indx = i
    exit
  end if
end do
!write(u6,'(a,i8)') 'peek_dscalar: indx is ',indx
!----------------------------------------------------------------------*
! Get data from buffer.                                                *
!----------------------------------------------------------------------*
if (indx == -1) then
  if (ds_no >= size(ds_value)) then
    call SysAbendMsg('Peek_dScalar','Too many fields','Increase nTabDS and recompile')
  end if
  ds_no = ds_no+1
  indx = ds_no
  call Qpg_dScalar(Label,Found)
  if (Found) then
    call Get_dScalar(Label,val)
  else
    call SysAbendMsg('Peek_dScalar','Field not found',Label)
  end if
  ds_label(indx) = Label
  ds_value(indx) = val
else
  val = ds_value(indx)
end if
!write(u6,'(a,e20.8)') 'peek_dscalar: Val is ',Val
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
return

end subroutine Peek_dScalar
