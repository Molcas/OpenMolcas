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
! This routine put scalar double data to the peek/poke buffer.         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: May 2008                                                    *
!                                                                      *
!***********************************************************************
!  Poke_dScalar
!
!> @brief
!>   Put scalar double data to the peek/poke buffer
!> @author Per-Olof Widmark
!>
!> @details
!> This routine is used to put scalar data of type ``Real*8``
!> in the peek/poke buffer. The data items are identified
!> by a text label.
!>
!> @param[in] Label Name of field
!> @param[in] Val   Data to put on runfile
!***********************************************************************

subroutine Poke_dScalar(Label,val)

use peekpoke, only: ds_label, ds_value, ds_no
use Definitions, only: wp, iwp

implicit none
!----------------------------------------------------------------------*
! Arguments                                                            *
!----------------------------------------------------------------------*
character(len=*), intent(in) :: Label
real(kind=wp), intent(in) :: val
!----------------------------------------------------------------------*
! Define local variables                                               *
!----------------------------------------------------------------------*
integer(kind=iwp) :: indx, i

!----------------------------------------------------------------------*
! Initialize local variables                                           *
!----------------------------------------------------------------------*
!write(u6,'(2a)') 'poke_dscalar: Label is ',Label
!write(u6,'(a,es20.8)') 'poke_dscalar: Val is ',Val
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
!write(u6,'(a,i8)') 'poke_dscalar: indx is ',indx
!----------------------------------------------------------------------*
! Save data in buffer.                                                 *
!----------------------------------------------------------------------*
if (indx == -1) then
  if (ds_no >= size(ds_value)) then
    call SysAbendMsg('Poke_dScalar','Too many fields','Increase nTabDS and recompile')
  end if
  ds_no = ds_no+1
  indx = ds_no
end if
ds_label(indx) = Label
ds_value(indx) = val
!write(u6,'(a,i8)') 'poke_dscalar: ds_no is ',ds_no
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
return

end subroutine Poke_dScalar
