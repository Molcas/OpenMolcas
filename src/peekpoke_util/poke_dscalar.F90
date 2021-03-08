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
!> @param[in] Data  Data to put on runfile
!***********************************************************************

subroutine Poke_dScalar(Label,data)

implicit none
#include "pp_ds_info.fh"
!----------------------------------------------------------------------*
! Arguments                                                            *
!----------------------------------------------------------------------*
character*(*) Label
real*8 data
!----------------------------------------------------------------------*
! Define local variables                                               *
!----------------------------------------------------------------------*
integer indx
integer i

!----------------------------------------------------------------------*
! Initialize local variables                                           *
!----------------------------------------------------------------------*
!write(6,'(2a)') 'poke_dscalar: Label is ',Label
!write(6,'(a,e20.8)') 'poke_dscalar: Data is ',Data
!----------------------------------------------------------------------*
! Locate item                                                          *
!----------------------------------------------------------------------*
indx = -1
do i=1,ds_no
  if (ds_label(i) == Label) indx = i
end do
!write(6,'(a,i8)') 'poke_dscalar: indx is ',indx
!----------------------------------------------------------------------*
! Save data in buffer.                                                 *
!----------------------------------------------------------------------*
if (indx == -1) then
  if (ds_no >= nTabDS) then
    call SysAbendMsg('Poke_dScalar','Too many fields','Increase nTabDS and recompile')
  end if
  ds_no = ds_no+1
  indx = ds_no
end if
ds_label(indx) = Label
ds_value(indx) = data
!write(6,'(a,i8)') 'poke_dscalar: ds_no is ',ds_no
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
return

end subroutine Poke_dScalar
