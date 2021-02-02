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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine dumps characters in xml format.                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University, Sweden                                     *
!                                                                      *
!***********************************************************************

subroutine xml_cDump(TagName,Appear,Units,Level,Content,nx,ny)

use Definitions, only: iwp

implicit none
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
character(len=*), intent(in) :: TagName, Appear, Units, Content(*)
integer(kind=iwp), intent(in) :: nx, ny, Level
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
integer(kind=iwp) :: ix, iy, ind, opt
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
if ((ny == 1).and.(nx < 5)) then
  call xml_cDumpa(TagName,len(TagName),Appear,len(Appear),Units,len(Units),Level,nx,ny,0)
  do ix=1,nx
    call xml_cDumpb(Content(ix),len(Content(ix)),0)
  end do
else
  call xml_cDumpa(TagName,len(TagName),Appear,len(Appear),Units,len(Units),Level,nx,ny,1)
  do iy=1,ny
    do ix=1,nx
      ind = (ix-1)*ny+iy
      if ((mod(ix,10) == 0).or.(ix == nx)) then
        opt = 1
      else
        opt = 0
      end if
      call xml_cDumpb(Content(ind),len(Content(ind)),opt)
    end do
  end do
end if
call xml_cDumpc(TagName,len(TagName))
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
return

end subroutine xml_cDump
