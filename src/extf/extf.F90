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
! Copyright (C) 2015, Luis Manuel Frutos                               *
!               2015, Alessio Valentini                                *
!***********************************************************************

! This module calculates an external force applied to the system.

subroutine extf(ireturn)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, auToN
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
integer(kind=iwp) :: i, atomNumberx3, nsAtom, efatom1, efatom2, ext, LuSpool
real(kind=wp) :: efmodul, efmodulAU, norm, posvect12(3)
real(kind=wp), allocatable :: gradient(:,:), ExtGrad(:,:), modgrad(:,:), coord(:,:)
real(kind=wp), parameter :: nnewt = auToN*1.0e9_wp
logical(kind=iwp) :: linear
character(len=180) :: Key, Line
character(len=180), external :: Get_Ln
integer(kind=iwp) :: isfreeunit

! get initial values

!call Mem_Info('EXTF')

call Get_iScalar('Unique atoms',nsAtom)

atomNumberx3 = 3*nsAtom
call mma_allocate(gradient,3,nsAtom,label='gradient')
call mma_allocate(ExtGrad,3,nsAtom,label='ExtGrad')
call mma_allocate(modgrad,3,nsAtom,label='modgrad')
call mma_allocate(coord,3,nsAtom,label='coord')

call Get_dArray('GRAD',gradient,atomNumberx3)
call Get_dArray('Unique Coordinates',coord,atomNumberx3)
! read molcas input

LuSpool = isfreeunit(21)
call SpoolInp(LuSpool)

rewind(LuSpool)
call RdNLst(LuSpool,'extf')
do
  Key = Get_Ln(LuSpool)
  Line = Key
  call UpCase(Line)
  if (Line(1:3) == 'END') exit
  if (Line(1:4) == 'LINE') then
!>>>>>>>>>>>>>>>>>>>> FORCe <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    linear = .true.
    write(u6,*) 'Linear forces between two atoms selected'
    Line = Get_Ln(LuSpool)
    call Get_I1(1,efatom1)
    Line = Get_Ln(LuSpool)
    call Get_I1(1,efatom2)
    Line = Get_Ln(LuSpool)
    call Get_F1(1,efmodul)
    Line = Get_Ln(LuSpool)
    call Get_I1(1,ext)
    write(u6,*) 'atom1:',efatom1
    write(u6,*) 'atom2:',efatom2
    write(u6,*) 'Force:',efmodul,' nN'
    if (ext == 1) then
      write(u6,*) 'Compression force'
    else
      write(u6,*) 'Extension force'
    end if
  end if
!>>>>>>>>>>>>>>>>>>>> END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end do
!>>>>>>>>>>>>>>>>>>>>> LINEAR CODE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if (linear) then
  write(u6,*) 'Gradient Found:'
  do i=1,nsAtom
    write(u6,*) i,gradient(:,i)
  end do
  write(u6,*)

  ! from nN to atomic units
  efmodulAU = efmodul/nnewt

  ! initialize the external gradient vector
  ExtGrad(:,:) = Zero

  ! posvect is the position vector from atom2 respect atom1
  posvect12(:) = coord(:,efatom2)-coord(:,efatom1)

  ! creating "norm": the norm of the posvect12 vector
  norm = sqrt(posvect12(1)**2+posvect12(2)**2+posvect12(3)**2)
  posvect12(:) = posvect12(:)/norm

  ExtGrad(:,efatom1) = posvect12(:)*efmodulAU
  ExtGrad(:,efatom2) = -posvect12(:)*efmodulAU
  ! Extension gradient vector created.

  ! Checking if it is a compression or extension force. In the former case
  ! the direction of the ExtGrad vector will be changed
  if (ext == 1) then
    ExtGrad(:,efatom1) = -ExtGrad(:,efatom1)
    ExtGrad(:,efatom2) = -ExtGrad(:,efatom2)
  end if

  write(u6,*)
  write(u6,*) 'External Force'
  do i=1,nsAtom
    write(u6,*) i,-ExtGrad(:,i)
  end do
  write(u6,*)

  ! Creating the final modified external gradient vector (modgrad)
  modgrad(:,:) = gradient(:,:)+ExtGrad(:,:)

end if
!>>>>>>>>>>>>>>>>>>>>> end of linear code <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

write(u6,*)
write(u6,*) 'Gradient after force application:'
do i=1,nsAtom
  write(u6,*) i,modgrad(:,i)
end do
write(u6,*)

call Put_dArray('GRAD',modgrad,atomNumberx3)

close(LuSpool)

call mma_deallocate(gradient)
call mma_deallocate(ExtGrad)
call mma_deallocate(modgrad)
call mma_deallocate(coord)

ireturn = 0

return

end subroutine extf
