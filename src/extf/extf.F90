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
integer(kind=iwp) :: i, j, atomNumberx3, nsAtom, efatom1, efatom2, ext, LuSpool
integer(kind=iwp) :: efatom3, efatom4
real(kind=wp) :: efmodul, efmodulAU, norm, posvect12(3)
real(kind=wp), allocatable :: gradient(:,:), ExtGrad(:,:), modgrad(:,:), coord(:,:)
real(kind=wp), parameter :: nnewt = auToN*1.0e9_wp
logical(kind=iwp) :: linear, Found
character(len=180) :: Key, Line
character(len=180), external :: Get_Ln
integer(kind=iwp) :: isfreeunit
integer(kind=iwp) :: nCent
real(kind=wp) :: Tau, gau_sigma, gau_t0, mdtime, time_scaling
real(kind=wp) :: Bt(3,4), dBt(3,4,3,4), fourAtoms(3,4)
logical(kind=iwp) :: torsional, lWrite, lWarn, ldB, gaussian_force
character(len=180) :: Label

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

linear=.false.
torsional=.false.

rewind(LuSpool)
call RdNLst(LuSpool,'extf')
do
  Key = Get_Ln(LuSpool)
  Line = Key
  call UpCase(Line)
  if (Line(1:3) == 'END') exit
  if (Line(1:4) == 'LINE') then
    !>>> LINEAR FORCe <<<
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
!>>>>>>>>>>>>>>>>>>>> LINEAR END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 if (Line(1:4) == 'TORS') then
!>>>>>>>>>>>>>>>>>>>> torsional FORCE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    torsional = .true.
    write(u6,*) 'Torsional force on a dihedral selected'
    Line = Get_Ln(LuSpool)
    call Get_I1(1,efatom1)
    Line = Get_Ln(LuSpool)
    call Get_I1(1,efatom2)
    Line = Get_Ln(LuSpool)
    call Get_I1(1,efatom3)
    Line = Get_Ln(LuSpool)
    call Get_I1(1,efatom4)
    Line = Get_Ln(LuSpool)
    call Get_F1(1,efmodul)
    Line = Get_Ln(LuSpool)
    call Get_I1(1,ext)
    write(u6,*) 'atom1:',efatom1
    write(u6,*) 'atom2:',efatom2
    write(u6,*) 'atom3:',efatom3
    write(u6,*) 'atom4:',efatom4
    write(u6,*) 'Force:',efmodul,' nN'
    if (ext == 1) then
      write(u6,*) 'Closing force'
    else
      write(u6,*) 'Opening force'
    end if
 end if

 if (Line(1:4) == 'GAUS') then
    gaussian_force = .true.
    Line = Get_Ln(LuSpool)
    call Get_F1(1,gau_sigma)
    Line = Get_Ln(LuSpool)
    call Get_F1(1,gau_t0)
 end if
!
!>>>>>>>>>>>>>>>>>>>> torsional END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end do

write(u6,*) 'Coordinates Found:'
do i=1,nsAtom
write(u6,*) i, coord(:,i)
end do
write(u6,*)

write(u6,*) 'Gradient Found:'
do i=1,nsAtom
write(u6,*) i, gradient(:,i)
end do
write(u6,*)

!>>>>>>>>>>>>>>>>>>>>> LINEAR CODE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if (linear) then

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



end if
!>>>>>>>>>>>>>>>>>>>>> end of linear code <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if (torsional) then
  write(u6, *) 'Torsional force applied'
  ! from nN to atomic units
  efmodulAU = efmodul/nnewt

  ! initialize the external gradient vector
  ExtGrad(:,:) = Zero

  fourAtoms(:,:) = Zero ! initialization matrix
  fourAtoms(:,1) = coord(:,efatom1)
  fourAtoms(:,2) = coord(:,efatom2)
  fourAtoms(:,3) = coord(:,efatom3)
  fourAtoms(:,4) = coord(:,efatom4)

  nCent = 4
  Tau = 0.0
  Bt(:,:) = Zero ! initialization matrix
  lWrite = .true. ! write something in output file
  lWarn = .true. ! write warnings in output file
  Label = 'Dihedral'
  dBt(:,:,:,:) = Zero ! second derivative
  ldB = .false. ! this will NOT calculate the second derivative
  Call Trsn(fourAtoms,nCent,Tau,Bt,lWrite,lWarn,Label,dBt,ldB)

  write(u6,*) "Bt vector:"
  do i=1,4
    write(u6,*) i, Bt(:,i)
  end do

  norm = 0.0
  do i=1,4
    do j=1,3
      norm = norm + Bt(j,i)**2
    end do
  end do
  norm = sqrt(norm)
  write(u6, *) 'Bt norm before normalization:', norm


  do i=1,4
    do j=1,3
      Bt(j,i) = Bt(j,i)/norm
    end do
  end do

  ExtGrad(:,efatom1) = Bt(:,1)*efmodulAU
  ExtGrad(:,efatom2) = Bt(:,2)*efmodulAU
  ExtGrad(:,efatom3) = Bt(:,3)*efmodulAU
  ExtGrad(:,efatom4) = Bt(:,4)*efmodulAU

end if

if (gaussian_force) then

  ! AV: This is called before dynamix, so at fist step it would not find
  ! the time value in the runfile
  call Qpg_dScalar('MD_Time',Found)
  if (found) then
    call Get_dScalar('MD_Time', mdtime)
  else
    mdtime = 0.0
  end if

  write(u6,*) 'Gaussian shaped force'
  write(u6,*) 'sigma', gau_sigma, 't0', gau_t0
  time_scaling = exp(-(mdtime-gau_t0)**2/(2*gau_sigma)**2)
  write(u6,*) 'gaussian:', mdtime, time_scaling

  ExtGrad(:,:) = ExtGrad(:,:) * time_scaling

end if

! Creating the final modified external gradient vector (modgrad)
modgrad(:,:) = gradient(:,:)+ExtGrad(:,:)

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
