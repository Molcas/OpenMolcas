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
!               2023, Ignacio Fdez. Galvan                             *
!***********************************************************************

! This module calculates an external force applied to the system.

subroutine extf(ireturn)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, auToN
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
integer(kind=iwp) :: atomNumberx3, cnt, efatom1, efatom2, efatom3, efatom4, i, LuSpool, nCent, nRoots, nsAtom
real(kind=wp) :: Bt(3,4), dBt(3,4,3,4), Dum, efmodul, efmodulAU, fourAtoms(3,4), gau_sigma, gau_t0, mdtime, norm, time_scaling
logical(kind=iwp) :: bending, Found, gaussian_force, ldB, linear, lWarn, lWrite, torsional
character(len=180) :: Key, Line
character(len=8) :: Label
real(kind=wp), allocatable :: coord(:,:), ExtGrad(:,:), gradient(:,:), modgrad(:,:)
real(kind=wp), parameter :: nnewt = auToN*1.0e9_wp
integer(kind=iwp), external :: isfreeunit, Read_Grad
character(len=180), external :: Get_Ln

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

efmodul = Zero
linear = .false.
bending = .false.
torsional = .false.
gaussian_force = .false.

rewind(LuSpool)
call RdNLst(LuSpool,'extf')
do
  Key = Get_Ln(LuSpool)
  Line = Key
  call UpCase(Line)
  if (Line(1:3) == 'END') then
    exit
  else if (Line(1:4) == 'MODU') then
    Line = Get_Ln(LuSpool)
    call Get_F1(1,efmodul)
    write(u6,*) 'Force:',efmodul,' nN'
  else if (Line(1:4) == 'LINE') then
    !>>> LINEAR FORCE <<<
    linear = .true.
    write(u6,*) 'Linear forces between two atoms selected'
    Line = Get_Ln(LuSpool)
    call Get_I1(1,efatom1)
    call Get_I1(2,efatom2)
    write(u6,*) 'atom1:',efatom1
    write(u6,*) 'atom2:',efatom2
  else if (Line(1:4) == 'BEND') then
    !>>> BENDING FORCE <<<
    bending = .true.
    write(u6,*) 'Bending force on an angle selected'
    Line = Get_Ln(LuSpool)
    call Get_I1(1,efatom1)
    call Get_I1(2,efatom2)
    call Get_I1(3,efatom3)
    write(u6,*) 'atom1:',efatom1
    write(u6,*) 'atom2:',efatom2
    write(u6,*) 'atom3:',efatom3
  else if (Line(1:4) == 'TORS') then
    !>>> TORSIONAL FORCE <<<
    torsional = .true.
    write(u6,*) 'Torsional force on a dihedral selected'
    Line = Get_Ln(LuSpool)
    call Get_I1(1,efatom1)
    call Get_I1(2,efatom2)
    call Get_I1(3,efatom3)
    call Get_I1(4,efatom4)
    write(u6,*) 'atom1:',efatom1
    write(u6,*) 'atom2:',efatom2
    write(u6,*) 'atom3:',efatom3
    write(u6,*) 'atom4:',efatom4
  else if (Line(1:4) == 'GAUS') then
    gaussian_force = .true.
    Line = Get_Ln(LuSpool)
    call Get_F1(1,gau_t0)
    call Get_F1(2,gau_sigma)
  else
    call WarningMessage(2,'Error in EXTF input: Unknown keyword '//trim(Line))
    call Quit_OnUserError()
  end if
end do

close(LuSpool)

cnt = 0
if (linear) cnt = cnt+1
if (bending) cnt = cnt+1
if (torsional) cnt = cnt+1
if (cnt == 0) then
  call WarningMessage(2,'Error in FALSE input: One of LINE, BEND, TORS must be given')
  call Quit_OnUserError()
else if (cnt > 1) then
  call WarningMessage(2,'Error in FALSE input: Only one of LINE, BEND, TORS can be given')
  call Quit_OnUserError()
end if

write(u6,*) 'Coordinates Found:'
do i=1,nsAtom
  write(u6,*) i,coord(:,i)
end do
write(u6,*)

write(u6,*) 'Gradient Found:'
do i=1,nsAtom
  write(u6,*) i,gradient(:,i)
end do
write(u6,*)

! from nN to atomic units
efmodulAU = efmodul/nnewt

! initializations
ExtGrad(:,:) = Zero
fourAtoms(:,:) = Zero
Bt(:,:) = Zero
dBt(:,:,:,:) = Zero ! second derivative
lWrite = .true. ! write something in output file
lWarn = .true. ! write warnings in output file
ldB = .false. ! this will NOT calculate the second derivative

!>>>>>>>>>>>>>>>>>>>>> LINEAR CODE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if (linear) then

  if (efmodul > Zero) then
    write(u6,*) 'Compression force'
  else
    write(u6,*) 'Extension force'
  end if

  fourAtoms(:,1) = coord(:,efatom1)
  fourAtoms(:,2) = coord(:,efatom2)

  nCent = 2
  Label = 'Line'
  call Strtch(fourAtoms,nCent,Dum,Bt,lWrite,Label,dBt,ldB)

  write(u6,*) "Bt vector:"
  do i=1,nCent
    write(u6,*) i,Bt(:,i)
  end do

  norm = sqrt(sum(Bt(:,1:nCent)**2))
  write(u6,*) 'Bt norm before normalization:',norm
  Bt(:,1:nCent) = Bt(:,1:nCent)/norm

  ExtGrad(:,efatom1) = Bt(:,1)*efmodulAU
  ExtGrad(:,efatom2) = Bt(:,2)*efmodulAU

end if
!>>>>>>>>>>>>>>>>>>>>> end of linear code <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!>>>>>>>>>>>>>>>>>>>>> BENDING CODE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if (bending) then

  if (efmodul > Zero) then
    write(u6,*) 'Closing force'
  else
    write(u6,*) 'Opening force'
  end if

  fourAtoms(:,1) = coord(:,efatom1)
  fourAtoms(:,2) = coord(:,efatom2)
  fourAtoms(:,3) = coord(:,efatom3)

  nCent = 3
  Label = 'Angle'
  call Bend(fourAtoms,nCent,Dum,Bt,lWrite,lWarn,Label,dBt,ldB)

  write(u6,*) "Bt vector:"
  do i=1,nCent
    write(u6,*) i,Bt(:,i)
  end do

  norm = sqrt(sum(Bt(:,1:nCent)**2))
  write(u6,*) 'Bt norm before normalization:',norm
  Bt(:,1:nCent) = Bt(:,1:nCent)/norm

  ExtGrad(:,efatom1) = Bt(:,1)*efmodulAU
  ExtGrad(:,efatom2) = Bt(:,2)*efmodulAU
  ExtGrad(:,efatom3) = Bt(:,3)*efmodulAU

end if
!>>>>>>>>>>>>>>>>>>>>> end of bending code <<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!>>>>>>>>>>>>>>>>>>>>> TORSIONAL CODE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if (torsional) then

  if (efmodul > Zero) then
    write(u6,*) 'Closing force'
  else
    write(u6,*) 'Opening force'
  end if

  fourAtoms(:,1) = coord(:,efatom1)
  fourAtoms(:,2) = coord(:,efatom2)
  fourAtoms(:,3) = coord(:,efatom3)
  fourAtoms(:,4) = coord(:,efatom4)

  nCent = 4
  Label = 'Dihedral'
  call Trsn(fourAtoms,nCent,Dum,Bt,lWrite,lWarn,Label,dBt,ldB)

  write(u6,*) "Bt vector:"
  do i=1,nCent
    write(u6,*) i,Bt(:,i)
  end do

  norm = sqrt(sum(Bt(:,1:nCent)**2))
  write(u6,*) 'Bt norm before normalization:',norm
  Bt(:,1:nCent) = Bt(:,1:nCent)/norm

  ExtGrad(:,efatom1) = Bt(:,1)*efmodulAU
  ExtGrad(:,efatom2) = Bt(:,2)*efmodulAU
  ExtGrad(:,efatom3) = Bt(:,3)*efmodulAU
  ExtGrad(:,efatom4) = Bt(:,4)*efmodulAU

end if
!>>>>>>>>>>>>>>>>>>>>> end of torsional code <<<<<<<<<<<<<<<<<<<<<<<<<<<

!>>>>>>>>>>>>>>>>>>>>> GAUSSIAN CODE <<<<_<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if (gaussian_force) then
  ! AV: This is called before dynamix, so at fist step it would not find the time value in the runfile
  mdtime = Zero
  call Qpg_dScalar('MD_Time',Found)
  if (Found) call Get_dScalar('MD_Time',mdtime)

  write(u6,*) 'Gaussian shaped force'
  write(u6,*) 'sigma',gau_sigma,'t0',gau_t0
  time_scaling = exp(-(mdtime-gau_t0)**2/(2*(gau_sigma**2)))
  write(u6,*) 'gaussian:',mdtime,time_scaling

  ExtGrad(:,:) = ExtGrad(:,:)*time_scaling

end if
!>>>>>>>>>>>>>>>>>>>>> end of gaussian code <<<<<<<<<<<<<<<<<<<<<<<<<<<<

write(u6,*)
write(u6,*) 'External Force'
do i=1,nsAtom
  write(u6,*) i,-ExtGrad(:,i)
end do
write(u6,*)

! Creating the final modified external gradient vector (modgrad)
modgrad(:,:) = gradient(:,:)+ExtGrad(:,:)

write(u6,*)
write(u6,*) 'Gradient after force application:'
do i=1,nsAtom
  write(u6,*) i,modgrad(:,i)
end do
write(u6,*)

call Put_dArray('GRAD',modgrad,atomNumberx3)

! Modify the GRADS file too
call Get_iScalar('Number of roots',nRoots)
do i=1,nRoots
  if (Read_Grad(gradient,3*nsAtom,i,0,0) == 1) then
    modgrad(:,:) = gradient(:,:)+ExtGrad(:,:)
    call Store_Grad(modgrad,3*nsAtom,i,0,0)
  end if
end do

call mma_deallocate(gradient)
call mma_deallocate(ExtGrad)
call mma_deallocate(modgrad)
call mma_deallocate(coord)

ireturn = 0

return

end subroutine extf
