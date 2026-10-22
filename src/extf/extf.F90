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

use spool, only: Spoolinp
use stdalloc, only: mma_allocate, mma_deallocate
use Symmetry_Info, only: Symmetry_Info_Dmp, Symmetry_Info_Free, Symmetry_Info_Get, VarR
use Constants, only: Zero, One, Pi, Angstrom, auToN, deg2rad
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
integer(kind=iwp) :: atomNumberx3, cnt, efatom1, efatom2, efatom3, efatom4, i, istatus, k, l, LuSpool, nback, nCent, nfor, nRoots, &
                     nsAtom, nSym, ntr
real(kind=wp) :: Bt(3,4), com(3), dBt(3,4,3,4), dperp, dperpmax, Dum, efmodul, fourAtoms(3,4), ftot, gau_sigma, gau_t0, mdtime, &
                 nhat(3), norm, ovlp, proj, rmaxAng, rnorm, time_scaling, totmass, vnorm, wf_phi, wf_theta
logical(kind=iwp) :: bending, Found, gaussian_force, keep_rot, ldB, linear, lWarn, lWrite, torsional, wall, wall_vdw
character(len=180) :: Key, Line
character(len=8) :: Label
real(kind=wp), allocatable :: charge(:), coord(:,:), ExtGrad(:,:), gradient(:,:), mass(:), modgrad(:,:), trvec(:,:,:)
real(kind=wp), parameter :: nnewt = auToN*1.0e9_wp
integer(kind=iwp), external :: isfreeunit, Read_Grad
real(kind=wp), external :: vdwRad
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
wall = .false.
wall_vdw = .false.
keep_rot = .false.
wf_theta = Zero
wf_phi = Zero
nhat(:) = Zero

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
  else if (Line(1:4) == 'WALL') then
    !>>> WALL FORCE (BALL MILLING) <<<
    wall = .true.
    write(u6,*) 'Wall-type (ball-milling) compression forces selected'
    Line = Get_Ln(LuSpool)
    call UpCase(Line)
    read(Line,*,iostat=istatus) nhat(1),nhat(2),nhat(3)
    if (istatus == 0) then
      write(u6,*) 'n:',nhat(:)
    else
      read(Line,*) wf_theta,wf_phi
      if (index(Line,'DEGREE') > 0) then
        write(u6,*) 'theta (deg):',wf_theta
        write(u6,*) 'phi (deg):',wf_phi
        wf_theta = wf_theta*deg2rad
        wf_phi = wf_phi*deg2rad
      else
        write(u6,*) 'theta (rad):',wf_theta
        write(u6,*) 'phi (rad):',wf_phi
      end if
    end if
  else if (Line(1:4) == 'VDWR') then
    ! >>> NEED TO ADD VAN DER WAALS RADII TO CALCULATE EFFECTIVE RADIUS <<<
    wall_vdw = .true.
    write(u6,*) 'van der Waals radii included in the effective radius'
  else if (Line(1:4) == 'ROTA') then
    !>>> keep the overall rotation (do not remove the net torque) <<<
    keep_rot = .true.
    write(u6,*) 'Overall rotation kept: net rotational wall force will NOT be removed'
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
if (wall) cnt = cnt+1
if (cnt == 0) then
  call WarningMessage(2,'Error in FALSE input: One of LINE, BEND, TORS, WALL must be given')
  call Quit_OnUserError()
else if (cnt > 1) then
  call WarningMessage(2,'Error in FALSE input: Only one of LINE, BEND, TORS, WALL can be given')
  call Quit_OnUserError()
end if

if (wall) then
  write(u6,*) 'Pressure:',efmodul,' GPa'
else
  write(u6,*) 'Force:',efmodul,' nN'
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
efmodul = efmodul/nnewt

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

  ExtGrad(:,efatom1) = Bt(:,1)*efmodul
  ExtGrad(:,efatom2) = Bt(:,2)*efmodul

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

  ExtGrad(:,efatom1) = Bt(:,1)*efmodul
  ExtGrad(:,efatom2) = Bt(:,2)*efmodul
  ExtGrad(:,efatom3) = Bt(:,3)*efmodul

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

  ExtGrad(:,efatom1) = Bt(:,1)*efmodul
  ExtGrad(:,efatom2) = Bt(:,2)*efmodul
  ExtGrad(:,efatom3) = Bt(:,3)*efmodul
  ExtGrad(:,efatom4) = Bt(:,4)*efmodul

end if
!>>>>>>>>>>>>>>>>>>>>> end of torsional code <<<<<<<<<<<<<<<<<<<<<<<<<<<

!>>>>>>>>>>>>>>>>>>>>> WALL FORCES CODE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Wall-type compression forces simulating ball-milling impacts
! (R. De Armas, M. Temprado, L. M. Frutos, J. Chem. Theory Comput. 2025, 21, 10353)
! Two virtual walls perpendicular to the direction n(theta,phi) compress
! the molecule. The total force is F = P*pi*R(theta,phi)^2, with R the
! effective radius of the molecular cross section perpendicular to n,
! and each wall exerts it on the atoms of its hemisphere. The geometry is
! shifted to its centre of mass and the molecule is divided in two halves
! by the plane perpendicular to n and passing through the centre of mass:
! every atom of a hemisphere is pushed towards that dividing plane with
! the same force, ftot divided by the number of atoms of its hemisphere.
if (wall) then

  ! WALL forces act on all atoms (centre of mass, effective radius, force distribution),
  ! so a symmetry-reduced list of unique atoms would give wrong results. Require no symmetry (C1).
  call Get_iScalar('nSym',nSym)
  if (nSym > 1) then
    call WarningMessage(2,'Error in EXTF: WALL (ball-milling) forces require no symmetry; please run in C1')
    call Quit_OnUserError()
  end if

  call mma_allocate(mass,nsAtom,label='mass')
  call mma_allocate(charge,nsAtom,label='charge')
  call Get_Mass(mass,nsAtom)
  call Get_dArray('Nuclear charge',charge,nsAtom)

  if (wall_vdw) then
    do i=1,nsAtom
      if (vdwRad(nint(charge(i))) <= Zero) then
        call WarningMessage(2,'Error in EXTF: van der Waals radius not available for some atom')
        call Quit_OnUserError()
      end if
    end do
  end if

  ! shift the geometry to its centre of mass
  totmass = sum(mass(:))
  com(1) = sum(coord(1,:)*mass(:))/totmass
  com(2) = sum(coord(2,:)*mass(:))/totmass
  com(3) = sum(coord(3,:)*mass(:))/totmass
  do i=1,nsAtom
    coord(:,i) = coord(:,i)-com(:)
  end do
  write(u6,*) 'Centre of mass (bohr):',com(:)

  ! direction of the compression, n(theta,phi)
  rnorm = sqrt(sum(nhat(:)**2))
  if (rnorm > Zero) then
    nhat(:) = nhat(:)/rnorm
  else
    nhat(1) = sin(wf_theta)*cos(wf_phi)
    nhat(2) = sin(wf_theta)*sin(wf_phi)
    nhat(3) = cos(wf_theta)
  end if
  write(u6,*) 'Compression direction:',nhat(:)

  ! effective radius: largest distance from the compression axis to the
  ! atomic surfaces. If requested, the van der Waals radius is added to
  ! every atom before taking the maximum (the outter atom is not
  ! necessarily the one further located from the axis)
  dperpmax = Zero
  do i=1,nsAtom
    proj = sum(coord(:,i)*nhat(:))
    dperp = sqrt(sum((coord(:,i)-proj*nhat(:))**2))
    if (wall_vdw) dperp = dperp+vdwRad(nint(charge(i)))/Angstrom
    if (dperp > dperpmax) dperpmax = dperp
  end do
  rmaxAng = dperpmax*Angstrom
  write(u6,*) 'Effective radius:',rmaxAng,' angstrom'

  ! total compression force: F = P*pi*R(theta,phi)^2
  ! (1 GPa * 1 Angstrom^2 = 0.01 nN)
  ftot = efmodul*0.01_wp*Pi*rmaxAng**2
  write(u6,*) 'Total wall force:',ftot*nnewt,' nN'

  call mma_deallocate(mass)
  call mma_deallocate(charge)

  ! distribute the force among the atoms: each wall exerts the total
  ! force ftot on the atoms of its hemisphere, shared equally among them.
  ! An atom in the "forward" hemisphere (r.n > 0) is pushed along -n
  ! with |F| = ftot/nfor, and an atom in the "backward" hemisphere
  ! (r.n < 0) is pushed along +n with |F| = ftot/nback, so that the
  ! moduli of the forces of each hemisphere sum up the wall force
  ! and the net force is zero.
  nfor = 0
  nback = 0
  do i=1,nsAtom
    proj = sum(coord(:,i)*nhat(:))
    if (proj > Zero) then
      nfor = nfor+1
    else if (proj < Zero) then
      nback = nback+1
    end if
  end do
  if ((nfor == 0) .or. (nback == 0)) then
    call WarningMessage(2,'Error in EXTF: all atoms lie on one side of the dividing plane, WALL forces are undefined')
    call Quit_OnUserError()
  end if
  do i=1,nsAtom
    proj = sum(coord(:,i)*nhat(:))
    if (proj > Zero) then
      ExtGrad(:,i) = (ftot/real(nfor,kind=wp))*nhat(:)
    else if (proj < Zero) then
      ExtGrad(:,i) = -(ftot/real(nback,kind=wp))*nhat(:)
    end if
  end do

  ! Remove the net force (always) and the net torque (unless ROTA is given)
  ! that the wall forces exert on the molecule in a net way, so that the
  ! compression forces does not provoke an artificial translation of the centre of mass
  ! or an overall rotation during a dynamics. This is done by projecting
  ! ExtGrad out of the space spanned by translations and rotations
  ntr = 3
  if (.not. keep_rot) then
    ntr = 6
    call Symmetry_Info_Get()
    VarR = .true.
    call Symmetry_Info_Dmp()
    call Symmetry_Info_Free()
  end if
  call mma_allocate(trvec,3,nsAtom,ntr,label='trvec')
  trvec(:,:,:) = Zero
  trvec(1,:,1) = One
  trvec(2,:,2) = One
  trvec(3,:,3) = One
  if (.not. keep_rot) then
    do i=1,nsAtom
      trvec(:,i,4) = [Zero,-coord(3,i),coord(2,i)]
      trvec(:,i,5) = [coord(3,i),Zero,-coord(1,i)]
      trvec(:,i,6) = [-coord(2,i),coord(1,i),Zero]
    end do
  end if
  ! orthonormalise, dropping null vectors such as
  ! the degenerate rotation of a linear molecule
  do k=1,ntr
    do l=1,k-1
      ovlp = sum(trvec(:,:,l)*trvec(:,:,k))
      trvec(:,:,k) = trvec(:,:,k)-ovlp*trvec(:,:,l)
    end do
    vnorm = sqrt(sum(trvec(:,:,k)**2))
    if (vnorm > 1.0e-8_wp) then
      trvec(:,:,k) = trvec(:,:,k)/vnorm
    else
      trvec(:,:,k) = Zero
    end if
  end do
  ! project ExtGrad onto the orthogonal complement
  do k=1,ntr
    ovlp = sum(trvec(:,:,k)*ExtGrad(:,:))
    ExtGrad(:,:) = ExtGrad(:,:)-ovlp*trvec(:,:,k)
  end do
  call mma_deallocate(trvec)

end if
!>>>>>>>>>>>>>>>>>>>>> end of wall forces code <<<<<<<<<<<<<<<<<<<<<<<<<

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
