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
integer(kind=iwp) :: atomNumberx3, cnt, efatom1, efatom2, efatom3, efatom4, i, istatus, iwmax, k, l, LuSpool, nCent, nRoots, &
                     nsAtom, nSym, ntr
real(kind=wp) :: Bt(3,4), com(3), dBt(3,4,3,4), dperp, dperpmax, Dum, efmodul, fourAtoms(3,4), ftot, gau_sigma, gau_t0, mdtime, &
                 nhat(3), norm, ovlp, proj, rmaxAng, rnorm, time_scaling, totmass, vdwtab(118), vnorm, wf_phi, wf_theta, wproj
logical(kind=iwp) :: bending, Found, gaussian_force, keep_rot, ldB, linear, lWarn, lWrite, torsional, wall, wall_vdw
character(len=180) :: Key, Line
character(len=8) :: Label
real(kind=wp), allocatable :: charge(:), coord(:,:), ExtGrad(:,:), gradient(:,:), mass(:), modgrad(:,:), trvec(:,:,:)
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
! and it is distributed among all the atoms. The geometry is shifted to its
! centre of mass and the molecule is divided in two halves by the plane
! perpendicular to n and passing through the centre of mass:
! each atom is pushed in the direction of that dividing plane,
! weighted by the projection of its position unit vector on n.
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

  ! van der Waals radii in Angstrom (zero means not available). Data taken from (except transition metals:
  !     "Consistent van der Waals Radii for the Whole Main Group"
  !     M. Mantina, A. C. Chamberlin, R. Valero, C. J. Cramer, D. G. Truhlar
  !     J. Phys. Chem. A 2009, 113, 19, 5806-5812
  !
  ! van der Waals radii in Angstrom for transition metals:
  !     "A cartography of the van der Waals territories"
  !     S. Alvarez
  !     Dalton Trans., 2013, 42, 8617-8636

  vdwtab(:) = Zero

  vdwtab(1) = 1.10_wp  ! H
  vdwtab(2) = 1.40_wp  ! He

  vdwtab(3) = 1.81_wp  ! Li
  vdwtab(4) = 1.53_wp  ! Be
  vdwtab(5) = 1.92_wp  ! B
  vdwtab(6) = 1.70_wp  ! C
  vdwtab(7) = 1.55_wp  ! N
  vdwtab(8) = 1.52_wp  ! O
  vdwtab(9) = 1.47_wp  ! F
  vdwtab(10) = 1.54_wp  ! Ne

  vdwtab(11) = 2.27_wp  ! Na
  vdwtab(12) = 1.73_wp  ! Mg
  vdwtab(13) = 1.84_wp  ! Al
  vdwtab(14) = 2.10_wp  ! Si
  vdwtab(15) = 1.80_wp ! P
  vdwtab(16) = 1.80_wp ! S
  vdwtab(17) = 1.75_wp ! Cl
  vdwtab(18) = 1.88_wp ! Ar

  vdwtab(19) = 2.75_wp ! K
  vdwtab(20) = 2.31_wp ! Ca
      vdwtab(21) = 2.58_wp ! Sc
      vdwtab(22) = 2.46_wp ! Ti
      vdwtab(23) = 2.42_wp ! V
      vdwtab(24) = 2.45_wp ! Cr
      vdwtab(25) = 2.45_wp ! Mn
      vdwtab(26) = 2.44_wp ! Fe
      vdwtab(27) = 2.40_wp ! Co
      vdwtab(28) = 2.40_wp ! Ni
      vdwtab(29) = 2.38_wp ! Cu
      vdwtab(30) = 2.39_wp ! Zn
  vdwtab(31) = 1.87_wp ! Ga
  vdwtab(32) = 2.11_wp ! Ge
  vdwtab(33) = 1.85_wp ! As
  vdwtab(34) = 1.90_wp ! Se
  vdwtab(35) = 1.83_wp ! Br
  vdwtab(36) = 2.02_wp ! Kr


  vdwtab(37) = 3.03_wp ! Rb
  vdwtab(38) = 2.49_wp ! Sr
      vdwtab(39) = 2.75_wp ! Y
      vdwtab(40) = 2.52_wp ! Zr
      vdwtab(41) = 2.56_wp ! Nb
      vdwtab(42) = 2.45_wp ! Mo
      vdwtab(43) = 2.44_wp ! Tc
      vdwtab(44) = 2.46_wp ! Ru
      vdwtab(45) = 2.44_wp ! Rh
      vdwtab(46) = 2.15_wp ! Pd
      vdwtab(47) = 2.53_wp ! Ag
      vdwtab(48) = 2.49_wp ! Cd
  vdwtab(49) = 1.93_wp ! In
  vdwtab(50) = 2.17_wp ! Sn
  vdwtab(51) = 2.06_wp ! Sb
  vdwtab(52) = 2.06_wp ! Te
  vdwtab(53) = 1.98_wp ! I
  vdwtab(54) = 2.16_wp ! Xe

  if (wall_vdw) then
    do i=1,nsAtom
      if (vdwtab(nint(charge(i))) <= Zero) then
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

  ! effective radius: largest atomic distance to the compression axis
  ! (plus the van der Waals radius of that atom, if requested)
  dperpmax = Zero
  iwmax = 1
  do i=1,nsAtom
    proj = sum(coord(:,i)*nhat(:))
    dperp = sqrt(sum((coord(:,i)-proj*nhat(:))**2))
    if (dperp > dperpmax) then
      dperpmax = dperp
      iwmax = i
    end if
  end do
  rmaxAng = dperpmax*Angstrom
  if (wall_vdw) rmaxAng = rmaxAng+vdwtab(nint(charge(iwmax)))
  write(u6,*) 'Effective radius:',rmaxAng,' angstrom'

  ! total compression force: F = P*pi*R(theta,phi)^2
  ! (1 GPa * 1 Angstrom^2 = 0.01 nN)
  ftot = efmodul*0.01_wp*Pi*rmaxAng**2
  write(u6,*) 'Total wall force:',ftot*nnewt,' nN'

  ! distribute the force among the atoms: an atom in the "forward"
  ! hemisphere (r.n > 0) is pushed along -n, and an atom in the
  ! "backward" hemisphere (r.n < 0) is pushed along +n
  do i=1,nsAtom
    rnorm = sqrt(sum(coord(:,i)**2))
    if (rnorm > Zero) then
      wproj = sum(coord(:,i)*nhat(:))/rnorm
    else
      wproj = Zero
    end if
    ExtGrad(:,i) = wproj*(ftot/real(nsAtom,kind=wp))*nhat(:)
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

  call mma_deallocate(mass)
  call mma_deallocate(charge)

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
