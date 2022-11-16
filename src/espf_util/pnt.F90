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
! Copyright (C) Christophe Chipot                                      *
!               Janos G. Angyan                                        *
!***********************************************************************

subroutine PNT(IOut,nAtoms,CO,IRMax,DeltaR,IAn,nPts,Grid,IsMM,Process)
! THIS SUBROUTINE CREATES A GRID OF POINTS AROUND A MOLECULE.
!
! THE POINTS ARE INITIALLY SELECTED IN A CUBE, THEN RE-SCALED TO THE
! SIZE OF THE MOLECULE.
!
! THE POINTS LYING WITHIN THE VAN DER WAALS ENVELOPE OF THE
! MOLECULE ARE EXCLUDED, BECAUSE OF THE PENETRATION EFFECTS.
! THIS PROCEDURE IS BASED ON THE SCALED VDW RADII.
! THE SCALING FACTOR : VDWFact
!
! THE POINTS LYING OUTSIDE THE EXCLUSION RADIUS ARE EXCLUDED TOO.
! THE EXCLUSION RADIUS  :  r_max
! THE GRID STEP         :  delta_r
!
! This subroutine comes from the Grid_3.1.2 program (C.CHIPOT, J.G.ANGYAN,
! Universite' Henri Poincare' - Nancy 1, France)

implicit real*8(A-H,O-Z)
logical Retain1, Retain2, Process
integer IsMM(nAtoms)
dimension CO(3,nAtoms), IAn(nAtoms), Grid(3,*), RVDW(0:54)
save Zero, Two, Huge, RVDW
data Zero/0.0d0/,Two/2.0d0/,Huge/1.0d8/,RVDW/2.00d0,1.20d0,1.40d0,1.82d0,0.00d0,0.00d0,1.70d0,1.55d0,1.50d0,1.47d0,1.54d0,2.27d0, &
     1.73d0,2.30d0,2.10d0,1.80d0,1.80d0,1.75d0,1.88d0,2.75d0,0.00d0,0.00d0,0.00d0,1.22d0,0.00d0,0.00d0,0.00d0,0.00d0,1.63d0, &
     1.40d0,1.39d0,1.87d0,0.00d0,1.85d0,1.90d0,1.85d0,2.02d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.0d00,0.00d0, &
     1.63d0,1.72d0,1.58d0,1.93d0,2.17d0,0.00d0,2.06d0,1.98d0,2.16d0/

! Criteria ?

iPL = iPL_espf()-1
RMax = dble(IRMax)
if (Process .and. (iPL >= 3)) then
  write(IOut,'(A,I2,A)') ' Max : ',IRMax,' van der Waals radii'
  write(IOut,'(A,F4.2,A)') ' ... with ',DeltaR,' angstroms between grid points.'
end if

! The VDW scaling factor is currently hard-coded. Too be changed.

VDWFact = 1.0d0

! SEARCH FOR THE EXTREMA OF THE MOLECULAR GEOMETRY WITHIN A CUBE.

xmax = -Huge
xmin = Huge
ymax = -Huge
ymin = Huge
zmax = -Huge
zmin = Huge

RVDWMX = Zero
do iAtom=1,nAtoms
  if (IsMM(iAtom) /= 0) then
    if (Process .and. (iPL >= 3)) then
      write(6,11) iAtom
    end if
    cycle
  end if
  xmax = max(xmax,CO(1,iAtom))
  xmin = min(xmin,CO(1,iAtom))
  ymax = max(ymax,CO(2,iAtom))
  ymin = min(ymin,CO(2,iAtom))
  zmax = max(zmax,CO(3,iAtom))
  zmin = min(zmin,CO(3,iAtom))
  RVDWMX = max(RVDWMX,RVDW(IAn(iAtom)))
end do
rSpher = RMax*RVDWMX

! DETERMINE THE MINIMUM PARALLELEPIPED DIMENSIONS REQUIRED TO
! CONTAIN THE MOLECULE, INCLUDING A MAXIMUM SELECTION RADIUS
! ON BOTH SIDES.

xrange = xmax-xmin+Two*rSpher
yrange = ymax-ymin+Two*rSpher
zrange = zmax-zmin+Two*rSpher

nbptx = int(xrange/DeltaR)
nbpty = int(yrange/DeltaR)
nbptz = int(zrange/DeltaR)

! DETERMINE THE MAXIMUM NUMBER OF POINTS FOR THE FITTING PROCEDURE.

maxpt = nbptx*nbpty*nbptz

! PRINT OUT INTERMEDIATE RESULTS, IF SO REQUESTED

if (iPL > 2) then
  write(IOut,'(" EXTREMA OF THE MOLECULAR GEOMETRY : ")')
  write(IOut,'(a)')
  write(IOut,'(" X_min = ",f8.4,4x,"X_max = ",f8.4)') xmin,xmax
  write(IOut,'(" Y_min = ",f8.4,4x,"Y_max = ",f8.4)') ymin,ymax
  write(IOut,'(" Z_min = ",f8.4,4x,"Z_max = ",f8.4)') zmin,zmax
  write(IOut,'(a)')
  write(IOut,'(" RE-SCALED PARALLELEPIPED SIZE : ")')
  write(IOut,'(a)')
  write(IOut,'(" X = ",f8.4)') xrange
  write(IOut,'(" Y = ",f8.4)') yrange
  write(IOut,'(" Z = ",f8.4)') zrange
  write(IOut,'(a)')
  write(IOut,'(" NUMBER OF POINTS PER DIRECTION : ")')
  write(IOut,'(a)')
  write(IOut,'(" Nb_x = ",i5)') nbptx
  write(IOut,'(" Nb_y = ",i5)') nbpty
  write(IOut,'(" Nb_z = ",i5)') nbptz
  write(IOut,'(a)')
  write(IOut,'(" MAXIMUM POSSIBLE NUMBER OF POINTS : ")')
  write(IOut,'(a)')
  write(IOut,'(" Nb_max = ",i10)') maxpt
  write(IOut,'(a)')
end if

! LOOP OVER ALL POSSIBLE POINTS IN THE PARALLELEPIPED.

nPts = 0

do Ipx=0,nbptx
  pointa = xmin-rSpher+dble(Ipx)*DeltaR

  do Ipy=0,nbpty
    pointb = ymin-rSpher+dble(Ipy)*DeltaR

    do Ipz=0,nbptz
      pointc = zmin-rSpher+dble(Ipz)*DeltaR

      ! IS THIS POINT WITHIN THE VAN DER WAALS ENVELOPE ?
      !      Retain1 => false
      ! IS THIS POINT BEYOND THE EXCLUSION RADIUS FROM ALL ATOMS ?
      !      Retain2 => false

      Retain1 = .true.
      Retain2 = .false.

      iAtom = 0
      do while (Retain1)
        iAtom = iAtom+1
        if (iAtom > nAtoms) exit
        if (IsMM(iAtom) /= 0) cycle
        RZ = RVDW(IAn(iAtom))
        VDWEnv = RZ*VDWFact
        Distance = (pointa-CO(1,iAtom))**2+(pointb-CO(2,iAtom))**2+(pointc-CO(3,iAtom))**2
        Distance = sqrt(Distance)
        if (Retain1) Retain1 = Distance > VDWEnv
        Retain2 = Retain2 .or. (Distance <= rmax*RZ)
      end do

      ! STORE THE POINTS REQUIRED FOR THE FITTING PROCEDURE.

      if (Retain1 .and. Retain2) then
        nPts = nPts+1
        if (Process) then
          Grid(1,nPts) = pointa
          Grid(2,nPts) = pointb
          Grid(3,nPts) = pointc
        end if
      end if
    end do
  end do
end do

return

11 format(' MM atom',I3,' is ignored in the grid construction')

end subroutine PNT
