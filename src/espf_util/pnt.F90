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

use Constants, only: Zero, One, Two, Angstrom
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IOut, nAtoms, IRMax, IAn(nAtoms), IsMM(nAtoms)
real(kind=wp), intent(in) :: CO(3,nAtoms), DeltaR
integer(kind=iwp), intent(out) :: nPts
real(kind=wp), intent(inout) :: Grid(3,*)
logical(kind=iwp), intent(in) :: Process
real(kind=wp) :: Distance, pointa, pointb, pointc, RMax, rSpher, RVDWMX, RZ, VDWEnv, VDWFact, xmax, xmin, xrange, ymax, ymin, &
                 yrange, zmax, zmin, zrange
integer(kind=iwp) :: iAtom, iPL, Ipx, Ipy, Ipz, maxpt, nbptx, nbpty, nbptz
logical(kind=iwp) :: Retain1, Retain2
real(kind=wp), parameter :: rHuge = 1.0e8_wp, &
                            ! These are mostly Bondi radii (in angstrom), but there are some differences
                            RVDW(0:54) = [2.00_wp, &                                                                    ! X
                                          1.20_wp,1.40_wp, &                                                            ! H-He
                                          1.82_wp,Zero,Zero,1.70_wp,1.55_wp,1.50_wp,1.47_wp,1.54_wp, &                  ! Li-Ne
                                          2.27_wp,1.73_wp,2.30_wp,2.10_wp,1.80_wp,1.80_wp,1.75_wp,1.88_wp, &            ! Na-Ar
                                          2.75_wp,Zero,Zero,Zero,1.22_wp,Zero,Zero,Zero,Zero,1.63_wp,1.40_wp,1.39_wp, & ! K-Zn
                                            1.87_wp,Zero,1.85_wp,1.90_wp,1.85_wp,2.02_wp, &                             ! Ga-Kr
                                          Zero,Zero,Zero,Zero,Zero,Zero,Zero,Zero,Zero,1.63_wp,1.72_wp,1.58_wp, &       ! Rb-Cd
                                            1.93_wp,2.17_wp,Zero,2.06_wp,1.98_wp,2.16_wp]                               ! In-Xe
integer(kind=iwp), external :: iPL_espf

! Criteria ?

iPL = iPL_espf()-1
RMax = real(IRMax,kind=wp)
if (Process .and. (iPL >= 3)) then
  write(IOut,'(A,I2,A)') ' Max : ',IRMax,' van der Waals radii'
  write(IOut,'(A,F4.2,A)') ' ... with ',DeltaR*Angstrom,' angstroms between grid points.'
end if

! The VDW scaling factor is currently hard-coded. Too be changed.

VDWFact = One

! SEARCH FOR THE EXTREMA OF THE MOLECULAR GEOMETRY WITHIN A CUBE.

xmax = -rHuge
xmin = rHuge
ymax = -rHuge
ymin = rHuge
zmax = -rHuge
zmin = rHuge

RVDWMX = Zero
do iAtom=1,nAtoms
  if (IsMM(iAtom) /= 0) then
    if (Process .and. (iPL >= 3)) then
      write(u6,11) iAtom
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
rSpher = RMax*RVDWMX/Angstrom

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
  write(IOut,'(" X_min = ",f8.4,4x,"X_max = ",f8.4)') xmin*Angstrom,xmax*Angstrom
  write(IOut,'(" Y_min = ",f8.4,4x,"Y_max = ",f8.4)') ymin*Angstrom,ymax*Angstrom
  write(IOut,'(" Z_min = ",f8.4,4x,"Z_max = ",f8.4)') zmin*Angstrom,zmax*Angstrom
  write(IOut,'(a)')
  write(IOut,'(" RE-SCALED PARALLELEPIPED SIZE : ")')
  write(IOut,'(a)')
  write(IOut,'(" X = ",f8.4)') xrange*Angstrom
  write(IOut,'(" Y = ",f8.4)') yrange*Angstrom
  write(IOut,'(" Z = ",f8.4)') zrange*Angstrom
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
  pointa = xmin-rSpher+real(Ipx,kind=wp)*DeltaR

  do Ipy=0,nbpty
    pointb = ymin-rSpher+real(Ipy,kind=wp)*DeltaR

    do Ipz=0,nbptz
      pointc = zmin-rSpher+real(Ipz,kind=wp)*DeltaR

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
        RZ = RVDW(IAn(iAtom))/Angstrom
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
