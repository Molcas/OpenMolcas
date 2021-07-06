!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine GVWrite(idx,nTs,NEsfP,nAt,AtmC,IAt,Coor_Sph,Tessera,NVert,Vert,ISphe,Q,ivts,MxVert)

use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: idx, nTs, NEsfP, nAt, IAt(nAt), NVert(*), ISphe(nTs), MxVert
real(kind=wp), intent(in) :: AtmC(3,nAt), Coor_Sph(4,*), Tessera(4,nTs), Vert(3,MxVert,*), Q(*)
integer(kind=iwp), intent(_OUT_) :: ivts(MxVert,*)
integer(kind=iwp) :: i, j, jcord, K, Last, Lu, N, NumV
real(kind=wp) :: C1, C2, C3, qmax, qmin, qt
integer(kind=iwp), external :: IsFreeUnit

! Prepare the input file for GeomView (coloured polyhedra)

Lu = IsFreeUnit(21)
if (idx == 1) call molcas_open(Lu,'GV.off')
if (idx == 2) call molcas_open(Lu,'GV1.off')
NumV = 0
do i=1,NTs
  NumV = NumV+NVert(i)
end do
write(Lu,500) 'COFF'
write(Lu,1000) NumV,NTs,NumV
K = 0
Last = 0
qt = Zero
if (idx == 2) then
  qmax = Zero
  qmin = Zero
  do i=1,NTs
    qt = q(i)/Tessera(4,i)
    if (qt >= qmax) qmax = qt
    if (qt <= qmin) qmin = qt
  end do
  write(Lu,1100) QMin,QMax
end if
do i=1,NTs
  if (idx == 2) qt = q(i)/Tessera(4,i)
  N = ISphe(i)
  if (N /= Last) write(Lu,1500) N
  Last = N
  if (idx == 1) call Colour(NEsfP,NAt,AtmC,IAt,Coor_Sph,N,C1,C2,C3)
  if (idx == 2) call Colchg(i,qt,qmax,qmin,C1,C2,C3)
  do j=1,NVert(i)
    IVTS(j,i) = K
    K = K+1
    write(Lu,2001) (VERT(jcord,j,i),jcord=1,3),C1,C2,C3,0.75_wp,i
  end do
end do
do i=1,NTs
  write(Lu,3000) NVert(i),(IVTS(j,i),j=1,nvert(i))
end do
close(Lu)

return

500 format(1x,a)
1000 format(3i10)
1100 format('# Minimum and maximum charge density ',2f12.6)
1500 format('# Sphere number ',i4)
2001 format('  ',3f16.9,4f5.2,' # Tess. ',i4)
3000 format('  ',14i10)

end subroutine GVWrite
!====
subroutine Colour(NesfP,NAt,AtmC,IAt,Coor_Sph,N,C1,C2,C3)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NesfP, NAt, IAt(*), N
real(kind=wp), intent(in) :: AtmC(3,*), Coor_Sph(4,*)
real(kind=wp), intent(out) :: C1, C2, C3
integer(kind=iwp) :: I, J
real(kind=wp) :: Diff
character(len=20) :: Col
real(kind=wp), parameter :: Delta = 1.0e-3_wp

! Assign tesserae colours for GeomView:
! Carbon: green, nitrogen: blue, oxygen: red, hydrogen: light blue,
! others: violet, added spheres: gray.

Col = ' '
if (N > NESFP) then
  Col = 'Gray'
  call ColTss(u6,Col,C1,C2,C3)
  return
end if
I = 0
do J=1,NAt
  I = I+1
  Diff = sqrt((AtmC(1,J)-Coor_Sph(1,N))**2+(AtmC(2,J)-Coor_Sph(2,N))**2+(AtmC(3,J)-Coor_Sph(3,N))**2)
  if (Diff < Delta) then
    if (IAt(I) == 6) then
      Col = 'Green'
    elseif (IAt(I) == 7) then
      Col = 'Blue'
    elseif (IAt(I) == 8) then
      Col = 'Red'
    elseif (IAt(I) == 1) then
      Col = 'Light Blue'
    else
      Col = 'Fuchsia'
    end if
  end if
end do
call ColTss(u6,Col,C1,C2,C3)

return

end subroutine Colour
!====
subroutine ColTss(IOut,Colour,C1,C2,C3)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IOut
character(len=*), intent(in) :: Colour
real(kind=wp), intent(out) :: C1, C2, C3

! Assign tesserae colours for GeomView:

if (Colour == 'White') then
  C1 = 1.0_wp
  C2 = 1.0_wp
  C3 = 1.0_wp
else if (Colour == 'Gray') then
  C1 = 0.66_wp
  C2 = 0.66_wp
  C3 = 0.66_wp
else if ((Colour == 'Blue') .or. (Colour == 'Dark Blue')) then
  C1 = 0.0_wp
  C2 = 0.0_wp
  C3 = 1.0_wp
else if (Colour == 'Light Blue') then
  C1 = 0.0_wp
  C2 = 1.0_wp
  C3 = 1.0_wp
else if (Colour == 'Green') then
  C1 = 0.0_wp
  C2 = 1.0_wp
  C3 = 0.0_wp
else if (Colour == 'Yellow') then
  C1 = 1.0_wp
  C2 = 1.0_wp
  C3 = 0.0_wp
else if (Colour == 'Orange') then
  C1 = 1.0_wp
  C2 = 0.5_wp
  C3 = 0.0_wp
else if (Colour == 'Violet') then
  C1 = 0.6_wp
  C2 = 0.0_wp
  C3 = 1.0_wp
else if ((Colour == 'Pink') .or. (Colour == 'Light Red')) then
  C1 = 1.0_wp
  C2 = 0.5_wp
  C3 = 1.0_wp
else if (Colour == 'Fuchsia') then
  C1 = 1.0_wp
  C2 = 0.0_wp
  C3 = 1.0_wp
else if ((Colour == 'Red') .or. (Colour == 'Dark Red')) then
  C1 = 1.0_wp
  C2 = 0.0_wp
  C3 = 0.0_wp
else if (Colour == 'Black') then
  C1 = 0.0_wp
  C2 = 0.0_wp
  C3 = 0.0_wp
else
  C1 = 0.0_wp  ! dummy assignement
  C2 = 0.0_wp  ! dummy assignement
  C3 = 0.0_wp  ! dummy assignement
  write(IOut,'(a)') 'Unrecognized colour in ColTss'
  call Abend()
end if

return

end subroutine ColTss
!====
subroutine ColChg(i,q,QMAX,QMIN,C1,C2,C3)

use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: i
real(kind=wp), intent(in) :: q, QMAX, QMIN
real(kind=wp), intent(out) :: C1, C2, C3
real(kind=wp) :: DNeg, DPos
character(len=20) :: Colour

! Assign tesserae colours for GeomView:

DPos = QMax*Half
DNeg = QMin*Half
! Total charge density < 0
if (q < DNeg) then
  Colour = 'Dark Blue'
else if (q < Zero) then
  Colour = 'Light Blue'
! Total charge density > 0
else if (q < DPos) then
  Colour = 'Pink'
else
  Colour = 'Red'
end if
call ColTss(u6,Colour,C1,C2,C3)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(i)

end subroutine ColChg
