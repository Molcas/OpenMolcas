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

subroutine GVWrite(Index,NTs,NEsfP,NAt,AtmC,IAt,Coor_Sph,Tessera,NVert,Vert,ISphe,Q,ivts,MxVert)

implicit real*8(A-H,O-Z)
dimension NVert(*), Vert(3,MxVert,*), IVTS(MxVert,*), Tessera(4,nTs), Q(*)
dimension AtmC(3,nAt), IAt(nAt), Coor_Sph(4,*), ISphe(nTs)

! Prepare the input file for GeomView (coloured polyhedra)

Lu = 21
Lu = IsFreeUnit(Lu)
if (Index == 1) call molcas_open(Lu,'GV.off')
if (Index == 2) call molcas_open(Lu,'GV1.off')
NumV = 0
do i=1,NTs
  NumV = NumV+NVert(i)
end do
write(Lu,500) 'COFF'
write(Lu,1000) NumV,NTs,NumV
K = 0
Last = 0
qt = 0.0d0
if (Index == 2) then
  qmax = 0.0d0
  qmin = 0.0d0
  do i=1,NTs
    qt = q(i)/Tessera(4,i)
    if (qt >= qmax) qmax = qt
    if (qt <= qmin) qmin = qt
  end do
  write(Lu,1100) QMin,QMax
end if
do i=1,NTs
  if (Index == 2) qt = q(i)/Tessera(4,i)
  N = ISphe(i)
  if (N /= Last) write(Lu,1500) N
  Last = N
  if (Index == 1) call Colour(NEsfP,NAt,AtmC,IAt,Coor_Sph,N,C1,C2,C3)
  if (Index == 2) call Colchg(i,qt,qmax,qmin,C1,C2,C3)
  do j=1,NVert(i)
    IVTS(j,i) = K
    K = K+1
    write(Lu,2001) (VERT(jcord,j,i),jcord=1,3),C1,C2,C3,0.75,i
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

implicit real*8(A-H,O-Z)
character*20 Col
dimension IAt(*), AtmC(3,*), Coor_Sph(4,*)

! Assign tesserae colours for GeomView:
! Carbon: green, nitrogen: blue, oxygen: red, hydrogen: light blue,
! others: violet, added spheres: gray.

Col = ' '
IOut = 6
Delta = 1.d-03
if (N > NESFP) then
  Col = 'Gray'
  call ColTss(IOut,Col,C1,C2,C3)
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
call ColTss(IOut,Col,C1,C2,C3)

return

end subroutine Colour
!====
subroutine ColTss(IOut,Colour,C1,C2,C3)

implicit real*8(A-H,O-Z)
character*20 Colour

! Assign tesserae colours for GeomView:

if (Colour == 'White') then
  C1 = 1.0d0
  C2 = 1.0d0
  C3 = 1.0d0
else if (Colour == 'Gray') then
  C1 = 0.66d0
  C2 = 0.66d0
  C3 = 0.66d0
else if ((Colour == 'Blue') .or. (Colour == 'Dark Blue')) then
  C1 = 0.0d0
  C2 = 0.0d0
  C3 = 1.0d0
else if (Colour == 'Light Blue') then
  C1 = 0.0d0
  C2 = 1.0d0
  C3 = 1.0d0
else if (Colour == 'Green') then
  C1 = 0.0d0
  C2 = 1.0d0
  C3 = 0.0d0
else if (Colour == 'Yellow') then
  C1 = 1.0d0
  C2 = 1.0d0
  C3 = 0.0d0
else if (Colour == 'Orange') then
  C1 = 1.0d0
  C2 = 0.5d0
  C3 = 0.0d0
else if (Colour == 'Violet') then
  C1 = 0.6d0
  C2 = 0.0d0
  C3 = 1.0d0
else if ((Colour == 'Pink') .or. (Colour == 'Light Red')) then
  C1 = 1.0d0
  C2 = 0.5d0
  C3 = 1.0d0
else if (Colour == 'Fuchsia') then
  C1 = 1.0d0
  C2 = 0.0d0
  C3 = 1.0d0
else if ((Colour == 'Red') .or. (Colour == 'Dark Red')) then
  C1 = 1.0d0
  C2 = 0.0d0
  C3 = 0.0d0
else if (Colour == 'Black') then
  C1 = 0.0d0
  C2 = 0.0d0
  C3 = 0.0d0
else
  C1 = 0.0d0  ! dummy assignement
  C2 = 0.0d0  ! dummy assignement
  C3 = 0.0d0  ! dummy assignement
  write(IOut,'(a)') 'Unrecognized colour in ColTss'
  call Abend()
end if

return

end subroutine ColTss
!====
subroutine ColChg(i,q,QMAX,QMIN,C1,C2,C3)

implicit real*8(A-H,O-Z)
character*20 Colour
save Colour
data Colour/'                    '/

! Assign tesserae colours for GeomView:

Zero = dble(0)
DPos = QMax/2.d0
DNeg = QMin/2.d0
! Total charge density < 0
if (q < DNeg) Colour = 'Dark Blue'
if ((q >= DNeg) .and. (q < Zero)) Colour = 'Light Blue'
! Total charge density > 0
if ((q >= Zero) .and. (q < DPos)) Colour = 'Pink'
if (q >= DPos) Colour = 'Red'
call ColTss(IOut,Colour,C1,C2,C3)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(i)

end subroutine ColChg
