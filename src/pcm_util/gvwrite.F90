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

implicit none
integer(kind=iwp), intent(in) :: idx, nTs, NEsfP, nAt, IAt(nAt), NVert(nTs), ISphe(nTs), MxVert
real(kind=wp), intent(in) :: AtmC(3,nAt), Coor_Sph(4,NEsfP), Tessera(4,nTs), Vert(3,MxVert,nTs), Q(nTs)
integer(kind=iwp), intent(out) :: ivts(MxVert,nTs)
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
  if (idx == 2) call Colchg(qt,qmax,qmin,C1,C2,C3)
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
