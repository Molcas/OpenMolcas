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

subroutine SGPrint(SGS)

use definitions, only: iwp, u6
use gugx, only: SGStruct

implicit none
type(SGStruct), intent(in) :: SGS
integer(kind=iwp) nLev, nVert, MidLev, MVSta, MVEnd, i, iv, ic

! Unpack structure SGS:
nLev = SGS%nLev
nVert = SGS%nVert
MidLev = SGS%MidLev
MVSta = SGS%MVSta
MVEnd = SGS%MVEnd

write(u6,*) ' Split-Graph UGA. Graph description:'
write(u6,*) ' Nr of levels:',nLev
write(u6,*) ' Orbital symmetry labels:'
write(u6,'(1x,30i2)') (SGS%ISm(i),i=1,nLev)
write(u6,*) ' Nr of vertices:',nVert
write(u6,*)
write(u6,*) ' Vertex    L  N    A  B  C      Downchain table        Upchain table'
write(u6,*)
do iv=1,nVert
  write(u6,'(1x,i4,5x,2i3,2x,3i3,5x,4i4,5x,4i4)') iv,(SGS%DRT(iv,i-1),i=1,5),(SGS%Down(iv,ic),ic=0,3),(SGS%Up(iv,ic),ic=0,3)
end do
write(u6,*)
write(u6,*) ' Mid Level:',MidLev
write(u6,*) ' Mid Vertices:',MVSta,'...',MVEnd
write(u6,*)
write(u6,*) ' Modified Arc Weight table:'
write(u6,*) '           Coupling case number'
write(u6,*) ' Vertex      0    1    2    3'
write(u6,*)
do iv=1,nVert
  write(u6,'(1x,i4,5x,4i5)') iv,(SGS%MAW(iv,ic),ic=0,3)
end do

end subroutine SGPrint
