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

subroutine SG_Print(iState)

use Definitions, only: iwp, u6
use sguga_states, only: SGS

implicit none
integer(kind=iwp), intent(in) :: iState
integer(kind=iwp) :: i, ic, iv

write(u6,*) ' Split-Graph UGA. Graph description:'
write(u6,*) ' Nr of levels:',SGS%nLev
write(u6,*) ' Orbital symmetry labels:'
write(u6,'(1x,30i2)') (SGS(istate)%ISm(i),i=1,SGS(istate)%nLev)
write(u6,*) ' Nr of vertices:',SGS(istate)%nVert
write(u6,*)
write(u6,*) ' Vertex    L  N    A  B  C      Downchain table        Upchain table'
write(u6,*)
do iv=1,SGS(istate)%nVert
  write(u6,'(1x,i4,5x,2i3,2x,3i3,5x,4i4,5x,4i4)') iv,(SGS(istate)%DRT(iv,i-1),i=1,5),(SGS(istate)%Down(iv,ic),ic=0,3), &
       (SGS(istate)%Up(iv,ic),ic=0,3)
end do
write(u6,*)
write(u6,*) ' Mid Level:',SGS(istate)%MidLev
write(u6,*) ' Mid Vertices:',SGS(istate)%MVSta,'...',SGS(istate)%MVEnd
write(u6,*)
write(u6,*) ' Modified Arc Weight table:'
write(u6,*) '           Coupling case number'
write(u6,*) ' Vertex      0    1    2    3'
write(u6,*)
do iv=1,SGS(istate)%nVert
  write(u6,'(1x,i4,5x,4i5)') iv,(SGS(istate)%MAW(iv,ic),ic=0,3)
end do

end subroutine SG_Print
