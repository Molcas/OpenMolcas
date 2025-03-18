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

subroutine MMSort(A,B,ldisp1)

use MCLR_Data, only: lDisp
use input_mclr, only: nSYm, nTPert

implicit none
real*8 A(*), B(*)
integer ldisp1(nsym)
logical geomi, geomj
integer ijD, iG, ijG, iii, iSym, iDisp, jDisp, ijD1, jG, ijG1
! Statement function
integer i, j, itri
itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

ijG = 0
ijD = 0
iii = 0
call icopy(nsym,[0],0,ldisp1,1)
do iSym=1,nsym
  iG = 0
  do idisp=1,ldisp(isym)
    geomi = iand(ntpert(iii+idisp),16) == 16
    if (geomi) then
      ldisp1(isym) = ldisp1(isym)+1
      iG = iG+1
      jG = 0
      do jdisp=1,idisp
        geomj = iand(ntpert(jdisp+iii),16) == 16
        if (geomj) then
          jG = jG+1
          ijg1 = ijG+itri(ig,jg)
          ijd1 = ijD+itri(idisp,jdisp)
          B(ijG1) = A(ijD1)
        end if
      end do
    end if
  end do
  ijG = ijG+iG*(iG+1)/2
  ijD = ijD+ldisp(isym)*(ldisp(isym)+1)/2
  iii = iii+ldisp(isym)
end do

end subroutine MMSort
