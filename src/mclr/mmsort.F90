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

use Index_Functions, only: iTri, nTri_Elem
use MCLR_Data, only: lDisp
use input_mclr, only: nSym, nTPert
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: A(*)
real(kind=wp), intent(_OUT_) :: B(*)
integer(kind=iwp), intent(out) :: ldisp1(nsym)
integer(kind=iwp) :: iDisp, iG, iii, ijD, ijD1, ijG, ijG1, iSym, jDisp, jG
logical(kind=iwp) :: geomi, geomj

ijG = 0
ijD = 0
iii = 0
ldisp1(:) = 0
do iSym=1,nsym
  iG = 0
  do idisp=1,ldisp(isym)
    geomi = btest(ntpert(iii+idisp),4)
    if (geomi) then
      ldisp1(isym) = ldisp1(isym)+1
      iG = iG+1
      jG = 0
      do jdisp=1,idisp
        geomj = btest(ntpert(jdisp+iii),4)
        if (geomj) then
          jG = jG+1
          ijg1 = ijG+iTri(ig,jg)
          ijd1 = ijD+iTri(idisp,jdisp)
          B(ijG1) = A(ijD1)
        end if
      end do
    end if
  end do
  ijG = ijG+nTri_Elem(iG)
  ijD = ijD+nTri_Elem(ldisp(isym))
  iii = iii+ldisp(isym)
end do

end subroutine MMSort
