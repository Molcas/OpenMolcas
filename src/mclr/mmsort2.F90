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

subroutine MMSort2(A,B,P,iel)

use Index_Functions, only: iTri, nTri_Elem
use MCLR_Data, only: DspVec, lDisp
use input_mclr, only: nSym, nTPert
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: A(*)
real(kind=wp), intent(_OUT_) :: B(*), P(*)
integer(kind=iwp), intent(out) :: iel(3)
integer(kind=iwp) :: iDisp, iG, iii, ijD, ijD1, ijG, ijP, iSym, jDisp
logical(kind=iwp) :: geomi, geomj

ijD = 0
iG = 0
ijG = 0
ijP = 0
iii = 0
iel(:) = 0
do iSym=1,nsym
  do idisp=1,ldisp(isym)
    geomi = btest(ntpert(idisp+iii),4)
    if (.not. geomi) then
      iG = iG+1
      iel(ig) = isym
      do jdisp=1,ldisp(isym)
        geomj = btest(ntpert(jdisp+iii),4)
        if (geomj) then
          ijg = ijg+1
          ijd1 = ijD+iTri(idisp,jdisp)
          B(ijG) = A(ijD1)
        else if (idisp <= jdisp) then
          ijP = iTri(dspvec(jdisp+iii),dspvec(idisp+iii))
          P(ijP) = A(ijD+iTri(idisp,jdisp))
        end if
      end do
    end if
  end do
  ijD = ijD+nTri_Elem(ldisp(isym))
  iii = iii+ldisp(isym)
end do

end subroutine MMSort2
