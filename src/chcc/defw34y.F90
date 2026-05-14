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

subroutine DefW34y(aGrp,bGrp,w3y,w4y,NaGrp)
! this routine does:
! define w3y and w4y keys, which indicates, if at least one
! W3/W4 file needs to be calculated on this node for given a',b'

#ifdef _MOLCAS_MPP_
use Index_Functions, only: nTri_Elem
use chcc_global, only: GrpaLow, GrpaUp, InqW3, InqW4
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: aGrp, bGrp, NaGrp
integer(kind=iwp), intent(out) :: w3y, w4y
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: abSGrp, aSGrp, bSGrp, bSGrpUp, cdSGrp, cSGrp, NSGrp

w3y = 0
w4y = 0

! cycle over a">=b" subgroups for a',b'
do aSGrp=GrpaLow(aGrp),GrpaUp(aGrp)
  if (aGrp == bGrp) then
    bSGrpUp = aSGrp
  else
    bSGrpUp = GrpaUp(bGrp)
  end if
  do bSGrp=GrpaLow(bGrp),bSGrpUp
    abSGrp = nTri_Elem(aSGrp-1)+bSGrp

    ! cycle over all c"
    NSGrp = GrpaUp(NaGrp)
    do cSGrp=1,NSGrp
      if (InqW3(abSGrp,cSGrp)) w3y = w3y+1
    end do

    ! cycle over all c">=d"
    do cdSGrp=1,nTri_Elem(NSGrp)
      if (InqW4(abSGrp,cdSGrp)) w4y = w4y+1
    end do

    ! end cycle over a">=b" subgroups
  end do
end do

#else

#include "macros.fh"
unused_var(aGrp)
unused_var(bGrp)
unused_var(NaGrp)

w3y = 1
w4y = 1
#endif

return

end subroutine DefW34y
