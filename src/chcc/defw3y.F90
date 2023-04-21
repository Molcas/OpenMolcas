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

subroutine DefW3y(aGrp,bGrp,cGrp,w3y)
! this routine does:
! define w3y key, which indicates, if at least one
! W3 file needs to be calculated on this node for given
! a',b',c'

#ifdef _MOLCAS_MPP_
use Index_Functions, only: nTri_Elem
use chcc_global, only: GrpaLow, GrpaUp, InqW3
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: aGrp, bGrp, cGrp
integer(kind=iwp), intent(out) :: w3y
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: abSGrp, aSGrp, bSGrp, bSGrpUp, cSGrp

w3y = 0

! cycle over a">=b" subgroups for a',b'
do aSGrp=GrpaLow(aGrp),GrpaUp(aGrp)
  if (aGrp == bGrp) then
    bSGrpUp = aSGrp
  else
    bSGrpUp = GrpaUp(bGrp)
  end if
  do bSGrp=GrpaLow(bGrp),bSGrpUp
    abSGrp = nTri_Elem(aSGrp-1)+bSGrp

    ! cycle over c" for c'
    do cSGrp=GrpaLow(cGrp),GrpaUp(cGrp)
      if (InqW3(abSGrp,cSGrp)) w3y = w3y+1
    end do

    ! end cycle over a">=b" subgroups
  end do
end do

#else

#include "macros.fh"
unused_var(aGrp)
unused_var(bGrp)
unused_var(cGrp)

w3y = 1
#endif

return

end subroutine DefW3y
