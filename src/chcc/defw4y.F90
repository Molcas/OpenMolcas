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

subroutine DefW4y(aGrp,bGrp,cGrp,dGrp,w4y)
! this routine does:
! define w4y key, which indicates, if at least one
! W4 file needs to be calculated on this node for
! given a',b',c',d'

#ifdef _MOLCAS_MPP_
use Index_Functions, only: nTri_Elem
use chcc_global, only: GrpaLow, GrpaUp, InqW4
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: aGrp, bGrp, cGrp, dGrp
integer(kind=iwp), intent(out) :: w4y
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: abSGrp, aSGrp, bSGrp, bSGrpUp, cdSGrp, cSGrp, dSGrp, dSGrpUp

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

    ! cycle over  c">=d" for c',d'
    do cSGrp=GrpaLow(cGrp),GrpaUp(cGrp)
      if (cGrp == dGrp) then
        dSGrpUp = cSGrp
      else
        dSGrpUp = GrpaUp(dGrp)
      end if
      do dSGrp=GrpaLow(dGrp),dSGrpUp
        cdSGrp = nTri_Elem(cSGrp-1)+dSGrp
        if (InqW4(abSGrp,cdSGrp)) w4y = w4y+1
      end do
    end do

    ! end cycle over a">=b" subgroups
  end do
end do

#else

#include "macros.fh"
unused_var(aGrp)
unused_var(bGrp)
unused_var(cGrp)
unused_var(dGrp)

w4y = 1
#endif

return

end subroutine DefW4y
