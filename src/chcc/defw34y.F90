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

implicit none
integer aGrp, bGrp, w3y, w4y, NaGrp
#include "chcc1.fh"
#include "o2v4.fh"
#ifdef _MOLCAS_MPP_
#include "parcc.fh"
! help variables
integer aSGrp, bSGrp, abSGrp, bSGrpUp, cSGrp, cdSGrp
integer NSGrp

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
    abSGrp = aSGrp*(aSGrp-1)/2+bSGrp

    ! cycle over all c"
    NSGrp = GrpaUp(NaGrp)
    do cSGrp=1,NSGrp
      if (InqW3(abSGrp,cSGrp)) w3y = w3y+1
    end do

    ! cycle over all c">=d"
    do cdSGrp=1,NSGrp*(NSGrp+1)/2
      if (InqW4(abSGrp,cdSGrp)) w4y = w4y+1
    end do

    ! end cycle over a">=b" subgroups
  end do
end do

#else
w3y = 1
w4y = 1
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(aGrp)
  call Unused_integer(bGrp)
  call Unused_integer(NaGrp)
end if
#endif

return

end subroutine DefW34y
