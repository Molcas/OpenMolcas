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

subroutine MakeWw(Ww,W1,W2,aSGrp,bSGrp,beSGrp,gaSGrp,key)
! this routine does:
! Make Ww(+,-)((ab)",(bega)") from W1(a",be",b",ga")
!                              and W2(b",be",a",ga")
!  where index ab" is always a">b" (not a"=b")
!  and   index bega" is be">=ga" for W+ and be">ga" for W-
! for all cases aSGrp >= bSGrp
! N.B. for a"=b" W2 is not used (although it is passed through
!      the header (because of uniform shape) but no matter what
!      is there)
! N.B. rutinky MakeWwHlpx niesu prilis vymakane @
!
! parameter description:
! Ww    - array for Ww+(-) (O)
! W1    - array for W1(a",be",b",ga") (I)
! W2    - array for W1(b",be",a",ga") (I)
! xSGrp - SubGroup of a",b",be",ga" (I)
! key   - 1 - calc Ww+, 2 - calc Ww- (I)

use Index_Functions, only: nTri_Elem
use chcc_global, only: DimSGrpa, DimSGrpbe
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: Ww(*)
real(kind=wp), intent(in) :: W1(*), W2(*)
integer(kind=iwp), intent(in) :: aSGrp, bSGrp, beSGrp, gaSGrp, key
integer(kind=iwp) :: dima, dimab, dimb, dimbe, dimbega, dimga

!1 def dimensions

dimbe = DimSGrpbe(beSGrp)
dimga = DimSGrpbe(gaSGrp)
dima = DimSGrpa(aSGrp)
dimb = DimSGrpa(bSGrp)

if (beSGrp == gaSGrp) then
  if (key == 1) then
    dimbega = nTri_Elem(dimbe)
  else
    dimbega = nTri_Elem(dimbe-1)
  end if
else
  dimbega = dimbe*dimga
end if

if (aSGrp == bSGrp) then
  dimab = nTri_Elem(dima-1)
else
  dimab = dima*dimb
end if

!2 Make Ww matrix

if (aSGrp == bSGrp) then
  if (beSGrp == gaSGrp) then
    call MakeWwHlp1(Ww,W1,dima,dimb,dimab,dimbe,dimga,dimbega,key)
  else
    call MakeWwHlp2(Ww,W1,dima,dimb,dimab,dimbe,dimga,key)
  end if
else
  if (beSGrp == gaSGrp) then
    call MakeWwHlp3(Ww,W1,W2,dima,dimb,dimbe,dimga,dimbega,key)
  else
    call MakeWwHlp4(Ww,W1,W2,dima,dimb,dimbe,dimga,key)
  end if
end if

return

end subroutine MakeWw
