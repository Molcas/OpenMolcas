!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine IFOCK(FC,NI,NJ,NK,FINI,II)

use Constants, only: Two
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: FC(*)
integer(kind=iwp), intent(in) :: NI, NJ, NK, II
real(kind=wp), intent(in) :: FINI
integer(kind=iwp) :: JKPOS

if (NI > 0) return
if ((NJ <= 0) .or. (NK <= 0)) return
JKPOS = NJ*(NJ-1)/2+NK
if (NK > NJ) JKPOS = NK*(NK-1)/2+NJ
if (II /= 0) then
  FC(JKPOS) = FC(JKPOS)+Two*FINI
else
  FC(JKPOS) = FC(JKPOS)-FINI
end if

return

end subroutine IFOCK
