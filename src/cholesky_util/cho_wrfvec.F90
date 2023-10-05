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

subroutine CHO_WRFVEC(VEC,ISYMA,ISYMB,IVEC1,NUMV)
!
! Purpose: write full storage vectors to disk.

use Cholesky, only: LUFV, NABPK
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_IN_) :: VEC(*)
integer(kind=iwp), intent(in) :: ISYMA, ISYMB, IVEC1, NUMV
integer(kind=iwp) :: IADR, IOPT, NTOT

IOPT = 1
IADR = NABPK(ISYMA,ISYMB)*(IVEC1-1)+1
NTOT = NABPK(ISYMA,ISYMB)*NUMV
call DDAFILE(LUFV(ISYMA,ISYMB),IOPT,VEC,NTOT,IADR)

end subroutine CHO_WRFVEC
