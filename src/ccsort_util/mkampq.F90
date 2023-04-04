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

subroutine mkampq(wrk,wrksize,a,ammap)
! this routine reconstructs #2 V2<_a,m|p,q> from corresponding TEMPDA2 file

use ccsort_global, only: lunda2, map2, mbas, NSYM, reclen
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize, a, ammap(mbas,8,8)
real(kind=wp), intent(_OUT_) :: wrk(wrksize)
integer(kind=iwp) :: iiv2, irec0, length, pos, symm, symp

! loops over symmetry combinations
do symm=1,nsym
  do symp=1,nsym

    ! def initial record position in TEMPDA2
    ! and corresponding position and length in wrk (#2)

    irec0 = ammap(a,symm,symp)
    iiv2 = map2%i(symm,symp,1)
    pos = map2%d(iiv2,1)
    length = map2%d(iiv2,2)

    if (length > 0) call daread(lunda2,irec0,wrk(pos),length,reclen)

  end do
end do

return

end subroutine mkampq
