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

integer function nMo(mIrr)

use Symmetry_Info, only: nIrrep

implicit integer(a-h,o-z)
#include "etwas.fh"

nInt = 0
nA = 0
do iS=0,nIrrep-1
  nA = nA+nAsh(is)
end do
NMM = nA*(nA+1)/2
nInt = nMM*(nMM+1)/2
nMo = nInt

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(mIrr)

end function nMo
