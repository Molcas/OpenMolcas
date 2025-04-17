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

subroutine ZSMCL(NSM,NOCTP,NSSO,ISTSM,ISTCL)
! set symmetry and class arrays for strings

use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NSM, NOCTP, NSSO(NOCTP,NSM)
integer(kind=iwp), intent(_OUT_) :: ISTSM(*), ISTCL(*)
integer(kind=iwp) :: ICL, IOFF, ISM

IOFF = 1
do ISM=1,NSM
  do ICL=1,NOCTP
    ISTSM(IOFF:IOFF+NSSO(ICL,ISM)-1) = ISM
    ISTCL(IOFF:IOFF+NSSO(ICL,ISM)-1) = ICL
    IOFF = IOFF+NSSO(ICL,ISM)
  end do
end do

return

end subroutine ZSMCL
