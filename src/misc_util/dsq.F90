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

subroutine DSQ(A,B,ICB,IRB,NROW)

use Constants, only: Two, Half
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: A(*)
real(kind=wp), intent(_OUT_) :: B(*)
integer(kind=iwp), intent(in) :: ICB, IRB, NROW
integer(kind=iwp) :: ICOL, IND, IROW

IND = 0
do IROW=0,NROW-1
  do ICOL=0,IROW
    IND = IND+1
    B(1+IROW*ICB+ICOL*IRB) = Half*A(IND)
    B(1+ICOL*ICB+IROW*IRB) = Half*A(IND)
  end do
  B(1+IROW*(ICB+IRB)) = Two*B(1+IROW*(ICB+IRB))
end do

return

end subroutine DSQ
