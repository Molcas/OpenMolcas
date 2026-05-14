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

subroutine FRMDSCN(VEC,NREC,LBLK,LU)
! Read VEC as multiple record file, NREC records read

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: VEC(*)
integer(kind=iwp), intent(in) :: NREC, LBLK, LU
integer(kind=iwp) :: IAMPACK, IMZERO, IOFF, IREC, LREC(1)

IOFF = 1
do IREC=1,NREC
  call IFRMDS(LREC(1),1,LBLK,LU)
  call FRMDSC(VEC(IOFF),LREC(1),LBLK,LU,IMZERO,IAMPACK)
  IOFF = IOFF+LREC(1)
end do

end subroutine FRMDSCN
