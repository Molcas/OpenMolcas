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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************

subroutine Get_Int_DCCD(rc,Xint,ipq,Nrs)

use Index_Functions, only: nTri_Elem
use GetInt_mod, only: nBas
use TwoDat, only: rcTwo
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp), intent(in) :: ipq, Nrs
real(kind=wp), intent(_OUT_) :: Xint(Nrs)

integer(kind=iwp) :: Npq

Npq = nTri_Elem(nBas(1))

if ((ipq < 1) .and. (ipq <= Npq)) then
  rc = rcTwo%RD10
  write(u6,*) 'ipq out of bounds: ',ipq
  call Abend()
end if

call GEN_INT_DCCD(rc,ipq,Xint)

end subroutine Get_Int_DCCD
