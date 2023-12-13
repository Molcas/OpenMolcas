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
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoRPA_MOTra_ReorderCMO(nSym,nBas,nOrb,CMOinp,CMOout)

#include "intent.fh"

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb(nSym)
real(kind=wp), intent(in) :: CMOinp(*)
real(kind=wp), intent(_OUT_) :: CMOout(*)
integer(kind=iwp) :: iSym, ip1, ip2, l

ip1 = 1
ip2 = 1
do iSym=1,nSym
  l = nBas(iSym)*nOrb(iSym)
  call dCopy_(l,CMOinp(ip1),1,CMOout(ip2),1)
  call fZero(CMOout(ip2+l),nBas(iSym)*nBas(iSym)-l)
  ip1 = ip1+l
  ip2 = ip2+nBas(iSym)*nBas(iSym)
end do

end subroutine ChoRPA_MOTra_ReorderCMO
