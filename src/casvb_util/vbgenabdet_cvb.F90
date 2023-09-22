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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine vbgenabdet_cvb(idetavb,idetbvb,iconfs,nconf,nconfion,ndetvb,nel,noe,nalf,nbet,norb)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ndetvb, idetavb(ndetvb), idetbvb(ndetvb), noe, nconf, iconfs(noe,nconf), nel, nconfion(0:nel), nalf, nbet, norb
#include "WrkSpc.fh"
integer(kind=iwp) :: i1, i2, i3, i4, i5, i6
integer(kind=iwp), external :: mstacki_cvb

i1 = mstacki_cvb((norb+1)*(nalf+1))
i2 = mstacki_cvb((norb+1)*(nbet+1))
i3 = mstacki_cvb(norb+1)
i4 = mstacki_cvb(norb+1)
i5 = mstacki_cvb(norb)
i6 = mstacki_cvb(norb)
call vbgenabdet2_cvb(idetavb,idetbvb,iconfs,nconf,nconfion,ndetvb,nel,noe,nalf,nbet,norb,iwork(i1),iwork(i2),iwork(i3),iwork(i4), &
                     iwork(i5),iwork(i6))
call mfreei_cvb(i1)

return

end subroutine vbgenabdet_cvb
