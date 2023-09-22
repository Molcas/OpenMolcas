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

subroutine psym1_cvb(civec1,civec2,osym,ientry)

use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: civec1(nda,ndb), civec2(nda,ndb), osym(mxirrep)
integer(kind=iwp) :: ientry
#include "WrkSpc.fh"
integer(kind=iwp) :: k1, k10, k11, k12, k13, k14, k15, k16, k17, k2, k3, k4, k5, k6, k7, k8, k9
integer(kind=iwp), external :: mstacki_cvb

k1 = mstacki_cvb(nda)
k2 = mstacki_cvb(ndb)
k3 = mstacki_cvb(mxirrep+1)
k4 = mstacki_cvb(mxirrep+1)
k5 = mstacki_cvb(nda)
k6 = mstacki_cvb(ndb)
k7 = mstacki_cvb(mxirrep)
k8 = mstacki_cvb(mxirrep)
k9 = mstacki_cvb(mxirrep)
k10 = mstacki_cvb(norb+1)
k11 = mstacki_cvb(norb+1)
k12 = mstacki_cvb(norb+1)
k13 = mstacki_cvb(norb+1)
k14 = mstacki_cvb(norb+1)
k15 = mstacki_cvb((norb+1)*(nalf+1))
k16 = mstacki_cvb((norb+1)*(nbet+1))
k17 = mstacki_cvb(mxirrep)

call symgen_cvb(nalf,nbet,nda,ndb,iwork(k1),iwork(k2),iwork(k3),iwork(k4),iwork(k5),iwork(k6),iwork(k7),iwork(k8),iwork(k9), &
                iwork(k10),iwork(k11),iwork(k12),iwork(k13),iwork(k14),iwork(k15),iwork(k16),iwork(k17))

call mfreei_cvb(k5)
call psym2_cvb(civec1,civec2,iwork(k1),iwork(k2),iwork(k3),iwork(k4),osym,ientry)
call mfreei_cvb(k1)

return

end subroutine psym1_cvb
