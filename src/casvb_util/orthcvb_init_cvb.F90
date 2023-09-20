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

subroutine orthcvb_init_cvb()

use casvb_global, only: cvbnrm_fr, nfrag, nvb_fr

implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "WrkSpc.fh"

if (nfrag <= 1) then
  cvbnrm = ddot_(nvb,work(lv(2)),1,work(lv(2)),1)
else
  ifr_off = 0
  do ifrag=1,nfrag
    cvbnrm_fr(ifrag) = ddot_(nvb_fr(ifrag),work(ifr_off+lv(2)),1,work(ifr_off+lv(2)),1)
    ifr_off = ifr_off+nvb_fr(ifrag)
  end do
end if

return

end subroutine orthcvb_init_cvb
