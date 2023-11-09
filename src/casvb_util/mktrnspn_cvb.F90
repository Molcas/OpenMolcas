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

subroutine mktrnspn_cvb()

use casvb_global, only: cvb, cvbdet, ipr, kbasis, kbasiscvb, nvb, spinb
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), external :: nvb_cvb

if (ipr(1) >= 1) write(u6,'(/,4a)') ' Changing spin basis : ',trim(spinb(kbasiscvb)),' --> ',trim(spinb(kbasis))
call str2vbc_cvb(cvb,cvbdet)
kbasiscvb = kbasis
nvb = nvb_cvb(kbasiscvb)
call vb2strc_cvb(cvbdet,cvb)

return

end subroutine mktrnspn_cvb
