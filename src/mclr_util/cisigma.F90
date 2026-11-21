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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine CISigma(iispin,iCsym,iSSym,Int1,nInt1,Int2s,nInt2s,Int2a,nInt2a,ipCI1,ipCI2,Have_2_el)
! sigma = H C through the CI interface (CI_MCLR_Sigma), for the vectors
! ipCI1 (symmetry iCSym) and ipCI2 (symmetry iSSym); iiSpin: 0 singlet,
! 1 triplet operator. Int2a is not used (it was for the retired spin
! polarization).

use ipPage, only: ipin, ipnout, opout, W
use MCLR_Data, only: nConf1
use CI_Interfaces_MCLR, only: CI_MCLR_Sigma
use input_mclr, only: nCSF, Page
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), parameter :: iState = 1
integer(kind=iwp), intent(in) :: iiSpin, iCSym, iSSym, nInt1, nInt2s, nInt2a, ipCI1, ipCI2
real(kind=wp), target, intent(in) :: Int1(nInt1), Int2s(nInt2s), Int2a(nInt2a)
logical(kind=iwp), intent(in) :: Have_2_el

#include "macros.fh"
unused_var(Int2a)

! yma: notice the nconf1 if DMRG
if (nconf1 == 0) return

if (page) call ipnout(ipci2)
call ipin(ipCI1)
call ipin(ipci2)
call CI_MCLR_Sigma(iState,iiSpin,iCSym,iSSym,nInt1,Int1,nInt2s,Int2s,Have_2_el,.false.,nCSF(iCSym),W(ipCI1)%A(1:nCSF(iCSym)), &
                   nCSF(iSSym),W(ipci2)%A(1:nCSF(iSSym)))
if (page) call opout(ipci1)

end subroutine CISigma
