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

subroutine CISigma_sa(iispin,iCsym,iSSym,Int1,nInt1,Int2s,nInt2s,Int2a,nInt2a,ipCI1,ipCI2,Have_2_el)
! sigma = H C for each root (weighted with the weight of the root) through
! the CI interface (CI_MCLR_Sigma); see CISigma. With paging only the first
! root is done, unweighted, as before.

use ipPage, only: ipin, ipnout, opout, W
use MCLR_Data, only: nConf1
use CI_Interfaces_MCLR, only: CI_MCLR_Sigma
use input_mclr, only: nCSF, nRoots, Page, Weight
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), parameter :: iState = 1
integer(kind=iwp), intent(in) :: iiSpin, iCSym, iSSym, nInt1, nInt2s, nInt2a, ipCI1, ipCI2
real(kind=wp), target, intent(in) :: Int1(nInt1), Int2s(nInt2s), Int2a(nInt2a)
logical(kind=iwp), intent(in) :: Have_2_el
integer(kind=iwp) :: i, nC, nS

#include "macros.fh"
unused_var(Int2a)

! yma: notice the nconf1 if DMRG
if (nconf1 == 0) return

nC = nCSF(iCSym)
nS = nCSF(iSSym)
if (.not. page) then
  call ipin(ipCI1)
  call ipin(ipci2)
  do i=1,nroots
    if (weight(i) == Zero) then
      W(ipci2)%A((i-1)*nS+1:i*nS) = Zero
    else
      call CI_MCLR_Sigma(iState,iiSpin,iCSym,iSSym,nInt1,Int1,nInt2s,Int2s,Have_2_el,.false.,nC,W(ipCI1)%A((i-1)*nC+1:i*nC),nS, &
                         W(ipci2)%A((i-1)*nS+1:i*nS))
      W(ipci2)%A((i-1)*nS+1:i*nS) = weight(i)*W(ipci2)%A((i-1)*nS+1:i*nS)
    end if
  end do
else
  call ipnout(ipci2)
  call ipin(ipCI1)
  call ipin(ipci2)
  call CI_MCLR_Sigma(iState,iiSpin,iCSym,iSSym,nInt1,Int1,nInt2s,Int2s,Have_2_el,.false.,nC,W(ipCI1)%A(1:nC),nS,W(ipci2)%A(1:nS))
  call opout(ipci1)
end if

end subroutine CISigma_sa
