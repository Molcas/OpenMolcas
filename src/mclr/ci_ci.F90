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

subroutine Ci_Ci(ipcid,ips2)

use ipPage, only: ipin, W
use MCLR_Data, only: INT2, FIMO
use MCLR_procedures, only: CISigma_sa
use input_mclr, only: ERASSCF, NCSF, nRoots, PotNuc, rIn_Ene, State_Sym, Weight
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ipCID, ipS2
integer(kind=iwp) :: i, n
real(kind=wp) :: EC, rDum(1)

call CISigma_sa(0,state_sym,state_sym,FIMO,size(FIMO),Int2,size(Int2),rDum,1,ipCId,ips2,.true.)
call ipin(ipCId)
call ipin(ipS2)
n = ncsf(State_Sym)
do i=1,nroots
  EC = (rin_ene+potnuc-ERASSCF(i))*Weight(i)
  W(ipS2)%A((i-1)*n+1:i*n) = W(ipS2)%A((i-1)*n+1:i*n)+EC*W(ipCId)%A((i-1)*n+1:i*n)
end do
W(ipS2)%A(1:nroots*ncsf(state_SYM)) = Two*W(ipS2)%A(1:nroots*ncsf(state_SYM))

end subroutine Ci_Ci
