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
use MCLR_Data, only: FIMO, INT2
use MCLR_procedures, only: CISigma_sa
use input_mclr, only: nRoots, rIn_Ene, PotNuc, ERASSCF, NCSF, Weight, State_Sym
use Constants, only: Two

implicit none
integer ipCID, ipS2
integer i
real*8 rDum(1), EC

call CISigma_sa(0,state_sym,state_sym,FIMO,size(FIMO),Int2,size(Int2),rDum,1,ipCId,ips2,.true.)
call ipin(ipCId)
call ipin(ipS2)
do i=0,nroots-1
  EC = (rin_ene+potnuc-ERASSCF(i+1))*Weight(i+1)
  call Daxpy_(ncsf(State_Sym),EC,W(ipCId)%A(1+i*ncsf(state_sym)),1,W(ipS2)%A(1+i*ncsf(state_sym)),1)
end do
call DSCAL_(nroots*ncsf(state_SYM),Two,W(ipS2)%A,1)

end subroutine Ci_Ci
