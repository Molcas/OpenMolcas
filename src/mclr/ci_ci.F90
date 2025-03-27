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

use ipPage, only: W
use Arrays, only: FIMO, INT2
use input_mclr, only: nRoots, rIn_Ene, PotNuc, ERASSCF, NCSF, Weight, State_Sym
use Constants, only: Two

implicit none
integer ipCID, ipS2
integer irc, i
integer, external :: ipIn
real*8 rDum(1), EC
!                                                                      *
!***********************************************************************
!                                                                      *
interface
  subroutine CISigma_sa(iispin,iCsym,iSSym,Int1,nInt1,Int2s,nInt2s,Int2a,nInt2a,ipCI1,ipCI2,Have_2_el)
    integer iispin, iCsym, iSSym
    integer nInt1, nInt2s, nInt2a
    real*8, target :: Int1(nInt1), Int2s(nInt2s), Int2a(nInt2a)
    integer ipCI1, ipCI2
    logical Have_2_el
  end subroutine CISigma_sa
end interface
!                                                                      *
!***********************************************************************
!                                                                      *

call CISigma_sa(0,state_sym,state_sym,FIMO,size(FIMO),Int2,size(Int2),rDum,1,ipCId,ips2,.true.)
irc = ipin(ipCId)
irc = ipin(ipS2)
do i=0,nroots-1
  EC = (rin_ene+potnuc-ERASSCF(i+1))*Weight(i+1)
  call Daxpy_(ncsf(State_Sym),EC,W(ipCId)%Vec(1+i*ncsf(state_sym)),1,W(ipS2)%Vec(1+i*ncsf(state_sym)),1)
end do
call DSCAL_(nroots*ncsf(state_SYM),Two,W(ipS2)%Vec,1)
#ifdef _WARNING_WORKAROUND_
if (.false.) call Unused_integer(irc)
#endif

end subroutine Ci_Ci
