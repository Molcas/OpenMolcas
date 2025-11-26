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

subroutine Make_Conn(F,Kappa,P,D)

use MCLR_Data, only: F0SQMO, ipMat, n2Dens, nDens
use PCM_grad, only: do_RF, PCM_grad_TimesE2
use input_mclr, only: nBas, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: F(*)
real(kind=wp), intent(in) :: Kappa(*), P(*), D(*)
integer(kind=iwp) :: iB, idSym, ijB, ipdum, ipTmp, ipTmp1, ipTmp2, iS, jB
real(kind=wp) :: dum(1), Fact
real(kind=wp), allocatable :: T1(:), T2(:), T3(:), T4(:)

! kappa=\bar{kappa}
! P = \bar{d}
! D = \bar{D}

call mma_allocate(T1,n2Dens,Label='MO')
call mma_allocate(T2,nDens,Label='F1')
call mma_allocate(T3,nDens,Label='F3')
call mma_allocate(T4,nDens,Label='F2')

! OIT of the Fock matrix --> Ftilde the orbital part of the effective Fock
! F = Ftilde, the active part.
! T1 = Dtilde

call Rint_generic(kappa,T1,dum,T2,T3,T4,F,1,-One,0)
if (do_RF) then
  idSym = 1
  ipdum = 1
  !! consider both MO and CI parts
  call PCM_grad_TimesE2(idSym,Kappa,F,ipdum)
end if
F(1:nDens) = Half*F(1:nDens)

! T2 = Fbar The ci part or the active Fock matrix

call FockGen(Zero,D,P,T2,T3,1)
Fact = One
do iS=1,nsym
  if (nBas(is) >= 1) then
    call DGEMM_('N','N',nBas(is),nBas(is),nBas(is),Fact,kappa(ipMat(is,is)),nBas(is),F0SQMO(ipMat(is,is)),nBas(is),One, &
                F(ipMat(is,is)),nBas(is))
    call DGEMM_('N','N',nBas(is),nBas(is),nBas(is),-Fact,F0SQMO(ipMat(is,is)),nBas(is),Kappa(ipmat(is,is)),nBas(is),One, &
                F(ipMat(is,is)),nBas(is))
  end if
end do
T2(:) = T2(:)+F(1:nDens)

call TCMO(T2,1,-2)
ijB = 1
do iS=1,nsym
  do iB=1,nBas(iS)
    iptmp = ipmat(iS,iS)
    iptmp1 = iptmp-nbas(iS)+iB-1
    iptmp2 = iptmp+nbas(iS)*(iB-1)-1
    do jB=1,iB-1
      iptmp1 = iptmp1+nbas(iS)
      iptmp2 = iptmp2+1
      F(ijB) = (T2(iptmp1)+T2(iptmp2))
      ijB = ijB+1
    end do
    F(ijB) = T2(iptmp2+1)
    ijB = ijB+1
  end do
end do

call mma_deallocate(T4)
call mma_deallocate(T3)
call mma_deallocate(T2)
call mma_deallocate(T1)

end subroutine Make_Conn
