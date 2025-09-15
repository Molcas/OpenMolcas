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
! Copyright (C) 2000, Jonna Stalring                                   *
!***********************************************************************

subroutine calcerr(kappa,iestate)
!-------------------------------------------------
! Jonna 000411
!
! Calculates the derivative d E_j
!                           ------ = 2\sum_pq Fpq kappa_pq
!                            d w
!
! which estimates the error in the SA.
!---------------------------------------------------

use Index_Functions, only: nTri_Elem
use MCLR_Data, only: ipMat, ISTATE, nDens
use input_mclr, only: nBas, ntAsh, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: kappa(nDens)
integer(kind=iwp), intent(in) :: ieState
integer(kind=iwp) :: iS, nG1
real(kind=wp) :: dejdw
real(kind=wp), allocatable :: G1q(:), G1r(:), G2r(:), Q(:), T(:)
real(kind=wp), external :: DDot_

ng1 = nTri_Elem(ntash)

call mma_allocate(G1Q,ng1,Label='G1Q')
call mma_allocate(G1R,ntash**2,Label='G1R')
call mma_allocate(G2R,ntash**4,Label='G2R')
call mma_allocate(T,nDens,Label='T')
call mma_allocate(Q,nDens,Label='Q')

! Get densities of iestate in triangular and rectangular form.

call rddj(G1R,G1Q,G2R,iestate)

! Get Fock matrix ipT

call FockGen(One,G1R,G2R,T,Q,1)

! Print energy

!call prnte(G1q,T)

!call Recprt('KAPP',' ',kappa,nDens,1)
!call Recprt('FOCK',' ',T,nDens,1)

dejdw = Zero
do iS=1,nsym
  dejdw = dejdw+ddot_(nBas(is)**2,T(ipmat(is,is)),1,kappa(ipmat(is,is)),1)
end do
dejdw = Two*dejdw
!is = 1
!write(u6,*) 'T(ipmat(is,is))',T(ipmat(is,is))
!write(u6,*) 'kappa(ipmat(is,is)+1)',kappa(ipmat(is,is)+1)

! d E_j / d w_i, without constraint
!write(u6,200) iestate,istate,dejdw

! d E_i / d w_i, with the constraint sum(w_i)=1,
! multiplied by 1-w_i... which happens to be equal to
! d E_i / d w_i, without constraint
! change sign, because this is the *error*
if (iestate == istate) write(u6,100) iestate,-dejdw

call mma_deallocate(Q)
call mma_deallocate(T)
call mma_deallocate(G2R)
call mma_deallocate(G1R)
call mma_deallocate(G1Q)

return

100 format(' **********'/'                  Estimated error in the energy of state ',I5,': ',ES12.5/' **********')
!200 format(' Derivative of the energy of state ',I5,' with the weight ',I5,': ',ES12.5)

end subroutine calcerr
