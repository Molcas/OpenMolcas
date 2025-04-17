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
! Copyright (C) Jonna Stalring                                         *
!***********************************************************************

subroutine RInt_td(ekappa,mkappa,isym)
!******************************************************
!            [2]                                      *
! Modifies  E    *kappa by subtracting Omega*S*kappa  *
! to make timedep PT possible.                        *
!                                                     *
!******************************************************
!
!  mkappa    is d/dx(kappa) in matrix form (derivative of orb rot mat w r t PT)
!  ekappa    is E*d/dx(kappa) in matrix form, REPLACED BY
!            the entire timedep expression
!  isym      Rinttd is called once for each symmtry
!
! Local variables
!
!  dens The density matrix
!  wDKt   Omega*(density matrix)*(kappa transposed)
!  wKtD   As above but different order

use Index_Functions, only: iTri
use Symmetry_Info, only: Mul
use MCLR_Data, only: G1t, ipCM, ipMat, nA, nDens
use input_mclr, only: nAsh, nBas, nIsh, nSym, Omega
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: ekappa(nDens)
real(kind=wp), intent(in) :: mkappa(nDens)
integer(kind=iwp), intent(in) :: iSym
integer(kind=iwp) :: iA, iB, ifin, iini, Inc, ip, ip2, ip3, iS, jA, jB, jS, lDens
real(kind=wp), allocatable :: Dens(:), wDKt(:), wKtD(:)

!-----------------------------------------------------------------------

! Allocate memory

lDens = sum(nBas(1:nSym)**2)
call mma_allocate(Dens,lDens,Label='Dens')
call mma_allocate(wDKt,nDens,Label='wDKt')
call mma_allocate(wKtD,nDens,Label='wKtD')

!******************************
! Construct the density matrix
!******************************

! For HF

Dens(:) = Zero
ip3 = 1
do iS=1,nSym
  inc = nBas(iS)+1
  Dens(ip3:ip3+nIsh(iS)*inc-1:inc) = Two
  ip3 = ip3+nBas(iS)*nBas(iS)
end do

!do is=1,nsym
!  size = nbas(is)*nbas(is)      ! The size of the sym i and j martix
!  do k=-1,(size-1)             ! Loop over the entire submatrix
!    residual = MOD((k+1),(nbas(is)+1))
!    if ((residual == zero) .and. (k < (nbas(is)*nIsh(is))) .and. (nIsh(is) /= 0)) then
!      Dens(1+ipCM(is)+k) = Two
!    else
!      Dens(1+ipCM(is)+k) = Zero
!    end if
!  end do
!end do

! For a CASSCF wavefunc. From Anders subrut r2elint
! Add the active active dens

do iS=1,nSym
  do iB=1,nAsh(iS)
    do jB=1,nAsh(iS)
      ip = ipCM(iS)+ib+nIsh(is)+(jB+nIsh(is)-1)*nBas(is)-1
      iA = nA(is)+ib
      jA = nA(is)+jb
      ip2 = iTri(iA,jA)
      Dens(ip) = G1t(ip2)
    end do
  end do
end do
!*************************************
! Multiply D and mkappa wDKt and wKtD
!*************************************
!call RecPrt('dens ',' ',dens,ldens,1)
!call RecPrt('ekappa ',' ',ekappa,nDens,1)
do is=1,nsym
  js = Mul(is,isym)
  ! wDKt
  if ((nBas(iS) > 0) .and. (nBas(jS) > 0)) then
    call DGEMM_('n','n',nbas(is),nbas(js),nbas(is),Two*Omega,Dens(ipCM(is)),nbas(is),mkappa(ipmat(is,js)),nbas(is),Zero, &
                wDKt(ipmat(is,js)),nbas(is))
    ! wKtD
    call DGEMM_('n','n',nbas(is),nbas(js),nbas(js),Two*Omega,mkappa(ipmat(is,js)),nbas(is),Dens(ipCM(js)),nbas(js),Zero, &
                wKtD(ipmat(is,js)),nbas(is))

    !****************************************
    ! Replace ekappa ekappa=ekappa-wDKt+wKtD
    !****************************************
    iini = ipmat(is,js)
    ifin = iini+nbas(is)*nbas(js)-1
    ekappa(iini:ifin) = ekappa(iini:ifin)+wDKt(iini:ifin)-wKtD(iini:ifin)
  end if
end do
!call RecPrt('wDKt ',' ',wDKt,nDens,1)
!call RecPrt('wKtD ',' ',wKtD,nDens,1)
!call RecPrt('ekappa ',' ',ekappa,nDens,1)

! Free memory

call mma_deallocate(Dens)
call mma_deallocate(wDKt)
call mma_deallocate(wKtD)

end subroutine RInt_td
