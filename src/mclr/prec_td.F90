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

subroutine Prec_td(pre2,DigPrec,isym)
! pre2      Preconditioner from Prec
! DigPrec Output - Diagonal of prec2
! isym      Symmetry of PT

use Index_Functions, only: iTri
use Symmetry_Info, only: Mul
use MCLR_Data, only: G1t, ipCM, ipMat, nA, nDens
use input_mclr, only: nAsh, nBas, nIsh, nSym, Omega
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: pre2(*)
real(kind=wp), intent(_OUT_) :: DigPrec(*)
integer(kind=iwp), intent(in) :: iSym
integer(kind=iwp) :: i, iA, iB, Inc, ip, ip1, ip2, ip3, iS, j, jA, jB, jS, k, l, nBasTot, nD
real(kind=wp) :: nonzero
real(kind=wp), allocatable :: Dens(:), PreTd(:), TempTd(:)

!                                                                      *
!***********************************************************************
!                                                                      *

!---------------------------------------------
! Construct one el density in MO Dens
!---------------------------------------------

nBasTot = sum(nBas(1:nSym)**2)
call mma_allocate(Dens,nBasTot,Label='Dens')

Dens(:) = Zero
ip3 = 1
do iS=1,nSym
  inc = nBas(iS)+1
  Dens(ip3:ip3+nIsh(iS)*inc-1:inc) = Two
  ip3 = ip3+nBas(iS)**2
end do

! For a CASSCF wavefunc. From Anders subrout r2elint
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

!call RECPRT('Dens',' ',Dens,nBasTot,1)
!write(u6,*) 'Diagnonal elements in D'
!do iS=1,nSym
!  do k=0,nBas(iS)-1
!    write(u6,*) Dens(ipCM(iS)+k*(nBas(iS)+1))
!  end do
!end do
!stop

!-------------------------------------------------------------------
! Construct the diagonal approximation to the orbital prec, PreTd
!-------------------------------------------------------------------

call mma_allocate(PreTd,nDens,Label='PreTd')
PreTd(:) = Zero
ip1 = 1
ip2 = 1
!ipsave = 0
do iS=1,nSym
  jS = Mul(iS,iSym)
  nD = nBas(jS)-nIsh(jS)
  do k=1,nIsh(iS)
    ip1 = ip1+nIsh(jS)
    PreTd(ip1:ip1+nD-1) = Pre2(ip2:ip2+nD-1)
    ip1 = ip1+nD
    ip2 = ip2+nD
  end do
  nD = nBas(jS)-nAsh(jS)
  do k=1,nAsh(iS)
    PreTd(ip1+1:ip1+nIsh(jS)) = Pre2(ip2+1:ip2+nIsh(jS))
    ip1 = ip1+nIsh(jS)+nAsh(jS)
    ip2 = ip2+nIsh(jS)
    PreTd(ip1+1:ip1+nD-nIsh(jS)) = Pre2(ip2+1:ip2+nD-nIsh(jS))
    ip1 = ip1+nD-nIsh(jS)
    ip2 = ip2+nD-nIsh(jS)
    if ((nBas(jS)-nAsh(jS)-nIsh(jS)) == 0) ip1 = ip1+nAsh(jS)
  end do
  !call RECPRT('PreTd',' ',PreTd(1+ipsave),nBas(jS),nBas(iS))
  ip1 = ip1+(nBas(iS)-nIsh(iS)-nAsh(iS))*nBas(jS)
  !ipsave = ip1
end do
!call RECPRT('PreTd',' ',PreTd,nDens,1)

!----------------------
! Symmetrize PreTd
!----------------------
call mma_allocate(TempTd,nDens,Label='TempTd')

do iS=1,nSym
  jS = Mul(iS,iSym)
  TempTd(:) = Zero
  call Trnsps(nBas(iS),nBas(jS),PreTd(ipMat(jS,iS)),TempTd)
  nD = nBas(iS)*nBas(jS)
  do i=0,nD-1
    nonzero = PreTd(ipMat(iS,jS)+i)
    if (nonzero /= Zero) TempTd(1+i) = PreTd(ipMat(iS,jS)+i)
  end do
  call Trnsps(nBas(jS),nBas(iS),TempTd,PreTd(ipMat(jS,iS)))

end do

!------------------------------------------------------------------
! Add the density part PreTd_at = PreTd_at - omega(D_aa - D_tt)
!------------------------------------------------------------------
i = 0
do iS=1,nSym
  jS = Mul(iS,iSym)
  nD = nBas(iS)*nBas(jS)
  j = 0
  l = 0
  do k=0,nD-1
    if (l == nBas(jS)) l = 0
    if (k == (j+1)*nBas(jS)) j = j+1

    i = i+1
    PreTd(i) = PreTd(i)+Two*Omega*(Dens(ipCM(iS)+j*(nBas(iS)+1))+Dens(ipCM(jS)+l*(nBas(jS)+1)))

    l = l+1
  end do

end do

!-----------------------------------------------------------------------
! Symmetry transpose PreTd - To get the same order fo sym as in b_x and
! as required by compress.
!-----------------------------------------------------------------------

TempTd(:) = Zero

do iS=1,nSym
  jS = Mul(iS,iSym)
  nD = nBas(iS)*nBas(jS)
  do k=0,nD-1
    TempTd(ipmat(iS,jS)+k) = PreTd(ipmat(jS,iS)+k)
  end do
end do

call Compress(TempTd,DigPrec,isym)

call mma_deallocate(TempTd)
call mma_deallocate(PreTd)
call mma_deallocate(Dens)

end subroutine Prec_td
