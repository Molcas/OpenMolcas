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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine CnstAntiC(DPT2Canti,UEFF,U0)

use sguga, only: SGS
use caspt2_global, only: iRoot1, iRoot2, OLagFull
use caspt2_module, only: ENERGY, NASH, NASHT, NBAS, NBAST, NBSQT, NCONF, NDEL, NFRO, NISH, NORB, NSTATE, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: DPT2Canti(*)
real(kind=wp), intent(in) :: UEFF(nState,nState), U0(*)
integer(kind=iwp) :: i, iMO1, iMO2, iOrb0, iOrb2, iStat, iSym, j, jOrb0, jOrb2, jStat, nLev, nOrbI1, nOrbI2
real(kind=wp) :: Scal
real(kind=wp), allocatable :: CI1(:), CI2(:), G1(:,:), SGM1(:), SGM2(:), TG1(:), WRK1(:), WRK2(:)

nLev = SGS%nLev

call mma_allocate(CI1,nConf,Label='CI1')
call mma_allocate(CI2,nConf,Label='CI2')
call mma_allocate(SGM1,nConf,Label='SGM1')
call mma_allocate(SGM2,nConf,Label='SGM2')
call mma_allocate(TG1,nAshT**2,Label='TG1')
call mma_allocate(G1,nAshT,nAshT,Label='G1')

G1(:,:) = Zero
do iStat=1,nState
  !! UEFF is in the CASSCF basis, so the CI coefficient
  !! has to be transformed(-back) accordingly
  call LoadCI_XMS('N',1,nConf,nState,CI1,iStat,U0)
  do jStat=1,nState
    if (iStat == jStat) cycle
    call LoadCI_XMS('N',1,nConf,nState,CI2,jStat,U0)

    call Dens1T_RPT2(CI1,CI2,SGM1,TG1,nLev)
    Scal = UEFF(iStat,iRoot2)*UEFF(jStat,iRoot1)-UEFF(jStat,iRoot2)*UEFF(iStat,iRoot1)
    Scal = Scal*Half
    !call sqprt(tg1,5)
    call DaXpY_(nAshT**2,Scal,TG1,1,G1,1)
  end do
end do
!call sqprt(g1,5)

call mma_deallocate(CI1)
call mma_deallocate(CI2)
call mma_deallocate(SGM1)
call mma_deallocate(SGM2)
call mma_deallocate(TG1)

!! The DPT2Canti computed so far has been doubled
DPT2Canti(1:nBasT**2) = DPT2Canti(1:nBasT**2)*Half

iMO1 = 1
iMO2 = 1
do iSym=1,nSym
  nOrbI1 = nOrb(iSym)
  nOrbI2 = nBas(iSym)-nDel(iSym)
  if (nOrbI2 > 0) then
    !! Add active orbital density
    !! Probably incorrect if symmetry
    do iOrb0=1,nAsh(iSym)
      !iOrb1 = nIsh(iSym)+iOrb0
      iOrb2 = nFro(iSym)+nIsh(iSym)+iOrb0
      do jOrb0=1,nAsh(iSym)
        !jOrb1 = nIsh(iSym)+jOrb0
        jOrb2 = nFro(iSym)+nIsh(iSym)+jOrb0
        DPT2Canti(iMO2+iOrb2-1+nOrbI2*(jOrb2-1)) = DPT2Canti(iMO2+iOrb2-1+nOrbI2*(jOrb2-1))+G1(iOrb0,jOrb0)
      end do
    end do
  end if
  iMO1 = iMO1+nOrbI1*nOrbI1
  iMO2 = iMO2+nOrbI2*nOrbI2
end do
call mma_deallocate(G1)

!! Add orbital response
call mma_allocate(WRK1,NBSQT,Label='WRK1')
call mma_allocate(WRK2,NBSQT,Label='WRK2')

!! Scale with the CASPT2 energy difference
Scal = ENERGY(iRoot2)-ENERGY(iRoot1)
WRK1(1:nBasT**2) = Scal*DPT2Canti(1:nBasT**2)
!! anti-symmetrize the orbital response
call DGeSub(WRK1,nBas(1),'N',WRK1,nBas(1),'T',WRK2,nBas(1),nBas(1),nBas(1))

!! Probably, the way CSF term is computed in MOLCAS is different
!! from in BAGEL; see Eqs.(51)--(53). In the active space, i.e.
!! for SA-CASSCF, the orbital response in the active space cancel
!! (see JCP 2004, 120, 7322), but not in off-diagonal blocks.
!! The way MOLCAS computes adds more than the one BAGEL does, so
!! the orbital response has to be subtracted?
OLagFull(1:nBasT**2) = OLagFull(1:nBasT**2)-WRK1(1:nBasT**2)

call mma_deallocate(WRK1)
call mma_deallocate(WRK2)

!! S[x]^A:D = S[x]:D^A, according to alaska_util/csfgrad
!! so antisymmetrize D
do i=1,nBasT
  do j=1,i-1
    Scal = DPT2Canti(i+nBasT*(j-1))-DPT2Canti(j+nBasT*(i-1))
    DPT2Canti(i+nBasT*(j-1)) = Scal*Half
    DPT2Canti(j+nBasT*(i-1)) = -Scal*Half
  end do
end do

end subroutine CnstAntiC
