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

subroutine Compute_Tr_Dab(nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,CMO,nCMO,OrbE,nOrbE,TrD)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nFro(nSym), nIsh(nSym), nAsh(nSym), nSsh(nSym), nDel(nSym), nCMO, nOrbE
real(kind=wp), intent(in) :: CMO(nCMO), OrbE(nOrbE)
real(kind=wp), intent(out) :: TrD(nSym)
integer(kind=iwp) :: iE, ifr, ioff, ip_Y, irc, iSkip, iSym, ito, iV, joff, kEOcc, kEVir, kfr, koff, kto, lnDel(8), lnFro(8), &
                     lnOcc(8), lnOrb(8), lnVir(8), nAct(8), nBB, nOA, nOrb, nVV
real(kind=wp) :: Dummy
real(kind=wp), allocatable :: CMON(:), DMat(:), EOrb(:)

nAct(:) = 0
nVV = sum(nSsh(:)**2)
nOrb = 0
do iSym=1,nSym
  iE = nOrb+nFro(iSym)+nIsh(iSym)
  nAct(iSym) = count(OrbE(iE+1:iE+nAsh(iSym)) < Zero)
  nOrb = nOrb+nBas(iSym)
end do

! setup info
lnOrb(1:nSym) = nBas(:)
lnFro(1:nSym) = nFro(:)
lnOcc(1:nSym) = nIsh(:)+nAct(1:nSym)
lnVir(1:nSym) = nSsh(:)
lnDel(1:nSym) = nDel(:)
nBB = sum(nBas(:)**2)
nOA = sum(lnOcc(1:nSym))

call mma_allocate(Eorb,2*nOrb,Label='EOrb')
kEOcc = 1
kEVir = kEOcc+nOrb
ioff = 0
joff = 0
koff = 0
do iSym=1,nSym
  ifr = 1+ioff+nFro(iSym)
  ito = kEOcc+joff
  EORb(ito:ito+lnOcc(iSym)-1) = OrbE(ifr:ifr+lnOcc(iSym)-1)
  ifr = 1+ioff+nFro(iSym)+nIsh(iSym)+nAsh(iSym)
  ito = kEVir+koff
  EORb(ito:ito+nSsh(iSym)-1) = OrbE(ifr:ifr+nSsh(iSym)-1)
  ioff = ioff+nBas(iSym)
  joff = joff+lnOcc(iSym)
  koff = koff+nSsh(iSym)
end do

call mma_allocate(Dmat,nVV+nOA,Label='DMat')
ip_Y = 1+nVV
DMAT(:) = Zero

call LovCASPT2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,.true.)
call mma_allocate(CMON,nBB,Label='CMON')
CMON(:) = Zero
iOff = 0
do iSym=1,nSym
  kfr = 1+iOff+nBas(iSym)*nFro(iSym)
  kto = 1+iOff+nBas(iSym)*lnFro(iSym)
  CMON(kto:kto+nBas(iSym)*lnOcc(iSym)-1) = CMO(kfr:kfr+nBas(iSym)*lnOcc(iSym)-1)
  kfr = 1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
  kto = kto+nBas(iSym)*lnOcc(iSym)
  CMON(kto:kto+nBas(iSym)*lnVir(iSym)-1) = CMO(kfr:kfr+nBas(iSym)*lnVir(iSym)-1)
  iOff = iOff+nBas(iSym)**2
end do

call Check_Amp(nSym,lnOcc,lnVir,iSkip)
if (iSkip > 0) then
  call ChoMP2_Drv(irc,Dummy,CMON,EOrb(kEOcc),Eorb(kEVir),DMAT(1:nVV),DMAT(ip_Y:))
  if (irc /= 0) then
    write(u6,*) 'MP2 pseudodensity calculation failed !'
    call Abend()
  end if
else
  write(u6,*)
  write(u6,*) 'There are ZERO amplitudes T(ai,bj) with the given '
  write(u6,*) 'combinations of inactive and virtual orbitals !! '
  write(u6,*) 'Check your input and rerun the calculation! Bye!!'
  call Abend()
end if
call mma_deallocate(CMON)

iV = 0
do iSym=1,nSym
  TrD(iSym) = sum(DMat(iV+1:iV+lnVir(iSym)**2:1+lnVir(iSym)))
  iV = iV+lnVir(iSym)**2
end do
call mma_deallocate(Dmat)
call mma_deallocate(Eorb)

end subroutine Compute_Tr_Dab
