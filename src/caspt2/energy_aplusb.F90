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
! Copyright (C) 2008, Francesco Aquilante                              *
!***********************************************************************

subroutine Energy_AplusB(nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,CMO,nCMO,OrbE,nOrbE,E2_ab)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nFro(nSym), nIsh(nSym), nAsh(nSym), nSsh(nSym), nDel(nSym), nCMO, nOrbE
real(kind=wp), intent(in) :: CMO(nCMO), OrbE(nOrbE)
real(kind=wp), intent(out) :: E2_ab
integer(kind=iwp) :: iE, ifr, ioff, irc, iSkip, iSym, ito, joff, k, kEOcc, kEVir, kfr, koff, kto, lnDel(8), lnFro(8), lnOcc(8), &
                     lnOrb(8), lnVir(8), nAct(8), nBB, nOA, nOrb, nVV
real(kind=wp) :: Dummy(1)
real(kind=wp), allocatable :: CMOX(:), Eorb(:)

call Izero(nAct,nSym)
nVV = 0
nOrb = 0
do iSym=1,nSym
  iE = 1+nOrb+nFro(iSym)+nIsh(iSym)
  do k=0,nAsh(iSym)-1
    if (OrbE(iE+k) < Zero) nAct(iSym) = nAct(iSym)+1
  end do
  nVV = nVV+nSsh(iSym)**2
  nOrb = nOrb+nBas(iSym)
end do

nBB = 0
nOA = 0
do iSym=1,nSym  ! setup info
  lnOrb(iSym) = nBas(iSym)
  lnFro(iSym) = nFro(iSym)
  lnOcc(iSym) = nIsh(iSym)+nAct(iSym)
  lnVir(iSym) = nSsh(iSym)
  lnDel(iSym) = nDel(iSym)
  nBB = nBB+nBas(iSym)**2
  nOA = nOA+lnOcc(iSym)
end do

call mma_allocate(Eorb,2*nOrb,Label='Eorb')
kEOcc = 1
kEVir = kEOcc+nOrb
ioff = 0
joff = 0
koff = 0
do iSym=1,nSym
  ifr = 1+ioff+nFro(iSym)
  ito = kEOcc+joff
  call dcopy_(lnOcc(iSym),OrbE(ifr),1,Eorb(ito),1)
  ifr = 1+ioff+nFro(iSym)+nIsh(iSym)+nAsh(iSym)
  ito = kEVir+koff
  call dcopy_(nSsh(iSym),OrbE(ifr),1,Eorb(ito),1)
  ioff = ioff+nBas(iSym)
  joff = joff+lnOcc(iSym)
  koff = koff+nSsh(iSym)
end do

call LovCASPT2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,.false.)
call mma_allocate(CMOX,nBB,Label='CMOX')
CMOX(:) = Zero
iOff = 0
do iSym=1,nSym
  kfr = 1+iOff+nBas(iSym)*nFro(iSym)
  kto = 1+iOff+nBas(iSym)*lnFro(iSym)
  call dcopy_(nBas(iSym)*lnOcc(iSym),CMO(kfr),1,CMOX(kto),1)
  kfr = 1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
  kto = kto+nBas(iSym)*lnOcc(iSym)
  call dcopy_(nBas(iSym)*lnVir(iSym),CMO(kfr),1,CMOX(kto),1)
  iOff = iOff+nBas(iSym)**2
end do

call Check_Amp(nSym,lnOcc,lnVir,iSkip)
if (iSkip > 0) then
  call ChoMP2_Drv(irc,E2_ab,CMOX,Eorb(kEOcc),Eorb(kEVir),Dummy,Dummy)
  if (irc /= 0) then
    write(u6,*) 'MP2 calculation failed in energy_AplusB !'
    call Abend()
  end if
else
  write(u6,*)
  write(u6,*) 'There are ZERO amplitudes T(ai,bj) with the given '
  write(u6,*) 'combinations of inactive and virtual orbitals !! '
  write(u6,*) 'Check your input and rerun the calculation! Bye!!'
  call Abend()
end if
call mma_deallocate(CMOX)

call mma_deallocate(Eorb)

end subroutine Energy_AplusB
