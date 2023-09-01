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

subroutine Get_Tr_Dab(nSym,nBas,nFro,nIsh,nSsh,nDel,CMO,EOcc,EVir,TrD)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nFro(nSym), nIsh(nSym), nSsh(nSym), nDel(nSym)
real(kind=wp), intent(in) :: CMO(*)
real(kind=wp), intent(inout) :: Eocc(*), EVir(*)
real(kind=wp), intent(out) :: TrD(nSym)
integer(kind=iwp) :: iOff, irc, iSkip, iSym, iV, kfr, kto, lnDel(8), lnFro(8), lnOcc(8), lnOrb(8), lnVir(8), nBB, nOA, nVV
real(kind=wp) :: Dummy
real(kind=wp), allocatable :: CMON(:), X(:)
real(kind=wp), external :: ddot_

nVV = 0
nBB = 0
nOA = 0
do iSym=1,nSym  ! setup info
  lnOrb(iSym) = nBas(iSym)
  lnFro(iSym) = nFro(iSym)
  lnOcc(iSym) = nIsh(iSym)
  lnVir(iSym) = nSsh(iSym)
  lnDel(iSym) = nDel(iSym)
  nVV = nVV+lnVir(iSym)**2
  nBB = nBB+nBas(iSym)**2
  nOA = nOA+lnOcc(iSym)
end do

call mma_allocate(X,nVV+nOA,label='Dmat')
X(:) = Zero

call LovMP2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,.true.)
call mma_allocate(CMON,nBB,label='CMON')
CMON(:) = Zero
iOff = 1
do iSym=1,nSym
  kfr = iOff+nBas(iSym)*nFro(iSym)
  kto = iOff+nBas(iSym)*lnFro(iSym)
  call dcopy_(nBas(iSym)*lnOcc(iSym),CMO(kfr),1,CMON(kto),1)
  kfr = iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
  kto = kto+nBas(iSym)*lnOcc(iSym)
  call dcopy_(nBas(iSym)*lnVir(iSym),CMO(kfr),1,CMON(kto),1)
  iOff = iOff+nBas(iSym)**2
end do

call Check_Amp2(nSym,lnOcc,lnVir,iSkip)
if (iSkip > 0) then
  call ChoMP2_Drv(irc,Dummy,CMON,EOcc,EVir,X(1:NVV),X(nVV+1:))
  if (irc /= 0) then
    write(u6,*) 'MP2 pseudodensity calculation failed !'
    call Abend()
  end if
else
  write(u6,*)
  write(u6,*) 'There are ZERO amplitudes T(ai,bj) with the given '
  write(u6,*) 'combinations of occupied and virtual orbitals !! '
  write(u6,*) 'Check your input and rerun the calculation! Bye!!'
  call Abend()
end if
call mma_deallocate(CMON)

iV = 1
do iSym=1,nSym
  TrD(iSym) = ddot_(lnVir(iSym),X(iV),1+lnVir(iSym),[One],0)
  iV = iV+lnVir(iSym)**2
end do
call mma_deallocate(X)

return

end subroutine Get_Tr_Dab
