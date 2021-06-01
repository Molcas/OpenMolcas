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

use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nFro(nSym), nIsh(nSym), nSsh(nSym), nDel(nSym)
real(kind=wp), intent(in) :: CMO(*)
real(kind=wp), intent(inout) :: Eocc(*), EVir(*)
real(kind=wp), intent(out) :: TrD(nSym)
integer(kind=iwp) :: iCMO, iOff, ip_X, ip_Y, irc, iSkip, iSym, iV, kfr, kto, lnDel(8), lnFro(8), lnOcc(8), lnOrb(8), lnVir(8), &
                     nBB, nOA, nVV
real(kind=wp) :: Dummy
real(kind=wp), external :: ddot_
#include "WrkSpc.fh"

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

call GetMem('Dmat','Allo','Real',ip_X,nVV+nOA)
ip_Y = ip_X+nVV
call FZero(Work(ip_X),nVV+nOA)

call LovMP2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,ip_X,ip_Y,.true.)
call GetMem('CMON','Allo','Real',iCMO,nBB)
call FZero(Work(iCMO),nBB)
iOff = 0
do iSym=1,nSym
  kfr = 1+iOff+nBas(iSym)*nFro(iSym)
  kto = iCMO+iOff+nBas(iSym)*lnFro(iSym)
  call dcopy_(nBas(iSym)*lnOcc(iSym),CMO(kfr),1,Work(kto),1)
  kfr = 1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
  kto = kto+nBas(iSym)*lnOcc(iSym)
  call dcopy_(nBas(iSym)*lnVir(iSym),CMO(kfr),1,Work(kto),1)
  iOff = iOff+nBas(iSym)**2
end do

call Check_Amp2(nSym,lnOcc,lnVir,iSkip)
if (iSkip > 0) then
  call ChoMP2_Drv(irc,Dummy,Work(iCMO),EOcc,EVir)
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
call GetMem('CMON','Free','Real',iCMO,nBB)

iV = ip_X
do iSym=1,nSym
  TrD(iSym) = ddot_(lnVir(iSym),Work(iV),1+lnVir(iSym),[One],0)
  iV = iV+lnVir(iSym)**2
end do
call GetMem('Dmat','Free','Real',ip_X,nVV+nOA)

return

end subroutine Get_Tr_Dab
