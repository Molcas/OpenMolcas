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

subroutine Cho_x_Loc(irc,Thrs,nSym,nBas,nFro,nIsh,nAsh,nSsh,CMO,nCMO)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(in) :: Thrs
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nFro(nSym), nIsh(nSym), nAsh(nSym), nSsh(nSym), nCMO
real(kind=wp), intent(inout) :: CMO(nCMO)
integer(kind=iwp) :: iSym, kOff1, kOffC, l_Dens
real(kind=wp) :: yNrm
real(kind=wp), allocatable :: Dens(:)

irc = 0
l_Dens = 0
do iSym=1,nSym
  l_Dens = max(l_Dens,nBas(iSym)**2)
end do
call mma_allocate(Dens,l_Dens,Label='Dens')
kOffC = 0
do iSym=1,nSym
  if (nIsh(iSym) > 0) then
    kOff1 = 1+kOffC+nBas(iSym)*nFro(iSym)
    call GetDens_Localisation(Dens,CMO(kOff1),nBas(iSym),nIsh(iSym))
    CMO(kOff1:kOff1+nBas(iSym)*nIsh(iSym)-1) = Zero
    call ChoLoc(irc,Dens,CMO(kOff1),Thrs,yNrm,nBas(iSym),nIsh(iSym))
    if (irc /= 0) then
      call mma_deallocate(Dens)
      irc = 1
      return
    end if
  end if
  if (nSsh(iSym) > 0) then
    kOff1 = 1+kOffC+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
    call GetDens_Localisation(Dens,CMO(kOff1),nBas(iSym),nSsh(iSym))
    CMO(kOff1:kOff1+nBas(iSym)*nSsh(iSym)-1) = Zero
    call ChoLoc(irc,Dens,CMO(kOff1),Thrs,yNrm,nBas(iSym),nSsh(iSym))
    if (irc /= 0) then
      call mma_deallocate(Dens)
      irc = 1
      return
    end if
  end if
  kOffC = kOffC+nBas(iSym)**2
end do

call mma_deallocate(Dens)

end subroutine Cho_x_Loc
