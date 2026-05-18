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

use definitions, only: iwp, wp
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nFro(nSym), nAsh(nSym), nIsh(nSym), nSsh(nSym)
real(kind=wp), intent(in) :: Thrs
integer(kind=iwp), intent(in) :: nCMO
real(kind=wp), intent(inout) :: CMO(nCMO)
real(kind=wp), allocatable :: Dens(:)
integer(kind=iwp) irc, iSym, kOff1, kOffC, l_Dens
real(kind=wp) yNrm

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
    call FZero(CMO(kOff1),nBas(iSym)*nIsh(iSym))
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
    call FZero(CMO(kOff1),nBas(iSym)*nSsh(iSym))
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
