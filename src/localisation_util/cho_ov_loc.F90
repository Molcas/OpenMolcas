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

subroutine Cho_ov_Loc(irc,Thrs,nSym,nBas,nFro,nIsh,nAsh,nSsh,CMO,SMAT,iD_vir)

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(in) :: Thrs, SMAT(*)
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nFro(nSym), nIsh(nSym), nAsh(nSym), nSsh(nSym)
real(kind=wp), intent(inout) :: CMO(*)
integer(kind=iwp), intent(_OUT_) :: iD_vir(*)
#include "WrkSpc.fh"
integer(kind=iwp) :: i, ii, ip_Dens, ipD2, iSym, jD, kOff1, kOff2, kOffC, l_Dens, nOcc
real(kind=wp) :: yNrm

irc = 0
l_Dens = 0
do iSym=1,nSym
  l_Dens = max(l_Dens,nBas(iSym)**2)
end do
call GetMem('Density','Allo','Real',ip_Dens,2*l_Dens)
ipD2 = ip_Dens+l_Dens
kOffC = 0
jD = 1
do iSym=1,nSym
  if (nIsh(iSym) > 0) then
    kOff1 = 1+kOffC+nBas(iSym)*nFro(iSym)
    call GetDens_Localisation(Work(ip_Dens),CMO(kOff1),nBas(iSym),nIsh(iSym))
    call FZero(CMO(kOff1),nBas(iSym)*nIsh(iSym))
    call ChoLoc(irc,Work(ip_Dens),CMO(kOff1),Thrs,yNrm,nBas(iSym),nIsh(iSym))
    if (irc /= 0) then
      call GetMem('Density','Free','Real',ip_Dens,2*l_Dens)
      irc = 1
      return
    end if
  end if
  call izero(iD_vir(jD),nBas(iSym))
  if (nSsh(iSym) > 0) then
    kOff1 = 1+kOffC
    nOcc = nFro(iSym)+nIsh(iSym)+nAsh(iSym)
    call GetDens_Localisation(Work(ip_Dens),CMO(kOff1),nBas(iSym),nOcc)
    if (nOcc+nSsh(iSym) < nBas(iSym)) then  ! nDel > 0
      write(u6,*) ' ******************************************'
      write(u6,*) ' Cho_ov_Loc found Deleted orbitals in your '
      write(u6,*) ' original MOs. She cannot properly handle  '
      write(u6,*) ' this situation. The program may crash !! '
      write(u6,*) ' ******************************************'
    end if
    ! compute -DS
    call DGEMM_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),-One,Work(ip_Dens),nBas(iSym),SMAT(kOff1),nBas(iSym),Zero,Work(ipD2), &
                nBas(iSym))
    ! compute 1-DS = 1 + (-DS)
    do i=0,nBas(iSym)-1
      ii = ipD2+nBas(iSym)*i+i
      Work(ii) = One+Work(ii)
    end do
    ! compute (1-DS)*(1-DS)'
    call GetDens_Localisation(Work(ip_Dens),Work(ipD2),nBas(iSym),nBas(iSym))
    kOff2 = kOff1+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
    call FZero(CMO(kOff2),nBas(iSym)*nSsh(iSym))
    call ChoLoc_xp(irc,Work(ip_Dens),CMO(kOff2),Thrs,yNrm,nBas(iSym),nSsh(iSym),iD_vir(jD))
    if (irc /= 0) then
      call GetMem('Density','Free','Real',ip_Dens,2*l_Dens)
      irc = 1
      return
    end if
  end if
  kOffC = kOffC+nBas(iSym)**2
  jD = jD+nBas(iSym)
end do

call GetMem('Density','Free','Real',ip_Dens,2*l_Dens)

return

end subroutine Cho_ov_Loc
