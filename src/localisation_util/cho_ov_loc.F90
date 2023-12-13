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

use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(in) :: Thrs, SMAT(*)
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nFro(nSym), nIsh(nSym), nAsh(nSym), nSsh(nSym)
real(kind=wp), intent(inout) :: CMO(*)
integer(kind=iwp), intent(_OUT_) :: iD_vir(*)
integer(kind=iwp) :: i, ii, iSym, jD, kOff, l_Dens, nOcc
real(kind=wp) :: yNrm
type(DSBA_Type) :: C, S
real(kind=wp), allocatable :: D2(:), Dens(:)

irc = 0
l_Dens = 0
do iSym=1,nSym
  l_Dens = max(l_Dens,nBas(iSym)**2)
end do
call mma_allocate(Dens,l_Dens,label='Density')
call mma_allocate(D2,l_Dens,label='D2')
call Allocate_DT(C,nBas,nBas,nSym,label='C',Ref=CMO)
call Allocate_DT(S,nBas,nBas,nSym,label='S',Ref=SMAT)
jD = 1
do iSym=1,nSym
  if (nIsh(iSym) > 0) then
    kOff = nFro(iSym)+1
    call GetDens_Localisation(Dens,C%SB(iSym)%A2(:,kOff:),nBas(iSym),nIsh(iSym))
    C%SB(iSym)%A2(:,kOff:kOff+nIsh(iSym)-1) = Zero
    call ChoLoc(irc,Dens,C%SB(iSym)%A2(:,kOff:),Thrs,yNrm,nBas(iSym),nIsh(iSym))
    if (irc /= 0) then
      irc = 1
      exit
    end if
  end if
  iD_vir(jD:jD+nBas(iSym)-1) = 0
  if (nSsh(iSym) > 0) then
    nOcc = nFro(iSym)+nIsh(iSym)+nAsh(iSym)
    call GetDens_Localisation(Dens,C%SB(iSym)%A2,nBas(iSym),nOcc)
    if (nOcc+nSsh(iSym) < nBas(iSym)) then  ! nDel > 0
      write(u6,*) ' ******************************************'
      write(u6,*) ' Cho_ov_Loc found Deleted orbitals in your '
      write(u6,*) ' original MOs. She cannot properly handle  '
      write(u6,*) ' this situation. The program may crash !! '
      write(u6,*) ' ******************************************'
    end if
    ! compute -DS
    call DGEMM_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),-One,Dens,nBas(iSym),S%SB(iSym)%A2,nBas(iSym),Zero,D2,nBas(iSym))
    ! compute 1-DS = 1 + (-DS)
    do i=0,nBas(iSym)-1
      ii = 1+nBas(iSym)*i+i
      D2(ii) = One+D2(ii)
    end do
    ! compute (1-DS)*(1-DS)'
    call GetDens_Localisation(Dens,D2,nBas(iSym),nBas(iSym))
    kOff = nFro(iSym)+nIsh(iSym)+nAsh(iSym)+1
    C%SB(iSym)%A2(:,kOff:kOff+nSsh(iSym)-1) = Zero
    call ChoLoc_xp(irc,Dens,C%SB(iSym)%A2(:,kOff:),Thrs,yNrm,nBas(iSym),nSsh(iSym),iD_vir(jD))
    if (irc /= 0) then
      irc = 1
      exit
    end if
  end if
  jD = jD+nBas(iSym)
end do

call mma_deallocate(Dens)
call mma_deallocate(D2)
call Deallocate_DT(C)
call Deallocate_DT(S)

return

end subroutine Cho_ov_Loc
