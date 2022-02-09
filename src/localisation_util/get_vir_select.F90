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

subroutine Get_Vir_Select(irc,CMO,XMO,Eorb,Smat,BName,NamAct,ind_V,nSym,nActa,mOrb,nBas,ortho,n_OK)

use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(in) :: ind_V(*), nSym, nActa, mOrb(nSym), nBas(nSym)
integer(kind=iwp), intent(out) :: irc, n_OK(nSym)
real(kind=wp), intent(inout) :: CMO(*), Eorb(*)
real(kind=wp), intent(in) :: XMO(*), Smat(*)
character(len=LenIn8), intent(in) :: BName(*)
character(len=LenIn), intent(in) :: NamAct(nActa)
logical(kind=iwp), intent(in) :: ortho
integer(kind=iwp) :: i, iOff, iSym, j, ja, km, kx, lOff, mOx, n_KO, nBmx, nBx, nOrbmx, nOx
character(len=LenIn) :: tmp
type(DSBA_Type) :: C, S, X
real(kind=wp), allocatable :: Fock(:,:), U(:,:)
real(kind=wp), allocatable, target :: Scr1(:,:)
real(kind=wp), pointer :: CC(:,:), Scr(:,:), X2(:,:), Z(:,:)
integer(kind=iwp), allocatable :: iD(:,:)

irc = 0

nBmx = 0
nOrbmx = 0
do iSym=1,nSym
  nBmx = max(nBmx,nBas(iSym))
  nOrbmx = max(nOrbmx,mOrb(iSym))
end do
call mma_allocate(iD,nOrbmx,2,label='iD')
call mma_allocate(Scr1,nBmx*nOrbmx,4,label='Scr1')
call mma_allocate(U,nOrbmx,nOrbmx,label='U')
call mma_allocate(Fock,nOrbmx,2,label='Fock')

call Allocate_DT(C,nBas,mOrb,nSym,label='C',Ref=CMO)
call Allocate_DT(X,nBas,mOrb,nSym,label='X',Ref=XMO)
call Allocate_DT(S,nBas,nBas,nSym,label='S',Ref=Smat)

iOff = 0
lOff = 0
do iSym=1,nSym

  CC(1:mOrb(iSym),1:nBas(iSym)) => Scr1(1:mOrb(iSym)*nBas(iSym),1)
  X2(1:nBas(iSym),1:mOrb(iSym)) => Scr1(1:nBas(iSym)*mOrb(iSym),2)
  Z(1:nBas(iSym),1:mOrb(iSym)) => Scr1(1:nBas(iSym)*mOrb(iSym),3)
  Scr(1:nBas(iSym),1:mOrb(iSym)) => Scr1(1:nBas(iSym)*mOrb(iSym),4)

  n_OK(iSym) = 0
  n_KO = 0
  do i=1,mOrb(iSym)
    ja = iOff+ind_V(i+iOff)
    tmp = BName(ja)(1:LenIn)
    kx = 0

    !write(u6,*) ' We simulate Afreeze with all Vir'
    !kx = 1

    do j=1,nActa
      if (NamAct(j) == tmp) kx = kx+1
    end do
    if (kx > 0) then
      n_OK(iSym) = n_OK(iSym)+1
      X2(:,n_OK(iSym)) = X%SB(iSym)%A2(:,i)
      iD(n_OK(iSym),1) = i
    else
      n_KO = n_KO+1
      Z(:,n_KO) = X%SB(iSym)%A2(:,i)
      iD(n_KO,2) = i
    end if
  end do

  if (.not. ortho) then

    call Ortho_orb(X2,S%SB(iSym)%A2,nBas(iSym),n_OK(iSym),2,.false.)
    call Ortho_orb(Z,S%SB(iSym)%A2,nBas(iSym),n_KO,2,.false.)
  end if

  nBx = max(1,nBas(iSym))
  mOx = max(1,mOrb(iSym))

  call DGEMM_('T','N',mOrb(iSym),nBas(iSym),nBas(iSym),One,C%SB(iSym)%A2,nBx,S%SB(iSym)%A2,nBx,Zero,CC,mOx)

  if (n_KO > 0) then
    call DGEMM_('N','N',mOrb(iSym),n_KO,nBas(iSym),One,CC,mOx,Z,nBx,Zero,U,mOx)

    call Get_Can_Lorb(Eorb(lOff+1),Fock(:,1),n_KO,mOrb(iSym),iD(:,2),U)

    call DGEMM_('N','N',nBas(iSym),n_KO,n_KO,One,Z,nBx,U,n_KO,Zero,Scr,nBx)

    ! Reorder the final MOs such that those of the active site come first
    km = n_OK(iSym)+1
    C%SB(iSym)%A2(:,km:km+n_KO-1) = Scr(:,1:n_KO)
    Fock(:,2) = Fock(:,1)
  end if

  call DGEMM_('N','N',mOrb(iSym),n_OK(iSym),nBas(iSym),One,CC,mOx,X2,nBx,Zero,U,mOx)

  call Get_Can_Lorb(Eorb(lOff+1),Fock(:,1),n_OK(iSym),mOrb(iSym),iD(:,1),U)

  nOx = max(1,n_OK(iSym))
  call DGEMM_('N','N',nBas(iSym),n_OK(iSym),n_OK(iSym),One,X2,nBx,U,nOx,Zero,Scr,nBx)

  do i=1,n_OK(iSym)
    Eorb(lOff+i) = Fock(iD(i,1),1)
  end do
  C%SB(iSym)%A2(:,1:n_OK(iSym)) = Scr(:,1:n_OK(iSym))

  do i=1,n_KO
    Eorb(lOff+n_OK(iSym)+i) = Fock(iD(i,2),2)
  end do

  iOff = iOff+nBas(iSym)
  lOff = lOff+mOrb(iSym)
end do

nullify(CC)
nullify(X2)
nullify(Z)
nullify(Scr)

call mma_deallocate(iD)
call mma_deallocate(Scr1)
call mma_deallocate(U)
call mma_deallocate(Fock)
call Deallocate_DT(C)
call Deallocate_DT(X)
call Deallocate_DT(S)

return

end subroutine Get_Vir_Select
