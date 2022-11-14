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

subroutine Get_Orb_Select(irc,CMO,XMO,Eorb,Smat,Saa,BName,NamAct,nSym,nActa,mOrb,nBas,ortho,ThrSel,n_OK)

use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(in) :: nSym, nActa, mOrb(nSym), nBas(nSym)
integer(kind=iwp), intent(out) :: irc, n_OK(nSym)
real(kind=wp), intent(inout) :: CMO(*), Eorb(*)
real(kind=wp), intent(in) :: XMO(*), Smat(*), Saa(*), ThrSel
character(len=LenIn8), intent(in) :: BName(*)
character(len=LenIn), intent(in) :: NamAct(nActa)
logical(kind=iwp), intent(in) :: ortho
integer(kind=iwp) :: i, ia, iOff, iSym, j, ja, km, lOff, mOx, n_KO, nBa, nBax, nBmx, nBx, nORbmx, nOx
real(kind=wp) :: ThrS
type(DSBA_Type) :: C, S, X
integer(kind=iwp), allocatable :: iD(:)
real(kind=wp), allocatable :: Fock(:,:), Q(:), U(:,:)
real(kind=wp), allocatable, target :: Scr1(:,:), SQt(:)
real(kind=wp), pointer :: C2(:,:), CC(:,:), Scr(:,:), SQ(:,:), X2(:,:), Z(:,:)
character(len=LenIn) :: tmp
real(kind=wp), external :: ddot_

irc = 0

nBmx = 0
nOrbmx = 0
do iSym=1,nSym
  nBmx = max(nBmx,nBas(iSym))
  nOrbmx = max(nOrbmx,mOrb(iSym))
end do
call mma_allocate(iD,2*nOrbmx+nBmx,label='iD')
call mma_allocate(SQt,nBmx**2,label='Smx')
call mma_allocate(Scr1,nBmx*nOrbmx,5,label='Scr1')
call mma_allocate(U,nOrbmx,nOrbmx,label='U')
call mma_allocate(Q,nOrbmx,label='Q')
call mma_allocate(Fock,nOrbmx,2,label='Fock')
Q(:) = Zero

call Allocate_DT(C,nBas,mOrb,nSym,label='C',Ref=CMO)
call Allocate_DT(X,nBas,mOrb,nSym,label='X',Ref=XMO)
call Allocate_DT(S,nBas,nBas,nSym,label='S',Ref=Smat)

iOff = 0
lOff = 0
do iSym=1,nSym

  nBa = 0
  do ia=1,nBas(iSym)
    ja = iOff+ia
    tmp = BName(ja)(1:LenIn)
    do j=1,nActa
      if (NamAct(j) == tmp) then
        nBa = nBa+1
        iD(nBa) = ia
      end if
    end do
  end do
  SQ(1:nBas(iSym),1:nBas(iSym)) => SQt(1:nBas(iSym)**2)
  C2(1:nBa,1:mOrb(iSym)) => Scr1(1:nBa*mOrb(iSym),1)
  CC(1:mOrb(iSym),1:nBas(iSym)) => Scr1(1:mOrb(iSym)*nBas(iSym),2)
  X2(1:nBas(iSym),1:mOrb(iSym)) => Scr1(1:nBas(iSym)*mOrb(iSym),3)
  Z(1:nBa,1:mOrb(iSym)) => Scr1(1:nBa*mOrb(iSym),4)
  Scr(1:nBas(iSym),1:mOrb(iSym)) => Scr1(1:nBas(iSym)*mOrb(iSym),5)
  do ia=1,nBa
    C2(ia,:) = X%SB(iSym)%A2(iD(ia),:)
  end do
  do ia=1,nBa
    SQ(:,ia) = S%SB(iSym)%A2(:,iD(ia))
  end do
  nBx = max(1,nBas(iSym))
  nBax = max(1,nBa)
  call DGEMM_('T','N',nBa,mOrb(iSym),nBas(iSym),One,SQ,nBx,X%SB(iSym)%A2,nBx,Zero,Z,nBax)
  do i=1,mOrb(iSym)
    Q(i) = ddot_(nBa,C2(:,i),1,Z(:,i),1)**2
  end do
  Z(1:nBas(iSym),1:mOrb(iSym)) => Scr1(1:nBas(iSym)*mOrb(iSym),4)
  n_OK(iSym) = 0
  n_KO = 0
  do i=1,mOrb(iSym)
    ThrS = ThrSel*Saa(lOff+i)
    if (sqrt(Q(i)) >= ThrS) then
      n_OK(iSym) = n_OK(iSym)+1
      X2(:,n_OK(iSym)) = X%SB(iSym)%A2(:,i)
      iD(nBmx+n_OK(iSym)) = i
    else
      n_KO = n_KO+1
      Z(:,n_KO) = X%SB(iSym)%A2(:,i)
      iD(nBmx+nOrbmx+n_KO) = i
    end if
  end do

  mOx = max(1,mOrb(iSym))

  if (.not. ortho) then

    call Ortho_orb(X2,S%SB(iSym)%A2,nBas(iSym),n_OK(iSym),2,.false.)
    call Ortho_orb(Z,S%SB(iSym)%A2,nBas(iSym),n_KO,2,.false.)
  end if

  call DGEMM_('T','N',mOrb(iSym),nBas(iSym),nBas(iSym),One,C%SB(iSym)%A2,nBx,S%SB(iSym)%A2,nBx,Zero,CC,mOx)

  if (n_KO > 0) then
    call DGEMM_('N','N',mOrb(iSym),n_KO,nBas(iSym),One,CC,mOx,Z,nBx,Zero,U,mOx)

    call Get_Can_Lorb(Eorb(lOff+1),Fock(:,1),n_KO,mOrb(iSym),iD(nBmx+nOrbmx+1:),U)

    call DGEMM_('N','N',nBas(iSym),n_KO,n_KO,One,Z,nBx,U,n_KO,Zero,Scr,nBx)

    ! Reorder the final MOs such that those of the active site come first
    km = n_OK(iSym)+1
    C%SB(iSym)%A2(:,km:km+n_KO-1) = Scr(:,1:n_KO)
    Fock(:,2) = Fock(:,1)
  end if

  call DGEMM_('N','N',mOrb(iSym),n_OK(iSym),nBas(iSym),One,CC,mOx,X2,nBx,Zero,U,mOx)

  call Get_Can_Lorb(Eorb(lOff+1),Fock(:,1),n_OK(iSym),mOrb(iSym),iD(nBmx+1:),U)

  nOx = max(1,n_OK(iSym))
  call DGEMM_('N','N',nBas(iSym),n_OK(iSym),n_OK(iSym),One,X2,nBx,U,nOx,Zero,Scr,nBx)

  do i=1,n_OK(iSym)
    Eorb(lOff+i) = Fock(iD(nBmx+i),1)
  end do
  C%SB(iSym)%A2(:,1:n_OK(iSym)) = Scr(:,1:n_OK(iSym))

  do i=1,n_KO
    Eorb(lOff+n_OK(iSym)+i) = Fock(iD(nBmx+nOrbmx+i),2)
  end do

  iOff = iOff+nBas(iSym)
  lOff = lOff+mOrb(iSym)
end do

nullify(SQ)
nullify(C2)
nullify(CC)
nullify(X2)
nullify(Z)
nullify(Scr)

call mma_deallocate(iD)
call mma_deallocate(SQt)
call mma_deallocate(Scr1)
call mma_deallocate(U)
call mma_deallocate(Q)
call mma_deallocate(Fock)
call Deallocate_DT(C)
call Deallocate_DT(X)
call Deallocate_DT(S)

return

end subroutine Get_Orb_Select
