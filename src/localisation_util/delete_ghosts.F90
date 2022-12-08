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
! Copyright (C) 2010, Francesco Aquilante                              *
!***********************************************************************

subroutine Delete_Ghosts(irc,nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,BName,nUniqAt,ThrS,isCASPT2,CMO,EOrb)
!***********************************************************************
!                                                                      *
! Purpose:  Eliminates MOs of ghost atoms from PT2 treatment           *
!                                                                      *
! Author:   F. Aquilante  (Geneva, July 2010)                          *
!                                                                      *
!***********************************************************************

use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use OneDat, only: sNoNuc, sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nFro(nSym), nIsh(nSym), nAsh(nSym), nUniqAt
integer(kind=iwp), intent(inout) :: nSsh(nSym), nDel(nSym)
character(len=LenIn8), intent(in) :: BName(*)
real(kind=wp), intent(in) :: ThrS
logical(kind=iwp), intent(in) :: isCASPT2
real(kind=wp), intent(inout) :: CMO(*), EOrb(*)
integer(kind=iwp) :: i, ia, iAt, iCMO, iComp, ik, iOff, iOpt, iSym, isymlbl, iv, j, jBas, jOff, jZ, kBas, kCMO, lBas, mAsh, n_KO, &
                     n_OK(nSym), nActa, nBa, nBasT, nBax, nBmx, nBx, NCMO, nOkk, nSmx
character(len=LenIn) :: tmp
character(len=8) :: Label
type(DSBA_Type) :: LCMO, S, SQ
integer(kind=iwp), allocatable :: iD(:), nBas_per_Atom(:), nBas_Start(:)
real(kind=wp), allocatable :: Q(:,:), Qa(:), Qt(:)
real(kind=wp), allocatable, target :: Ct(:), St(:), Xt(:), Zt(:)
real(kind=wp), pointer :: C(:,:), S2(:,:), X(:,:), Z(:,:)
character(len=LenIn), allocatable :: NamAct(:)
real(kind=wp), external :: ddot_

irc = 0

!----------------------------------------------------------------------*
!     GET THE TOTAL NUMBER OF BASIS FUNCTIONS, etc. AND CHECK LIMITS   *
!----------------------------------------------------------------------*

nBasT = 0
nBmx = 0
nSmx = 0
mAsh = 0
do i=1,nSym
  nBasT = nBasT+nBas(i)
  nBmx = max(nBmx,nBas(i))
  nSmx = max(nSmx,nSsh(i))
  mAsh = max(mAsh,nIsh(i)+nAsh(i))
end do
if (nBasT > mxBas) then
  write(u6,'(/6X,A)') 'The number of basis functions exceeds the present limit'
  call Abend()
end if

! nUniqAt = # of symm. unique atoms. Initialize NamAct to blanks.
! ---------------------------------------------------------------

if ((nUniqAt < 1) .or. (nUniqAt > MxAtom)) then
  write(u6,'(A,I9)') 'nUniqAt =',nUniqAt
  call Abend()
end if
call mma_allocate(NamAct,nUniqAt,label='NamAct')
NamAct(:) = ' '

! Allocate and get index arrays for basis functions per atom.
! -----------------------------------------------------------

call mma_allocate(nBas_per_Atom,nUniqAt,label='nB_per_Atom')
call mma_allocate(nBas_Start,nUniqAt,label='nB_Start')

!----------------------------------------------------------------------*
!     Read the overlap matrix                                          *
!----------------------------------------------------------------------*
call Allocate_DT(SQ,nBas,nBas,nSym,label='SMAT')
call Allocate_DT(S,nBas,nBas,nSym,aCase='TRI',label='SLT')
isymlbl = 1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
iComp = 1
Label = 'Mltpl  0'
call RdOne(irc,iOpt,Label,iComp,S%A0,isymlbl)
if (irc /= 0) return
do iSym=1,nSym
  call Square(S%SB(iSym)%A1,SQ%SB(iSym)%A2,1,nBas(iSym),nBas(iSym))
end do
call Deallocate_DT(S)

call Allocate_DT(LCMO,nBas,nBas,nSym,label='LCMO')
NCMO = size(LCMO%A0)
LCMO%A0(:) = CMO(1:NCMO)

!----------------------------------------------------------------------*
!     Compute Mulliken atomic charges of each occupied orbital         *
!             on each center to define the Active Site                 *
!----------------------------------------------------------------------*
call mma_allocate(Q,nUniqAt,mAsh,label='Qai')
call mma_allocate(Qa,nUniqAt,label='Qa')
Qa(:) = Zero
call mma_allocate(Zt,nBmx*mAsh,label='Zm')
lBas = 0
do iSym=1,nSym
  nOkk = nIsh(iSym)+nAsh(iSym)
  nBx = max(1,nBas(iSym))
  Z(1:nBas(iSym),1:nOkk) => Zt(1:nBas(iSym)*nOkk)
  call DGEMM_('N','N',nBas(iSym),nOkk,nBas(iSym),One,SQ%SB(iSym)%A2,nBx,LCMO%SB(iSym)%A2(:,nFro(iSym)+1:),nBx,Zero,Z,nBx)
  jBas = lBas+1
  kBas = lBas+nBas(iSym)
  call BasFun_Atom_Sym(nBas_per_Atom,nBas_Start,BName,jBas,kBas,nUniqAt,.false.)
  do ik=1,nOkk
    do iAt=1,nUniqAt
      Q(iAt,ik) = ddot_(nBas_per_Atom(iAt),LCMO%SB(iSym)%A2(nBas_Start(iAt):,nFro(iSym)+ik),1,Z(nBas_Start(iAt):,ik),1)
    end do
  end do
  do iAt=1,nUniqAt
    Qa(iAt) = Qa(iAt)+ddot_(nOkk,Q(iAt,:),1,Q(iAt,:),1)
    if (sqrt(Qa(iAt)) >= ThrS) then
      if (nBas_per_Atom(iAt) > 0) NamAct(iAt) = BName(lBas+nBas_Start(iAt))(1:LenIn)
    end if
  end do
  lBas = lBas+nBas(iSym)
end do
nullify(Z)
call mma_deallocate(Zt)
call mma_deallocate(Q)
call mma_deallocate(Qa)

! We have now completed the definition of the active site
!----------------------------------------------------------------------*
call mma_allocate(iD,nUniqAt,label='ID_A')
nActa = 0
do iAt=1,nUniqAt
  if (NamAct(iAt)(1:4) /= ' ') then
    nActa = nActa+1
    iD(nActa) = iAt
  end if
end do
do iAt=1,nActa
  NamAct(iAt) = NamAct(iD(iAt))
end do
do iAt=nActa+1,nUniqAt
  NamAct(iAt)(1:4) = ' '
end do
write(u6,*)
write(u6,'(A,F6.3)') ' Threshold for atom selection: ',ThrS
write(u6,*)
if (nActa /= 0) then
  write(u6,'(A,I3,A)') ' Selected ',nActa,' atoms: '
  write(u6,*)
  write(u6,*) NamAct(1:nActa)
  write(u6,*)
else
  write(u6,*) ' None of the occupied non-frozen orbitals has been '
  write(u6,*) ' assigned to the Active region of the molecule.    '
  write(u6,*) ' This is presumably NOT what you want !!!          '
  write(u6,*) ' I will Stop here. Bye Bye !! '
  write(u6,*)
  call Abend()
end if

call mma_deallocate(iD)
call mma_deallocate(nBas_per_Atom)
call mma_deallocate(nBas_Start)

!----------------------------------------------------------------------*
!     Virtual orbital selection                                        *
!----------------------------------------------------------------------*
call mma_allocate(iD,nBmx+2*nSmx,label='ID_vir')
call mma_allocate(St,nBmx**2,label='St')
call mma_allocate(Qt,nSmx,label='Qt')
call mma_allocate(Ct,nBmx*nSmx,label='Ct')
call mma_allocate(Zt,nBmx*nSmx,label='Zt')
call mma_allocate(Xt,nBmx*nSmx,label='Xt')

iOff = 0
do iSym=1,nSym

  nBa = 0
  do ia=1,nBas(iSym)
    tmp = BName(iOff+ia)(1:LenIn)
    do j=1,nActa
      if (NamAct(j) == tmp) then
        nBa = nBa+1
        iD(nBa) = ia
      end if
    end do
  end do

  C(1:nBa,1:nSsh(iSym)) => Ct(1:nBa*nSsh(iSym))
  S2(1:nBas(iSym),1:nBa) => St(1:nBas(iSym)*nBa)
  X(1:nBas(iSym),1:nSsh(iSym)) => Xt(1:nBa*nSsh(iSym))
  Z(1:nBas(iSym),1:nSsh(iSym)) => Zt(1:nBa*nSsh(iSym))

  iCMO = nFro(iSym)+nIsh(iSym)+nAsh(iSym)+1
  do ia=1,nBa
    C(ia,:) = LCMO%SB(iSym)%A2(iD(ia),iCMO:)
    S2(:,ia) = SQ%SB(iSym)%A2(:,iD(ia))
  end do

  nBx = max(1,nBas(iSym))
  nBax = max(1,nBa)
  call DGEMM_('T','N',nBa,nSsh(iSym),nBas(iSym),One,S2,nBx,LCMO%SB(iSym)%A2(:,iCMO:),nBx,Zero,Z,nBax)
  do i=1,nSsh(iSym)
    Qt(i) = ddot_(nBa,C(:,i),1,Z(:,i),1)**2
  end do
  n_OK(iSym) = 0
  n_KO = 0
  do i=1,nSsh(iSym)
    if (sqrt(Qt(i)) >= ThrS) then
      n_OK(iSym) = n_OK(iSym)+1
      X(:,n_OK(iSym)) = LCMO%SB(iSym)%A2(:,iCMO+i-1)
      iD(nBmx+n_OK(iSym)) = i
    else
      n_KO = n_KO+1
      Z(:,n_KO) = LCMO%SB(iSym)%A2(:,iCMO+i-1)
      iD(nBmx+nSmx+n_KO) = i
    end if
  end do

  LCMO%SB(iSym)%A2(:,iCMO:iCMO+n_OK(iSym)-1) = X(:,1:n_OK(iSym))
  kCMO = iCMO+n_OK(iSym)
  LCMO%SB(iSym)%A2(:,kCMO:) = Z(:,1:n_KO)
  if (.not. isCASPT2) then
    jZ = 1
    jOff = iOff+nFro(iSym)+nOkk
    do i=nBmx+1,nBmx+n_OK(iSym)
      iv = jOff+iD(i)
      Z(jZ,1) = EOrb(iv)
      jZ = jZ+1
    end do
    do i=nBmx+nSmx+1,nBmx+nSmx+n_KO
      iv = jOff+iD(i)
      Z(jZ,1) = EOrb(iv)
      jZ = jZ+1
    end do
    EOrb(jOff+1:jOff+nSsh(iSym)) = Z(1:nSsh(iSym),1)
  end if

  iOff = iOff+nBas(iSym)
end do
nullify(C)
nullify(S2)
nullify(X)
nullify(Z)
call mma_deallocate(St)
call mma_deallocate(Qt)
call mma_deallocate(Ct)
call mma_deallocate(Zt)
call mma_deallocate(Xt)
call mma_deallocate(iD)
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Update nSsh, nDel for the Active site PT2
nDel(:) = nDel(:)+nSsh(:)-n_OK(:)
nSsh(:) = n_OK(:)

CMO(1:NCMO) = LCMO%A0

call Deallocate_DT(SQ)
call Deallocate_DT(LCMO)
call mma_deallocate(NamAct)

return

end subroutine Delete_Ghosts
