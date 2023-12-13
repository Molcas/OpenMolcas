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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Sort_Localisation(CMO,nBas,nOcc,nFro,nSym)
! Thomas Bondo Pedersen, November 2005.
!
! Purpose: sort CMOs according to Cholesky orbital ordering.

use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use OneDat, only: sNoNuc, sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOcc(nSym), nFro(nSym)
real(kind=wp), intent(inout) :: CMO(*)
integer(kind=iwp) :: iComp, iOpt, irc, iSyLbl, iSym, kC
real(kind=wp) :: xNrm
character(len=80) :: Txt
character(len=8) :: Label
type(DSBA_Type) :: C, Oaux, Ovlp, X
real(kind=wp), allocatable :: Den(:,:), Scr(:,:), U(:,:)
real(kind=wp), parameter :: ThrCho = 1.0e-12_wp ! Static setting of decomposition threshold
character(len=*), parameter :: SecNam = 'Sort_Localisation'

! Get a copy of the occupied orbitals: X=CMO.
! -------------------------------------------

call Allocate_DT(X,nBas,nOcc,nSym,label='XCho')
call Allocate_DT(C,nBas,nBas,nSym,label='C',Ref=CMO)
do iSym=1,nSym
  kC = nFro(iSym)+1
  X%SB(iSym)%A2(:,:) = C%SB(iSym)%A2(:,kC:kC+nOcc(iSym)-1)
end do

! Get the overlap matrix.
! -----------------------

call Allocate_DT(Ovlp,nBas,nBas,nSym,label='Ovlp')
call Allocate_DT(Oaux,nBas,nBas,nSym,aCase='TRI',label='Ovlp')

irc = -1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
iComp = 1
iSyLbl = 1
Label = 'Mltpl  0'
call RdOne(irc,iOpt,Label,iComp,Oaux%A0,iSyLbl)
if (irc /= 0) then
  write(u6,*) SecNam,': RdOne returned ',irc
  write(u6,*) 'Label = ',Label,'  iSyLbl = ',iSyLbl
  call SysAbendMsg(SecNam,'I/O error in RdOne',' ')
end if

do iSym=1,nSym
  call Tri2Rec(Oaux%SB(iSym)%A1,Ovlp%SB(iSym)%A2,nBas(iSym),.false.)
end do
call Deallocate_DT(Oaux)

! Sort each symmetry block.
! -------------------------

do iSym=1,nSym

  ! Cycle loop for empty symmetry blocks.
  ! -------------------------------------

  if ((nBas(iSym) < 1) .or. (nOcc(iSym) < 1)) cycle

  ! Allocations.
  ! ------------

  call mma_allocate(Den,nBas(iSym),nBas(iSym),label='SrtDen')
  call mma_allocate(U,nOcc(iSym),nOcc(iSym),label='SrtU')
  call mma_allocate(Scr,nBas(iSym),nOcc(iSym),label='SrtScr')

  ! Cholesky decompose D=C^TC and thus define the ordering.
  ! At this stage, X contains the original MOs (CMO).
  ! After the decomposition, X contains the Cholesky MOs.
  ! -------------------------------------------------------

  call GetDens_Localisation(Den,X%SB(iSym)%A2,nBas(iSym),nOcc(iSym))
  irc = -1
  call ChoLoc(irc,Den,X%SB(iSym)%A2,ThrCho,xNrm,nBas(iSym),nOcc(iSym))
  if (irc /= 0) then
    write(u6,*) SecNam,': ChoLoc returned ',irc
    write(u6,*) 'Symmetry block: ',iSym
    write(u6,*) 'Unable to continue...'
    write(Txt,'(A,I6)') 'ChoLoc return code:',irc
    call SysAbendMsg(SecNam,'Density Cholesky decomposition failed!',Txt)
  end if

  ! Compute address in CMO, skipping frozen orbitals.
  ! -------------------------------------------------

  kC = nFro(iSym)+1

  ! Compute U = X^TSC.
  ! ------------------

  call GetUmat_Localisation(U,X%SB(iSym)%A2,Ovlp%SB(iSym)%A2,C%SB(iSym)%A2(:,kC:),Scr,nBas(iSym),nOcc(iSym))

  ! Sort.
  ! -----

  call Sort_Localisation_1(C%SB(iSym)%A2(:,kC:),U,nBas(iSym),nOcc(iSym))

  ! De-allocations.
  ! ---------------

  call mma_deallocate(Den)
  call mma_deallocate(U)
  call mma_deallocate(Scr)

end do

! De-allocations.
! ---------------

call Deallocate_DT(X)
call Deallocate_DT(C)
call Deallocate_DT(Ovlp)

end subroutine Sort_Localisation
