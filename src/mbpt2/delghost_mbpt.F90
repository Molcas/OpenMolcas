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

subroutine DelGHOST_MBPT()

use MBPT2_Global, only: CMO, CMO_Internal, DelGhost, EOrb, nBas, nDsto, nnB, Thr_ghs
use cOrbInf, only: nDel, nExt, nFro, nOcc, nOrb, nSym
use Molcas, only: LenIn
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, iLen, irc, iStart, iStart_t, iSym, nUniqAt, nZero(8)
real(kind=wp), allocatable :: CMO_t(:), EOrb_t(:)
character(len=LenIn+8), allocatable :: UBName(:)
logical(kind=iwp), parameter :: Debug = .false.

if (.not. DelGHOST) return

if (Debug) then
  write(u6,'(A,8I5)') 'nSym:',nSym
  write(u6,'(A,8I5)') 'nBas:',(nBas(i),i=1,nSym)
  write(u6,'(A,8I5)') 'nOrb:',(nOrb(i),i=1,nSym)
  write(u6,'(A,8I5)') 'nOcc:',(nOcc(i),i=1,nSym)
  write(u6,'(A,8I5)') 'nFro:',(nFro(i),i=1,nSym)
end if
do iSym=1,nSym
  nDel(iSym) = nBas(iSym)-nOrb(iSym)
  nExt(iSym) = nOrb(iSym)-nOcc(iSym)-nFro(iSym)
  nDsto(iSym) = nDel(iSym)
  nZero(iSym) = 0
end do

call move_alloc(CMO_Internal,CMO_t)
call move_alloc(EOrb,EOrb_t)

call mma_allocate(CMO_Internal,size(CMO_t),label='CMO')
CMO => CMO_Internal
call mma_allocate(EOrb,size(EOrb_t),label='EOrb')

write(u6,'(A)') '-------------------------------------------------------'
write(u6,'(A)') ' GHOST virtual space removal'
write(u6,'(A)') '-------------------------------------------------------'
write(u6,'(A,8I4)')
write(u6,'(A,8I4)') ' Secondary orbitals before selection:',(nExt(i),i=1,nSym)
write(u6,'(A,8I4)') ' Deleted orbitals before selection:  ',(nDel(i),i=1,nSym)

call Get_iScalar('Unique atoms',nUniqAt)
call mma_allocate(UBName,nnB,label='UBName')
call Get_cArray('Unique Basis Names',UBName,(LenIn+8)*nnB)
call Delete_GHOSTS(irc,nSym,nBas,nFro,nOcc,nZero,nExt,nDel,UBName,nUniqAt,thr_ghs,.false.,CMO_t,EOrb_t)
call mma_deallocate(UBName)

if (irc /= 0) then
  write(u6,*) 'Delete_GHOSTS returned rc= ',irc
  call abend()
end if
write(u6,'(A,8I4)')
write(u6,'(A)') '-------------------------------------------------------'
write(u6,'(A,8I4)')
write(u6,'(A,8I4)')

! set MO coefficients of the deleted orbitals to zero
! Observe that these are not included at all in the basis
iStart = 1
iStart_t = 1
do iSym=1,nSym
  iLen = nOrb(iSym)*nBas(iSym)
  CMO(iStart:iStart+iLen-1) = CMO_t(iStart_t:iStart_t+iLen-1)
  iStart = iStart+iLen
  iStart_t = iStart_t+iLen

  iLen = (nBas(iSym)-nOrb(iSym))*nBas(iSym)
  CMO(iStart:iStart+iLen-1) = Zero
  iStart = iStart+iLen
end do
call mma_deallocate(CMO_t)

! set energies of the deleted orbitals to zero

iStart = 1
iStart_t = 1
do iSym=1,nSym
  iLen = nOrb(iSym)
  EOrb(iStart:iStart+iLen-1) = EOrb_t(iStart_t:iStart+iLen-1)
  iStart = iStart+iLen
  iStart_t = iStart_t+iLen

  iLen = nBas(iSym)-nOrb(iSym)
  EOrb(iStart:iStart+iLen-1) = Zero
  iStart = iStart+iLen
end do
call mma_deallocate(EOrb_t)

return

end subroutine DelGHOST_MBPT
