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

subroutine Cho_SetAtomShl(irc,iAtomShl,n)
!
! Purpose: set mapping from shell to atom (i.e., center).

use ChoArr, only: iSOShl
use stdalloc

implicit real*8(a-h,o-z)
#include "Molcas.fh"
integer iAtomShl(n)
#include "cholesky.fh"
#include "choprint.fh"
#include "choorb.fh"
character(len=14), parameter :: SecNam = 'Cho_SetAtomShl'
character(len=LENIN8) AtomLabel(MxBas)
integer, parameter :: Info_Debug = 4
#ifdef _DEBUGPRINT_
logical, parameter :: Debug = .true.
#else
logical, parameter :: Debug = .false.
#endif
integer, allocatable :: nBas_per_Atom(:)
integer, allocatable :: nBas_Start(:)

if (Debug) write(Lupri,*) '>>> Enter ',SecNam

! Check.
! ------

irc = 0
if (nSym /= 1) then ! does not work with symmetry
  irc = 1
  if (Debug) write(Lupri,*) '>>> Exit ',SecNam,' (error exit: symmetry not allowed!)'
  return
end if
if (n < nShell) call Cho_Quit(SecNam//': iAtomShl not allocated correctly!',104)

! Get number of atoms.
! --------------------

!call Get_nAtoms_All(nAtom)
call Get_iScalar('Bfn Atoms',nAtom)

! Get atomic labels and basis function labels.
! --------------------------------------------

call Get_cArray('Unique Basis Names',AtomLabel,LENIN8*nBasT)

! Allocate and get index arrays for indexation of basis functions on
! each atom.
! ------------------------------------------------------------------

call mma_allocate(nBas_per_Atom,nAtom,Label='nBas_per_Atom')
call mma_allocate(nBas_Start,nAtom,Label='nBas_Start')
call BasFun_Atom(nBas_per_Atom,nBas_Start,AtomLabel,nBasT,nAtom,Debug)

! Set shell-to-atom mapping.
! --------------------------

do iAtom=1,nAtom
  i1 = nBas_Start(iAtom)
  i2 = i1+nBas_per_Atom(iAtom)-1
  do i=i1,i2
    iAtomShl(iSOShl(i)) = iAtom
  end do
end do

if ((iPrint >= Info_Debug) .or. Debug) then
  nErr = 0
  write(Lupri,*)
  write(Lupri,*) SecNam,': shell-to-atom mapping:'
  numSh = 7
  nBatch = (nShell-1)/numSh+1
  do iBatch=1,nBatch
    if (iBatch == nBatch) then
      nSh = nShell-numSh*(nBatch-1)
    else
      nSh = numSh
    end if
    iSh0 = numSh*(iBatch-1)
    iSh1 = iSh0+1
    iSh2 = iSh0+nSh
    write(Lupri,'(/,A,7(1X,I9))') 'Shell:',(iSh,iSh=iSh1,iSh2)
    write(Lupri,'(A,7(1X,I9))') 'Atom :',(iAtomShl(iSh),iSh=iSh1,iSh2)
    do iSh=iSh1,iSh2
      if ((iAtomShl(iSh) < 1) .or. (iAtomShl(iSh) > nAtom)) nErr = nErr+1
    end do
  end do
  if (nErr /= 0) call Cho_Quit(SecNam//': shell-to-atom init failed!',104)
end if

! Deallocations.
! --------------

call mma_deallocate(nBas_Start)
call mma_deallocate(nBas_per_Atom)

if (Debug) write(Lupri,*) '>>> Exit ',SecNam

end subroutine Cho_SetAtomShl
