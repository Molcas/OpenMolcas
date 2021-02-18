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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine perform some initialization after reading the first     *
! one electron file.                                                   *
! These parameters are checked when reading subsequent one electron    *
! files.                                                               *
!                                                                      *
!======================================================================*
!                                                                      *
! Author: Per-Olof Widmark                                             *
!         IBM Sweden                                                   *
!                                                                      *
!***********************************************************************

subroutine Init_GenANO()

use Genano_globals, only: MxLqn, nSym, nBas, nPrim, nDsym, iSymBk, pDsym, tDsym, Ssym, LenIn, Center, BasName, symlab
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: i, iBas, iComp, ind, iShell, iSym, l, length, nAtoms, next
character(len=LenIn+4), allocatable :: LblCnt(:)
logical(kind=iwp) :: Found

call Get_nAtoms_all(nAtoms)
call mma_allocate(LblCnt,nAtoms,label='LblCnt')
call Get_cArray('LP_L',LblCnt,(LenIn+4)*nAtoms)
Found = .false.
do i=1,nAtoms
  if (LblCnt(i)(1:len(Center)) == Center) Found = .true.
end do
call mma_deallocate(LblCnt)
if (.not. Found) then
  call WarningMessage(2,'Center '//Center//' not found')
  call Quit_OnUserError()
end if
ind = 0
do iSym=1,nSym
  !write(u6,'(a,i1,a)') ' <<< Symmetry ',iSym,' >>>'
  do iBas=1,nBas(iSym)
    ind = ind+1
    !write(u6,*) BasName(ind)
    do l=0,MxLqn
      if (BasName(ind)(1:len(Center)) == Center) then
        if (BasName(ind)(len(Center)+1:) == symlab(l*(l+1)+1)) then
          nPrim(l) = nPrim(l)+1
        end if
      end if
    end do
  end do
end do
write(u6,*)
write(u6,'(a,8i5)') 'Number of primitives per shell:',nPrim
nDsym = 0
do l=0,MxLqn
  !write(u6,'(1x,a,i5)') symlab(l*(l+1)+1),nPrim(l)
  nDsym = nDsym+(2*l+1)*nPrim(l)*(nPrim(l)+1)/2
end do
call mma_allocate(pDsym,nDsym,label='pDsym')
call mma_allocate(tDsym,nDsym,label='tDsym')
call mma_allocate(Ssym,nDsym,label='Ssym')
pDsym(:) = Zero
tDsym(:) = Zero
!-----------------------------------------------------------------------
ind = 0
length = 0
next = 1
do iShell=0,MxLqn
  length = nPrim(iShell)*(nPrim(iShell)+1)/2
  do iComp=-iShell,iShell
    ind = ind+1
    iSymBk(ind) = next
    next = next+length
  end do
end do
!write(u6,*) 'In Init: density block pointers'
!write(u6,'(1x,10i5)') (iSymBk(i),i=1,(MxLqn+1)**2

return

end subroutine Init_GenANO
