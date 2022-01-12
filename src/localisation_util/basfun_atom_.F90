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
! Copyright (C) Yannick Carissan                                       *
!               Thomas Bondo Pedersen                                  *
!               Francesco Aquilante                                    *
!***********************************************************************

subroutine BasFun_Atom_Sym(nBas_per_Atom,nBas_Start,BName,jBas,nBas,nAtoms,DoPrint)
! Author: Y. Carissan / T. B. Pedersen
!         [adapted to cases with symmetry by F. Aquilante]

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(in) :: nAtoms, jBas, nBas
integer(kind=iwp), intent(out) :: nBas_per_Atom(nAtoms), nBas_Start(nAtoms)
character(len=LenIn8), intent(in) :: BName(nBas)
logical(kind=iwp), intent(in) :: DoPrint
integer(kind=iwp) :: iAt, iCount, jCount, kBas
character(len=LenIn) :: Lbl
character(len=80) :: Txt, Formt
character(len=LenIn), allocatable :: AtName(:)
character(len=*), parameter :: SecNam = 'BasFun_Atom_Sym'

! Counters.
! ---------

! IFG: To count basis functions per atom, we need a list of atom names,
!      since there is no guarantee all atoms will be present in a give irrep
call mma_allocate(AtName,nAtoms,label='AtName')
call Get_cArray('Unique Atom Names',AtName,LenIn*nAtoms)

kBas = jBas
do iAt=1,nAtoms
  nBas_per_Atom(iAt) = 0
  Lbl = AtName(iAt)
  do while ((BName(kBas)(1:LenIn) == Lbl) .and. (kBas <= nBas))
    nBas_per_Atom(iAt) = nBas_per_Atom(iAt)+1
    kBas = kBas+1
  end do
end do
call mma_deallocate(AtName)

! Offsets.
! --------

iCount = 0
do iAt=1,nAtoms
  nBas_Start(iAt) = iCount+1
  iCount = iCount+nBas_per_Atom(iAt)
end do
jCount = iCount+jBas-1
if (jCount /= nBas) then
  write(Txt,'(A,I9,A,I9)') 'jCount =',jCount,'  nBas =',nBas
  call SysAbendMsg(SecNam,'jCount /= nBas',Txt)
end if

! Print.
! ------

if (DoPrint) then
  write(Formt,'(3(a6,i3,a5))') '(/,a6,',nAtoms,'i5,/,','   a6,',nAtoms,'i5,/,','   a6,',nAtoms,'i5)'
  write(u6,Formt) 'Atom  ',(iAt,iAt=1,nAtoms),'Start ',nBas_Start(:),'nBas  ',nBas_per_Atom(:)
end if

end subroutine BasFun_Atom_Sym
