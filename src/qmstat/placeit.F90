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

subroutine PlaceIt(Coord,iQ_Atoms,iCNum)

use qmstat_global, only: Cordst, iPrint, nCent, nPart
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iQ_Atoms, iCNum
real(kind=wp) :: Coord(3,iQ_Atoms)
integer(kind=iwp) :: i, iextr, ind, iTemp, iz, j, k
real(kind=wp) :: Atemp, S, Sbig
logical(kind=iwp) :: Changed
character(len=200) :: Head
integer(kind=iwp), allocatable :: IndexSet(:)
real(kind=wp), allocatable :: AvstPart(:), CordstTemp(:,:)

call mma_allocate(AvstPart,nPart,label='AvstPart')

!For each solvent particle, compute the smallest distance to any QM-atom from the oxygen of water.
do i=1,nPart
  Sbig = 1.0e20_wp
  do j=1,iQ_Atoms
    S = Zero
    do k=1,3
      S = S+(Coord(k,j)-Cordst(k,nCent*(i-1)+1))**2
    end do
    if (S <= Sbig) then
      Sbig = S
      AvstPart(i) = S
    end if
  end do
end do

call mma_allocate(IndexSet,nPart,label='IndexSet')
do i=1,nPart
  IndexSet(i) = i
end do

do
  ! Order the indices suchwise that smallest distance goes first. The sorting routine is blunt
  ! but at this stage of the execution time is not a problem.
  Changed = .false.
  do i=1,nPart-1
    if (AvstPart(i+1) < AvstPart(i)) then
      Atemp = AvstPart(i)
      AvstPart(i) = AvstPart(i+1)
      AvstPart(i+1) = Atemp
      iTemp = IndexSet(i)
      IndexSet(i) = IndexSet(i+1)
      IndexSet(i+1) = iTemp
      Changed = .true.
    end if
  end do
  if (.not. Changed) exit
end do

call mma_deallocate(AvstPart)

! Put coordinates of solvent suchwise that smallest distances goes first.
call mma_allocate(CordstTemp,3,nPart*nCent,label='CordstTemp')
do i=1,nPart
  do j=1,nCent
    CordstTemp(1,(i-1)*nCent+j) = Cordst(1,(i-1)*nCent+j)
    CordstTemp(2,(i-1)*nCent+j) = Cordst(2,(i-1)*nCent+j)
    CordstTemp(3,(i-1)*nCent+j) = Cordst(3,(i-1)*nCent+j)
  end do
end do
do i=1,nPart
  ind = IndexSet(i)
  do j=1,nCent
    Cordst(1,(i-1)*nCent+j) = CordstTemp(1,(ind-1)*nCent+j)
    Cordst(2,(i-1)*nCent+j) = CordstTemp(2,(ind-1)*nCent+j)
    Cordst(3,(i-1)*nCent+j) = CordstTemp(3,(ind-1)*nCent+j)
  end do
end do
call mma_deallocate(IndexSet)
call mma_deallocate(CordstTemp)

! Substitute the first coordinate slots with QM-molecule, or since we have
! ordered above, this is equivalent with removing closest solvents and there put QM-mol.
do iz=1,iQ_Atoms
  Cordst(1,iz) = Coord(1,iz)
  Cordst(2,iz) = Coord(2,iz)
  Cordst(3,iz) = Coord(3,iz)
end do
! Just dummy-coordinates added to empty slots.
do iextr=iQ_Atoms+1,iCnum*nCent
  Cordst(1,iextr) = Coord(1,1)
  Cordst(2,iextr) = Coord(2,1)
  Cordst(3,iextr) = Coord(3,1)
end do

if (iPrint >= 10) then !Optional printing.
  write(Head,*) 'Coordinates of the system after substitution and reordering of solvent molecules.'
  call Cooout(Head,Cordst,nPart,nCent)
end if

return

end subroutine PlaceIt
