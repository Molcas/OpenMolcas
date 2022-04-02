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
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iQ_Atoms, iCNum
real(kind=wp), intent(in) :: Coord(3,iQ_Atoms)
integer(kind=iwp) :: i, ind, iTemp, j, k
real(kind=wp) :: Atemp, S, Sbig
logical(kind=iwp) :: Changed
character(len=200) :: Head
integer(kind=iwp), allocatable :: IndexSet(:)
real(kind=wp), allocatable :: AvstPart(:), CordstTemp(:,:)

call mma_allocate(AvstPart,nPart,label='AvstPart')

!For each solvent particle, compute the smallest distance to any QM-atom from the oxygen of water.
do i=1,nPart
  k = (i-1)*nCent+1
  Sbig = 1.0e20_wp
  do j=1,iQ_Atoms
    S = (Coord(1,j)-Cordst(1,k))**2+(Coord(2,j)-Cordst(2,k))**2+(Coord(3,j)-Cordst(3,k))**2
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
CordstTemp(:,:) = Cordst(:,1:nPart*nCent)
do i=1,nPart
  k = (i-1)*nCent
  ind = (IndexSet(i)-1)*nCent
  Cordst(:,k+1:k+nCent) = CordstTemp(:,ind+1:ind+nCent)
end do
call mma_deallocate(IndexSet)
call mma_deallocate(CordstTemp)

! Substitute the first coordinate slots with QM-molecule, or since we have
! ordered above, this is equivalent with removing closest solvents and there put QM-mol.
Cordst(:,1:iQ_Atoms) = Coord
! Just dummy-coordinates added to empty slots.
Cordst(1,iQ_Atoms+1:iCnum*nCent) = Coord(1,1)
Cordst(2,iQ_Atoms+1:iCnum*nCent) = Coord(2,1)
Cordst(3,iQ_Atoms+1:iCnum*nCent) = Coord(3,1)

if (iPrint >= 10) then !Optional printing.
  write(Head,*) 'Coordinates of the system after substitution and reordering of solvent molecules.'
  call Cooout(Head,Cordst,nPart,nCent)
end if

return

end subroutine PlaceIt
