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

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "qminp.fh"
dimension Coord(MxAt*3)
dimension AvstPart(MxPut), IndexSet(MxPut)
dimension CordstTemp(MxPut*MxCen,3)
character Head*200
logical Changed

!For each solvent particle, compute the smallest distance to any QM-atom from the oxygen of water.
do i=1,nPart
  Sbig = 1D+20
  do j=1,iQ_Atoms
    S = 0
    do k=1,3
      S = S+(Coord((j-1)*3+k)-Cordst(nCent*(i-1)+1,k))**2
    end do
    if (S <= Sbig) then
      Sbig = S
      AvstPart(i) = S
    end if
  end do
end do

do i=1,MxPut
  IndexSet(i) = i
end do

31 continue
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
if (Changed) goto 31

! Put coordinates of solvent suchwise that smallest distances goes first.
do i=1,nPart
  do j=1,nCent
    CordstTemp((i-1)*nCent+j,1) = Cordst((i-1)*nCent+j,1)
    CordstTemp((i-1)*nCent+j,2) = Cordst((i-1)*nCent+j,2)
    CordstTemp((i-1)*nCent+j,3) = Cordst((i-1)*nCent+j,3)
  end do
end do
do i=1,nPart
  ind = IndexSet(i)
  do j=1,nCent
    Cordst((i-1)*nCent+j,1) = CordstTemp((ind-1)*nCent+j,1)
    Cordst((i-1)*nCent+j,2) = CordstTemp((ind-1)*nCent+j,2)
    Cordst((i-1)*nCent+j,3) = CordstTemp((ind-1)*nCent+j,3)
  end do
end do

! Substitute the first coordinate slots with QM-molecule, or since we have
! ordered above, this is equivalent with removing closest solvents and there put QM-mol.
do iz=1,iQ_Atoms
  Cordst(iz,1) = Coord((iz-1)*3+1)
  Cordst(iz,2) = Coord((iz-1)*3+2)
  Cordst(iz,3) = Coord((iz-1)*3+3)
end do
! Just dummy-coordinates added to empty slots.
do iextr=iQ_Atoms+1,iCnum*nCent
  Cordst(iextr,1) = Coord(1)
  Cordst(iextr,2) = Coord(2)
  Cordst(iextr,3) = Coord(3)
end do

if (iPrint >= 10) then !Optional printing.
  write(Head,*) 'Coordinates of the system after substitution and reordening of solvent molecules.'
  call Cooout(Head,Cordst,nPart,nCent)
end if

return

end subroutine PlaceIt
