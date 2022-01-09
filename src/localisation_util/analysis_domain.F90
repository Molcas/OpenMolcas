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
! Copyright (C) 2006, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Analysis_Domain(iDomain,QD,f,Coord,AtomLbl,nBas_Start,nAtom,nBas,nOcc)
! Thomas Bondo Pedersen, January 2006.
!
! Purpose: analyze orbital domains.

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtom, nOcc, iDomain(0:nAtom,nOcc), nBas_Start(nAtom), nBas
real(kind=wp), intent(in) :: QD(nOcc), f(nOcc), Coord(3,nAtom)
character(len=4), intent(in) :: AtomLbl(2,nBas)
integer(kind=iwp) :: i, iAt, iAtom, jAt, jAtom, nAt, nij
real(kind=wp) :: R, Rave, Rmax, Rmin

if ((nAtom < 1) .or. (nOcc < 1)) return

call Cho_Head('Orbital domain analysis','=',80,u6)
do i=1,nOcc
  nAt = iDomain(0,i)
  Rmin = huge(Rmin)
  Rmax = -huge(Rmax)
  Rave = Zero
  nij = 0
  do jAt=1,nAt-1
    jAtom = iDomain(jAt,i)
    do iAt=jAt+1,nAt
      iAtom = iDomain(iAt,i)
      R = sqrt((Coord(1,iAtom)-Coord(1,jAtom))**2+(Coord(2,iAtom)-Coord(2,jAtom))**2+(Coord(3,iAtom)-Coord(3,jAtom))**2)
      Rmin = min(Rmin,R)
      Rmax = max(Rmax,R)
      Rave = Rave+R
      nij = nij+1
    end do
  end do
  if (nij == 0) then
    Rmax = Zero
    Rmin = Zero
  else
    Rave = Rave/real(nij,kind=wp)
  end if
  write(u6,'(/,A,I6,A,I6)') 'Orbital domain',i,':  size:',nAt
  write(u6,'(A,1P,2(1X,D15.5))') '  Charge, completeness function:',QD(i),f(i)
  write(u6,'(A,1P,3(1X,D15.5))') '  Rmin, Rmax, Rave             :',Rmin,Rmax,Rave
  do iAt=1,nAt
    iAtom = iDomain(iAt,i)
    write(u6,'(A,I6,2X,A,1X,3(1X,F12.3))') '  Atom:',iAtom,AtomLbl(1,nBas_Start(iAtom)),Coord(:,iAtom)
  end do
end do

end subroutine Analysis_Domain
