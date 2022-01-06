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

subroutine CheckDomain(irc,iDomain,nAtom,nOcc)
! Thomas Bondo Pedersen, January 2006.
!
! Purpose: check domain definition.

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: nAtom, nOcc, iDomain(0:nAtom,nOcc)
integer(kind=iwp) :: i, iAt, iAtom

irc = 0
do i=1,nOcc
  if ((iDomain(0,i) < 1) .or. (iDomain(0,i) > nAtom)) then
    write(u6,*) 'Dimension of domain ',i,': ',iDomain(0,i)
    irc = irc+1
  else
    do iAt=1,iDomain(0,i)
      iAtom = iDomain(iAt,i)
      if ((iAtom < 1) .or. (iAtom > nAtom)) then
        write(u6,*) 'Atom ',iAt,' of domain ',i,': ',iAtom
        irc = irc+1
      end if
    end do
  end if
end do

end subroutine CheckDomain
