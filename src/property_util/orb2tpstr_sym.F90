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

subroutine orb2tpstr_sym(NF,NI,N1,N2,N3,NS,ND,TYPESTRING)
!SVC: convert orbital partition info to a typestring array
!     corresponding to a _specific_symmetry_ (so these variables are
!     scalars!!
!     A typestring array consists of characters of 'fi123sd'

use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NF, NI, N1, N2, N3, NS, ND
character, intent(_OUT_) :: typestring(*)
integer(kind=iwp) :: IB, iOff

iOff = 0

do IB=1,NF
  TYPESTRING(iOff+IB) = 'F'
end do
iOff = iOff+NF

do IB=1,NI
  TYPESTRING(iOff+IB) = 'I'
end do
iOff = iOff+NI

do IB=1,N1
  TYPESTRING(iOff+IB) = '1'
end do
iOff = iOff+N1

do IB=1,N2
  TYPESTRING(iOff+IB) = '2'
end do
iOff = iOff+N2

do IB=1,N3
  TYPESTRING(iOff+IB) = '3'
end do
iOff = iOff+N3

do IB=1,NS
  TYPESTRING(iOff+IB) = 'S'
end do
iOff = iOff+NS

do IB=1,ND
  TYPESTRING(iOff+IB) = 'D'
end do

end subroutine orb2tpstr_sym
