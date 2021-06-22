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

subroutine orb2tpidx_sym(NF,NI,N1,N2,N3,NS,ND,TYPEINDEX)
!SVC: convert orbital partition info to a typeindex array
!     corresponding to a _specific_symmetry_ (so these variables are
!     scalars!!
!     A typeindex array consists of integers with 7 possible values
!     corresponding to the types 'fi123sd' -> '1234567'.

use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NF, NI, N1, N2, N3, NS, ND
integer(kind=iwp), intent(_OUT_) :: typeindex(*)
integer(kind=iwp) :: IB, iOff

iOff = 0

do IB=1,NF
  TYPEINDEX(iOff+IB) = 1
end do
iOff = iOff+NF

do IB=1,NI
  TYPEINDEX(iOff+IB) = 2
end do
iOff = iOff+NI

do IB=1,N1
  TYPEINDEX(iOff+IB) = 3
end do
iOff = iOff+N1

do IB=1,N2
  TYPEINDEX(iOff+IB) = 4
end do
iOff = iOff+N2

do IB=1,N3
  TYPEINDEX(iOff+IB) = 5
end do
iOff = iOff+N3

do IB=1,NS
  TYPEINDEX(iOff+IB) = 6
end do
iOff = iOff+NS

do IB=1,ND
  TYPEINDEX(iOff+IB) = 7
end do

end subroutine orb2tpidx_sym
