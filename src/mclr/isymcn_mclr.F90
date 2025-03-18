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

integer function ISYMCN_MCLR(ICL,IOP,NCL,NOPEN)
! Master routine for symmetry of configuration
! with NCL doubly occupied orbitals and NOPEN singly occupied shells

use MCLR_Data, only: ISMFTO

implicit none
integer ICL(*), IOP(*)
integer NCL, NOPEN
integer IEL, IVV, JVV, KVV

ISYMCN_MCLR = 1
do IEL=1,NOPEN
  IVV = ISYMCN_MCLR-1
  JVV = ISMFTO(IOP(IEL))-1
  KVV = ieor(IVV,JVV)
  ISYMCN_MCLR = KVV+1
end do

! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(ICL)
  call Unused_integer(NCL)
end if

end function ISYMCN_MCLR
