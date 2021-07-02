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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************

subroutine Local_Triang(nRow,A)
! This routine is a modification of the Per-AAke's Triang routine
! found in src/caspt2/triang.f
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

implicit real*8(a-h,o-z)
implicit integer(i-n)
dimension A(nRow**2)

! Convert a square matrix to triangular in-place.
iFrom = 1+nRow
iTo = 2
do i=2,nRow
  call dCopy_(i,A(iFrom),1,A(iTo),1)
  iFrom = iFrom+nRow
  iTo = iTo+i
end do

return

end subroutine Local_Triang
