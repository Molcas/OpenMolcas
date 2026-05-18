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
! Copyright (C) 2006, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 2006  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine MKDREF_RPT2(N,G1,DREF,nDREF)

use definitions, only: iwp, wp

implicit none
integer(kind=iwp), intent(in) :: N, NDREF
real(kind=wp), intent(in) :: G1(N,N)
real(kind=wp), intent(out) :: DREF(NDREF)
integer(kind=iwp) I, J, IJ

! Compute DREF(PQ) = <0| Epq |0>
! from G1(P,Q) = <0| Epq |0>
! Storage differs: DREF is triangular.

do I=1,N
  do J=1,I
    IJ = (I*(I-1))/2+J
    DREF(IJ) = G1(I,J)
  end do
end do

end subroutine MKDREF_RPT2
