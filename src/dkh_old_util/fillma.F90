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

subroutine FILLMA(N,S,OVE)

implicit real*8(A-H,O-Z)
dimension S(*), OVE(N,N)

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    OVE(I,J) = S(IJ)
    OVE(J,I) = S(IJ)
  end do
end do

return

end subroutine FILLMA
