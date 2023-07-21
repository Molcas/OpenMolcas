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
! Copyright (C) 2007, Ten-no Research Group                            *
!               2012, Thomas Bondo Pedersen                            *
!***********************************************************************

function FindMx(n,X)
!-----------------------------------------------------------------------
! Function : Find Maximumnorm in X(IDim)
!-----------------------------------------------------------------------

implicit real*8(A-H,O-Z)
real*8 X(n)
XMax = 0.0D+00
do I=1,n
  XMax = max(XMax,abs(X(I)))
end do
FindMx = XMax

return

end function FindMx
