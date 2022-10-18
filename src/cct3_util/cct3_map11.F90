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

subroutine cct3_map11(a,b,dimp,nfact)
! mapping A(p) -> nfact*B(p)

integer dimp, nfact
real*8 a(1:dimp)
real*8 b(1:dimp)
integer pp

if (nfact == 1) then

  do pp=1,dimp
    b(pp) = a(pp)
  end do

else

  do pp=1,dimp
    b(pp) = -a(pp)
  end do

end if

return

end subroutine cct3_map11
