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

subroutine Ext_W4hlp1(V2,M1,nc,dimab,dimapp,dimabpp,addapp)
! this routine does:
! Extract M1(m,a"b") <- V2(m,a'b')

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nc, dimab, dimapp, dimabpp, addapp
real(kind=wp) :: V2(nc,dimab), M1(nc,dimabpp)
integer(kind=iwp) :: a, ab, abpp, app, bpp, m

abpp = 0
do app=1,dimapp
  a = addapp+app
  ab = a*(a-1)/2+addapp
  do bpp=1,app
    ab = ab+1
    abpp = abpp+1

    do m=1,nc
      M1(m,abpp) = V2(m,ab)
    end do

  end do
end do

return

end subroutine Ext_W4hlp1
