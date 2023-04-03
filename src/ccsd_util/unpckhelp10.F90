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

subroutine unpckhelp10(a,b,dimp,dimq,dime,dimf,eadd,noe,fadd,nof,bb,dimb)
! this routine does:
! b(e,f,_Bb) = a(pe,qf)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dimp, dimq, dime, dimf, eadd, noe, fadd, nof, bb, dimb
real(kind=wp) :: a(dimp,dimq), b(dime,dimf,dimb)
integer(kind=iwp) :: f, pe, qf

do qf=fadd+1,fadd+nof
  f = qf-fadd
  do pe=eadd+1,eadd+noe
    b(pe-eadd,f,bb) = a(pe,qf)
  end do
end do

return

end subroutine unpckhelp10
