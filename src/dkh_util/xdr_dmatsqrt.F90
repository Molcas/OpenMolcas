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

subroutine XDR_dmatsqrt(a,n)
! Compute the inverse square root of a real symmetric matrix : A**(-1/2)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(inout) :: a(n*n)
integer(kind=iwp) :: iCr, iW, iTmp, i, j, k, info
real(kind=wp) :: dia
#include "WrkSpc.fh"

call getmem('tmp ','ALLOC','REAL',iTmp,8*n)
call getmem('Cr  ','ALLOC','REAL',iCr,n*n+4)
call getmem('Eig ','ALLOC','REAL',iW,n+4)
do i=0,n*n-1
  Work(iCr+i) = a(i+1)
end do
call dsyev_('V','L',n,Work(iCr),n,Work(iW),Work(iTmp),8*n,info)
do i=0,n-1
  dia = One/sqrt(sqrt(Work(iW+i)))
  do j=0,n-1
    k = i*n+j
    Work(iCr+k) = Work(iCr+k)*dia
  end do
end do
call dmxma(n,'N','T',Work(iCr),Work(iCr),a,One)
call getmem('tmp ','FREE','REAL',iTmp,8*n)
call getmem('Cr  ','FREE','REAL',iCr,n*n+4)
call getmem('Eig ','FREE','REAL',iW,n+4)

return

end subroutine XDR_dmatsqrt
