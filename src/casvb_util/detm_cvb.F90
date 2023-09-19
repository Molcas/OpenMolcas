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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

function detm_cvb(a,n)

implicit real*8(a-h,o-z)
#include "WrkSpc.fh"
dimension a(n*n)
dimension det(2)
!start linpack_determinant
save zero, one
data zero/0d0/,one/1d0/

if (n == 0) then
  detm_cvb = one
  return
end if
i1 = mstackr_cvb(n*n)
i2 = mstacki_cvb(n)
ierr = 0
call fmove_cvb(a,work(i1),n*n)
call dgetrf_(n,n,work(i1),n,iwork(i2),ierr)
!start linpack_determinant
!call dgefa(work(i1),n,n,iwork(i2),ierr)
i3 = mstackr_cvb(n*n)
if (ierr /= 0) then
  detm_cvb = zero
  call mfreer_cvb(i1)
  return
end if
call dgedi(work(i1),n,n,iwork(i2),det,work(i3),10)
detm_cvb = det(1)*10d0**det(2)
call mfreer_cvb(i1)

return

end function detm_cvb
