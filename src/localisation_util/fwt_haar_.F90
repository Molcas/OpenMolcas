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
! Copyright (C) 2009, Francesco Aquilante                              *
!***********************************************************************

subroutine FWT_Haar_(n,m,B,X)

implicit real*8(a-h,o-z)
integer n, m
real*8 B(n,*), X(n,*)
parameter(hlf=0.5d0,one=1.0d0,xone=-1.0d0)

fac = sqrt(hlf)
nv = 2**m
mv = n*nv
do j=m,1,-1  ! A[m]:=X
  nv = nv/2
  kB = nv
  call dcopy_(n,X(1,1),1,B(1,kB),1)
  call daxpy_(n,xone,X(1,2),1,B(1,kB),1)
  call dscal_(n,fac,B(1,kB),1)
  call daxpy_(n,one,X(1,2),1,X(1,1),1)
  call dscal_(n,fac,X(1,1),1)
  do k=2,nv
    kB = kB+1
    jB = 2*k
    lB = jB-1
    call dcopy_(n,X(1,lB),1,B(1,kB),1)
    call daxpy_(n,xone,X(1,jB),1,B(1,kB),1)
    call dscal_(n,fac,B(1,kB),1)
    call daxpy_(n,one,X(1,jB),1,X(1,lB),1)
    call dcopy_(n,X(1,lB),1,X(1,k),1)
    call dscal_(n,fac,X(1,k),1)
  end do
end do
call dcopy_(mv-n,B(1,1),1,X(1,2),1) ! A[0] is already in place

return

end subroutine FWT_Haar_
