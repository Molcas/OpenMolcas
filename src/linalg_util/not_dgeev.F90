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

subroutine not_DGeEV(iOpt,a,lda,w,z,ldz,n,aux,naux)

use stdalloc, only: mma_allocate, mma_deallocate

implicit real*8(a-h,o-z)
real*8 a(lda,n), w(2,n), z(2*ldz*n), aux(naux)
real*8, allocatable :: w1(:)

if (iOpt == 2) then
  write(6,*) 'not_DGeEV: iOpt=2 is not implemented yet!'
  call Abend()
end if
if (ldz /= n) then
  write(6,*) 'not_DGeEV: ldz=/=n is not implemented yet!'
  call Abend()
end if
if (iOpt == 0) then
  write(6,*) 'not_DGeEV: iOpt=0 is not implemented yet!'
  call Abend()
end if
iOff = n+1
if (nAux < 2*n) then
  write(6,*) 'not_DGeEV: nAux is too small (naux<2*n)!'
  call Abend()
end if
call mma_allocate(w1,n,label='w1')
iErr = 0
call XEIGEN(iOpt,lda,n,a,w,w1,z,iErr)
if (iErr /= 0) then
  write(6,*) ' not_DGeEV: iErr=/= 0!'
  call Abend()
end if
!call RecPrt('w',' ',w,n,1)
!call RecPrt('w1',' ',w1,n,1)
!call RecPrt('z',' ',z,n,n)

!-----Order eigenvalues to ESSL standard

call dcopy_(n,w,1,aux,1)
do i=1,n
  w(1,i) = aux(i)    ! Real
  w(2,i) = w1(i)     ! Imaginary
end do
call mma_deallocate(w1)

!-----Order eigenvector to ESSL standard

i = N
do while (i >= 1)
  if (w(2,i) /= 0.0d0) then
    i = i-1
    iOff = (i-1)*n
    call dcopy_(2*N,Z(iOff+1),1,Aux,1)
    iOff = (i-1)*2*n
    call dcopy_(N,Aux(1),1,Z(iOff+1),2)
    call dcopy_(N,Aux(1+N),1,Z(iOff+2),2)
    iOff = i*2*n
    call dcopy_(N,Aux(1),1,Z(iOff+1),2)
    call dcopy_(N,Aux(1+N),1,Z(iOff+2),2)
    call DScal_(N,-1.0d0,Z(iOff+2),2)
  else
    iOff = (i-1)*n
    call dcopy_(N,Z(iOff+1),1,Aux,1)
    iOff = (i-1)*2*n
    call dcopy_(N,Aux(1),1,Z(iOff+1),2)
    call dcopy_(N,[0.0d0],0,Z(iOff+2),2)
  end if
  i = i-1
end do

return

end subroutine not_DGeEV
