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
! Copyright (C) 1990,1992, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

#ifdef _OLD_CODE_

subroutine NTMul(A,B,C,nRowA,nColA,nRowB)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nRowA, nColA, nRowB
real(kind=wp), intent(in) :: A(nRowA,nColA), B(nRowB,nColA)
real(kind=wp), intent(out) :: C(nRowA,nRowB)
integer(kind=iwp) :: i, Indj, j, jj, k, mCache, nCache, njVec

nCache = (64/8)*1024
mCache = (nCache*3)/4-nRowA*nColA
Incj = mCache/(nRowA+nColA)

! Sectioning of long index

do jj=1,nRowB,Incj
  njVec = min(Incj,nRowB-jj+1)

  do i=1,nRowA
    ! Set target to zero
    C(i,jj:jj+njVec-1) = Zero
    do k=1,nColA
      if (A(i,k) /= Zero) C(i,jj:jj+njVec-1) = C(i,jj:jj+njVec-1)+A(i,k)*B(jj:jj+njVec-1,k)
    end do
  end do

end do    ! End of sectioning

return

end subroutine NTMul

#else

subroutine ntmul(a,b,r,ncol,nlink,nrow)

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nCol, nLink, nRow
real(kind=wp), intent(in) :: a(ncol,*), b(nrow,*)
real(kind=wp), intent(out) :: r(ncol,nrow)
integer(kind=iwp), parameter :: mxind = 2000
integer(kind=iwp) :: i, ind(mxind), j, k, nnot, nr1
real(kind=wp) :: S(16)

do i=1,ncol

  nnot = 0
  do k=1,min(nlink,mxind)
    if (a(i,k) /= Zero) then
      nnot = nnot+1
      ind(nnot) = k
    end if
  end do

  do j=1,nrow-15,16
    s(:) = Zero
    do k=1,nnot
      s(:) = s(:)+a(i,ind(k))*b(j:j+15,ind(k))
    end do
    r(i,j:j+15) = s(:)
  end do

  nr1 = mod(nrow,16)
  if (nr1 == 0) cycle
  j = nrow-nr1+1

  if (nr1 >= 8) then
    s(1:8) = Zero
    do k=1,nnot
      s(1:8) = s(1:8)+a(i,ind(k))*b(j:j+7,ind(k))
    end do
    r(i,j:j+7) = s(1:8)
    nr1 = nr1-8
    j = j+8
  end if

  if (nr1 >= 4) then
    s(1:4) = Zero
    do k=1,nnot
      s(1:4) = s(1:4)+a(i,ind(k))*b(j:j+3,ind(k))
    end do
    r(i,j:j+3) = s(1:4)
    nr1 = nr1-4
    j = j+4
  end if

  if (nr1 == 3) then
    s(1:3) = Zero
    do k=1,nnot
      s(1:3) = s(1:3)+a(i,ind(k))*b(j:j+2,ind(k))
    end do
    r(i,j:j+2) = s(1:3)
  else if (nr1 == 2) then
    s(1:2) = Zero
    do k=1,nnot
      s(1:2) = s(1:2)+a(i,ind(k))*b(j:j+1,ind(k))
    end do
    r(i,j:j+1) = s(1:2)
  else if (nr1 == 1) then
    s(1) = Zero
    do k=1,nnot
      s(1) = s(1)+a(i,ind(k))*b(j,ind(k))
    end do
    r(i,j) = s(1)
  else if (nr1 > 3) then
    call WarningMessage(2,'nr1 > 3')
    call Abend()
  end if
end do

if (mxind >= nlink) return
call WarningMessage(2,'MxInd < nLink')
write(u6,*) 'mxind,nlink=',mxind,nlink
call Abend()

return

end subroutine ntmul

#endif
