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

implicit none
integer nRowA, nColA, nRowB
real*8 A(nRowA,nColA), B(nRowB,nColA), C(nRowA,nRowB)
integer nCache, mCache, Indj, jj, njVec, i, j, k

nCache = (64/8)*1024
mCache = (nCache*3)/4-nRowA*nColA
Incj = mCache/(nRowA+nColA)

! Sectioning of long index

do jj=1,nRowB,Incj
  njVec = min(Incj,nRowB-jj+1)

  do i=1,nRowA
    ! Set target to zero
    do j=jj,jj+njVec-1
      C(i,j) = Zero
    end do
    do k=1,nColA
      if (A(i,k) /= Zero) then
        do j=jj,jj+njVec-1
          C(i,j) = C(i,j)+A(i,k)*B(j,k)
        end do
      end if
    end do
  end do

end do    ! End of sectioning

return

end subroutine NTMul

#else
subroutine ntmul(a,b,r,ncol,nlink,nrow)

implicit none
integer nCol, nRow, nLink
real*8 r(ncol,*), a(ncol,*), b(nrow,*)
integer, parameter :: mxind = 2000
integer ind(mxind)
integer i, nnot, k, j, nr1
real*8 S1, S2, S3, S4, S5, S6, S7, S8, T1, T2, T3, T4, T5, T6, T7, T8

do i=1,ncol

  nnot = 0
  do k=1,min(nlink,mxind)
    if (a(i,k) /= 0.0d0) then
      nnot = nnot+1
      ind(nnot) = k
    end if
  end do

  do j=1,nrow-15,16
    s1 = 0.0d0
    s2 = 0.0d0
    s3 = 0.0d0
    s4 = 0.0d0
    s5 = 0.0d0
    s6 = 0.0d0
    s7 = 0.0d0
    s8 = 0.0d0
    t1 = 0.0d0
    t2 = 0.0d0
    t3 = 0.0d0
    t4 = 0.0d0
    t5 = 0.0d0
    t6 = 0.0d0
    t7 = 0.0d0
    t8 = 0.0d0
    do k=1,nnot
      s1 = s1+a(i,ind(k))*b(j,ind(k))
      s2 = s2+a(i,ind(k))*b(j+1,ind(k))
      s3 = s3+a(i,ind(k))*b(j+2,ind(k))
      s4 = s4+a(i,ind(k))*b(j+3,ind(k))
      s5 = s5+a(i,ind(k))*b(j+4,ind(k))
      s6 = s6+a(i,ind(k))*b(j+5,ind(k))
      s7 = s7+a(i,ind(k))*b(j+6,ind(k))
      s8 = s8+a(i,ind(k))*b(j+7,ind(k))
      t1 = t1+a(i,ind(k))*b(j+8,ind(k))
      t2 = t2+a(i,ind(k))*b(j+9,ind(k))
      t3 = t3+a(i,ind(k))*b(j+10,ind(k))
      t4 = t4+a(i,ind(k))*b(j+11,ind(k))
      t5 = t5+a(i,ind(k))*b(j+12,ind(k))
      t6 = t6+a(i,ind(k))*b(j+13,ind(k))
      t7 = t7+a(i,ind(k))*b(j+14,ind(k))
      t8 = t8+a(i,ind(k))*b(j+15,ind(k))
    end do
    r(i,j) = s1
    r(i,j+1) = s2
    r(i,j+2) = s3
    r(i,j+3) = s4
    r(i,j+4) = s5
    r(i,j+5) = s6
    r(i,j+6) = s7
    r(i,j+7) = s8
    r(i,j+8) = t1
    r(i,j+9) = t2
    r(i,j+10) = t3
    r(i,j+11) = t4
    r(i,j+12) = t5
    r(i,j+13) = t6
    r(i,j+14) = t7
    r(i,j+15) = t8
  end do

  nr1 = mod(nrow,16)
  if (nr1 == 0) goto 100
  j = nrow-nr1+1

  if (nr1 >= 8) then
    s1 = 0.0d0
    s2 = 0.0d0
    s3 = 0.0d0
    s4 = 0.0d0
    s5 = 0.0d0
    s6 = 0.0d0
    s7 = 0.0d0
    s8 = 0.0d0
    do k=1,nnot
      s1 = s1+a(i,ind(k))*b(j,ind(k))
      s2 = s2+a(i,ind(k))*b(j+1,ind(k))
      s3 = s3+a(i,ind(k))*b(j+2,ind(k))
      s4 = s4+a(i,ind(k))*b(j+3,ind(k))
      s5 = s5+a(i,ind(k))*b(j+4,ind(k))
      s6 = s6+a(i,ind(k))*b(j+5,ind(k))
      s7 = s7+a(i,ind(k))*b(j+6,ind(k))
      s8 = s8+a(i,ind(k))*b(j+7,ind(k))
    end do
    r(i,j) = s1
    r(i,j+1) = s2
    r(i,j+2) = s3
    r(i,j+3) = s4
    r(i,j+4) = s5
    r(i,j+5) = s6
    r(i,j+6) = s7
    r(i,j+7) = s8
    nr1 = nr1-8
    j = j+8
  end if

  if (nr1 >= 4) then
    s1 = 0.0d0
    s2 = 0.0d0
    s3 = 0.0d0
    s4 = 0.0d0
    do k=1,nnot
      s1 = s1+a(i,ind(k))*b(j,ind(k))
      s2 = s2+a(i,ind(k))*b(j+1,ind(k))
      s3 = s3+a(i,ind(k))*b(j+2,ind(k))
      s4 = s4+a(i,ind(k))*b(j+3,ind(k))
    end do
    r(i,j) = s1
    r(i,j+1) = s2
    r(i,j+2) = s3
    r(i,j+3) = s4
    nr1 = nr1-4
    j = j+4
  end if

  if (nr1 == 1) then
    s1 = 0.0d0
    do k=1,nnot
      s1 = s1+a(i,ind(k))*b(j,ind(k))
    end do
    r(i,j) = s1
  else if (nr1 == 2) then
    s1 = 0.0d0
    s2 = 0.0d0
    do k=1,nnot
      s1 = s1+a(i,ind(k))*b(j,ind(k))
      s2 = s2+a(i,ind(k))*b(j+1,ind(k))
    end do
    r(i,j) = s1
    r(i,j+1) = s2
  else if (nr1 == 3) then
    s1 = 0.0d0
    s2 = 0.0d0
    s3 = 0.0d0
    do k=1,nnot
      s1 = s1+a(i,ind(k))*b(j,ind(k))
      s2 = s2+a(i,ind(k))*b(j+1,ind(k))
      s3 = s3+a(i,ind(k))*b(j+2,ind(k))
    end do
    r(i,j) = s1
    r(i,j+1) = s2
    r(i,j+2) = s3
  else if (nr1 > 3) then
    call WarningMessage(2,'nr1 > 3')
    call Abend()
  end if
100 continue
end do

if (mxind >= nlink) return
call WarningMessage(2,'MxInd < nLink')
write(6,*) 'mxind,nlink=',mxind,nlink
call Abend()

return

end subroutine ntmul
#endif
