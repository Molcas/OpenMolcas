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

function regge3j(j1,j2,j3,m1,m2,m3)
!bs uses magic square of regge (see Lindner pp. 38-39)
!bs
!bs  ---                                            ---
!bs |                                                  |
!bs | -j1+j2+j3     j1-j2+j3         j1+j2-j3          |
!bs |                                                  |
!bs |                                                  |
!bs |  j1-m1        j2-m2            j3-m3             |
!bs |                                                  |
!bs |                                                  |
!bs |  j1+m1        j2+m2            j3+m3             |
!bs |                                                  |
!bs  ---                                            ---

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: regge3j
integer(kind=iwp), intent(in) :: j1, j2, j3, m1, m2, m3
integer(kind=iwp), parameter :: mxLinRE = 36, nprim = 11
integer(kind=iwp) :: I, ibm, icheck, icoeff, ICOL, icount, Idenom, IDUMMY, IFIRST, imaxi, imini, IROW, ISECOND, Isigma, isgn, &
                     isum, iwork(nprim), J, LIMIT, MAT(3,3)
real(kind=wp) :: down, factor, up
!BS logical(kind=iwp) :: testup,testdown
! decompose factorials into powers of prime numbers
integer(kind=iwp), parameter :: facul(nprim,0:mxLinRE) = reshape([0,0,0,0,0,0,0,0,0,0,0,   &
                                                                  0,0,0,0,0,0,0,0,0,0,0,   &
                                                                  1,0,0,0,0,0,0,0,0,0,0,   &
                                                                  1,1,0,0,0,0,0,0,0,0,0,   &
                                                                  3,1,0,0,0,0,0,0,0,0,0,   &
                                                                  3,1,1,0,0,0,0,0,0,0,0,   &
                                                                  4,2,1,0,0,0,0,0,0,0,0,   &
                                                                  4,2,1,1,0,0,0,0,0,0,0,   &
                                                                  7,2,1,1,0,0,0,0,0,0,0,   &
                                                                  7,4,1,1,0,0,0,0,0,0,0,   &
                                                                  8,4,2,1,0,0,0,0,0,0,0,   &
                                                                  8,4,2,1,1,0,0,0,0,0,0,   &
                                                                  10,5,2,1,1,0,0,0,0,0,0,  &
                                                                  10,5,2,1,1,1,0,0,0,0,0,  &
                                                                  11,5,2,2,1,1,0,0,0,0,0,  &
                                                                  11,6,3,2,1,1,0,0,0,0,0,  &
                                                                  15,6,3,2,1,1,0,0,0,0,0,  &
                                                                  15,6,3,2,1,1,1,0,0,0,0,  &
                                                                  16,8,3,2,1,1,1,0,0,0,0,  &
                                                                  16,8,3,2,1,1,1,1,0,0,0,  &
                                                                  18,8,4,2,1,1,1,1,0,0,0,  &
                                                                  18,9,4,3,1,1,1,1,0,0,0,  &
                                                                  19,9,4,3,2,1,1,1,0,0,0,  &
                                                                  19,9,4,3,2,1,1,1,1,0,0,  &
                                                                  22,10,4,3,2,1,1,1,1,0,0, &
                                                                  22,10,6,3,2,1,1,1,1,0,0, &
                                                                  23,10,6,3,2,2,1,1,1,0,0, &
                                                                  23,13,6,3,2,2,1,1,1,0,0, &
                                                                  25,13,6,4,2,2,1,1,1,0,0, &
                                                                  25,13,6,4,2,2,1,1,1,1,0, &
                                                                  26,14,7,4,2,2,1,1,1,1,0, &
                                                                  26,14,7,4,2,2,1,1,1,1,1, &
                                                                  31,14,7,4,2,2,1,1,1,1,1, &
                                                                  31,15,7,4,3,2,1,1,1,1,1, &
                                                                  32,15,7,4,3,3,1,1,1,1,1, &
                                                                  32,15,8,5,3,3,1,1,1,1,1, &
                                                                  34,17,8,5,3,3,1,1,1,1,1  &
                                                                 ],shape(facul)), &
                                ihigh(0:mxLinRE) = [0,0,1,2,2,3,3,4,4,4,4,5,5,6,6,6,6,7,7,8,8,8,8,9,9,9,9,9,9,10,10,11,11,11,11, &
                                                    11,11], &
                                prim(nprim) = [2,3,5,7,11,13,17,19,23,29,31]        !prime numbers
!bs facul,   integer array (nprim,0:mxLinRE) prime-expansion of factorials
!bs mxLinRE, integer max. number for facul is given
!bs nprim,   number of primes for expansion of factorials
!bs prim,    integer array with the first nprim prime numbers
!bs iwork)   integer array of size nprim

regge3j = Zero
!write(u6,'(A24,6I3)') '3J to be calculated for ',j1,j2,j3,m1,m2,m3
!bs quick check  if =/= 0 at all
icheck = m1+m2+m3
if (icheck /= 0) then
  !write(u6,*) 'sum over m =/= 0'
  return
end if
!bs check triangular relation (|j1-j2|<= j3 <= j1+j2 )
imini = abs(j1-j2)
imaxi = j1+j2
if ((j3 < imini) .or. (j3 > imaxi)) then
  !write(u6,*) 'triangular relation not fulfilled'
  return
end if
!bs quick check  if =/= 0 at all  end
!bs
!bs 3J-symbol is not zero by simple rules
!bs
!bs initialize MAT
MAT(1,1) = -j1+j2+j3
MAT(2,1) = j1-m1
MAT(3,1) = j1+m1
MAT(1,2) = j1-j2+j3
MAT(2,2) = j2-m2
MAT(3,2) = j2+m2
MAT(1,3) = j1+j2-j3
MAT(2,3) = j3-m3
MAT(3,3) = j3+m3
do I=1,3
  do J=1,3
    !bs check for even numbers (2*integer) and positive or zero
    if ((mod(MAT(J,I),2) /= 0) .or. (MAT(J,I) < 0)) then
      !write(u6,*) 'J,I,MAT(J,I): ',J,I,MAT(J,I)
      return
    end if
    MAT(J,I) = MAT(J,I)/2
    if (Mat(j,i) > mxLinRE) call SysAbendMsg('regge3j','increase mxLinRE for regge3j',' ')
  end do
end do
Isigma = (j1+j2+j3)/2
!bs check the magic sums
do I=1,3
  IROW = 0
  ICOL = 0
  do J=1,3
    IROW = IROW+MAT(I,J)
    ICOL = ICOL+MAT(J,I)
  end do
  if ((IROW /= Isigma) .or. (ICOL /= Isigma)) then
    !write(u6,*) 'I,IROW,ICOL ',I,IROW,ICOL
    return
  end if
end do
!bs if j1+j2+j3 is odd: check for equal rows or columns
isgn = 1
if (abs(mod(Isigma,2)) == 1) then
  isgn = -1
  do I=1,3
    do J=I+1,3
      if ((MAT(1,I) == MAT(1,J)) .and. (MAT(2,I) == MAT(2,J)) .and. (MAT(3,I) == MAT(3,J))) return
      if ((MAT(I,1) == MAT(J,1)) .and. (MAT(I,2) == MAT(J,2)) .and. (MAT(I,3) == MAT(J,3))) return
    end do
  end do
end if
!bs look for the lowest element indices: IFIRST,ISECOND
imini = MAT(1,1)
IFIRST = 1
ISECOND = 1
do I=1,3
  do J=1,3
    if (MAT(J,I) < imini) then
      IFIRST = J
      ISECOND = I
      imini = MAT(J,I)
    end if
  end do
end do
!write(u6,*) 'Matrix before commuting vectors'
!do ibm=1,3
!  write(u6,'(3I5)') (Mat(ibm,j),j=1,3)
!end do
if (IFIRST /= 1) then  !interchange rows
  !write(u6,*) 'IFIRST = ',ifirst
  do I=1,3
    IDUMMY = MAT(1,I)
    MAT(1,I) = MAT(IFIRST,I)
    MAT(IFIRST,I) = IDUMMY
  end do
end if
if (ISECOND /= 1) then  !interchange columns
  !write(u6,*) 'ISECOND = ',isecond
  do I=1,3
    IDUMMY = MAT(I,1)
    MAT(I,1) = MAT(I,ISECOND)
    MAT(I,ISECOND) = IDUMMY
  end do
end if
!bs lowest element is now on (1,1)
!write(u6,*) 'Matrix after commuting vectors'
!do ibm=1,3
!  write(u6,'(3I5)') (Mat(ibm,j),j=1,3)
!end do
!bs begin to calculate Sum over s_n
!bs first the simple cases
if (Mat(1,1) == 0) then
  isum = 1
else if (Mat(1,1) == 1) then
  isum = Mat(2,3)*Mat(3,2)-Mat(2,2)*Mat(3,3)
else if (Mat(1,1) == 2) then
  isum = Mat(2,3)*(Mat(2,3)-1)*Mat(3,2)*(Mat(3,2)-1)-2*Mat(2,3)*Mat(3,2)*Mat(2,2)*Mat(3,3)+Mat(2,2)*(Mat(2,2)-1)*Mat(3,3)* &
         (Mat(3,3)-1)
else !  all the cases with Mat(1,1) >= 3
  icoeff = 1
  do ibm=Mat(3,2)-Mat(1,1)+1,Mat(3,2)
    icoeff = icoeff*ibm
  end do
  do ibm=Mat(2,3)-Mat(1,1)+1,Mat(2,3)
    icoeff = icoeff*ibm
  end do
  isum = icoeff
  do icount=1,MAT(1,1)
    icoeff = -icoeff*(Mat(1,1)+1-icount)*(Mat(2,2)+1-icount)*(Mat(3,3)+1-icount)
    Idenom = icount*(Mat(2,3)-Mat(1,1)+icount)*(Mat(3,2)-Mat(1,1)+icount)
    icoeff = icoeff/Idenom
    isum = isum+icoeff
  end do
end if
!bs additional sign from interchanging rows or columns
if (ifirst /= 1) isum = isum*isgn
if (isecond /= 1) isum = isum*isgn
!write(u6,*) 'isum = ',isum
!bs Mat(2,3)+Mat(3,2)
!bs (-)
if (abs(mod((Mat(2,3)+Mat(3,2)),2)) == 1) isum = -isum
!bs final factor
LIMIT = ihigh(max(Mat(1,1),Mat(1,2),Mat(1,3),Mat(2,1),Mat(2,2),Mat(2,3),Mat(3,1),Mat(3,2),Mat(3,3),(Isigma+1)))
do I=1,LIMIT
  iwork(I) = facul(I,Mat(1,2))+facul(I,Mat(2,1))+facul(I,Mat(3,1))+facul(I,Mat(1,3))-facul(I,Mat(1,1))-facul(I,Mat(2,2))- &
             facul(I,Mat(3,3))-facul(I,(Isigma+1))-facul(I,Mat(2,3))-facul(I,Mat(3,2))
end do
!write(u6,*) 'Iwork: ',(iwork(i),i=1,LIMIT)
factor = One
!bs iup = 1
!BS idown = 1
!BS testup = .true.
!BS testdown = .true.
!BS do I=1,LIMIT
!BS   do J=1,iwork(I)
!BS     iup = iup*prim(i)
!BS     if (iup < 0) testup = .false. !check for Integer overflow
!BS   end do
!BS end do
!BS up = real(iup,kind=wp)
!BS if  (.not. testup) then ! if the integers did not run correctly
up = One
do I=1,LIMIT
  do J=1,iwork(I)
    up = up*real(prim(i),kind=wp)
  end do
end do
!BS endif
!BS do I=1,LIMIT
!BS   do J=1,-iwork(I)
!BS     idown = idown*prim(i)
!BS     if (idown < 0) testdown = .false.
!BS   end do
!BS end do
!BS down = real(idown,kind=wp)
!BS if (.not. testdown) then
down = One
do I=1,LIMIT
  do J=1,-iwork(I)
    down = down*real(prim(i),kind=wp)
  end do
end do
!BS endif
!if (.not. (testup .and. testdown)) then
!  write(u6,*) 'j1,j2,j3,m1,m2,m3 ',j1,j2,j3,m1,m2,m3
!  write(u6,*) 'iup,idown ',iup,idown,'up,down ',up,down
!end if
factor = factor*up/down
!bs final result
regge3j = sqrt(factor)*real(isum,kind=wp)

return

end function regge3j
