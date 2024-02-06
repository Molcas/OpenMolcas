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

subroutine atens(moment,dim,gtens,maxes,iprint)

use Constants, only: Zero, One, Twelve, Half, cZero
use Definitions, only: wp, u6

implicit none
! Calling variables:
integer, intent(in) :: dim, iprint
complex(kind=8), intent(in) :: moment(3,dim,dim)
real(kind=8), intent(out) :: gtens(3)
real(kind=8), intent(out) :: maxes(3,3)
!-----------------------------------------------------------------------
!  dim    -- size of the magnetic moment
!            dim = muliplicity of the pseuDospin ( 2*S+1, where S is the pseuDospin);
!  moment -- matrix of size (3,dim,dim) of the moment (magnetic, spin or angular)
!  gtens  -- array of size (3) keeping the main values of the A tensor ( sqrt(main_values) )
!  maxes  -- array of size (3,3) keeping the main axes of the A tensor Writen in
!            the right coordinate system (Determinant = +1)
!  iprint -- the print level of the Subroutine
!            iprint = 1 => no output
!            iprint = 2 => standard
!            iprint = 3 => print for debug
!-----------------------------------------------------------------------
! local variables
integer :: ic1, ic2, i, j, info
real(kind=8) :: A_TENS_TERM(3,3), W(3), MAIN(3), Z(3,3), factor, Det_gtens, diff12, diff23, ZR(3,3), dnorm
real(kind=8) :: dznrm2, FindDetR
complex(kind=8) :: AC_TENS(3,3), trace
external :: dznrm2, FindDetR, trace

dnorm = Zero
dnorm = dznrm2(3*dim*dim,moment,1)

if (dnorm <= tiny(dnorm)) then
  write(u6,'(A)') 'Norm of the magnetic moment is zero.'
  write(u6,'(A)') 'Returning the default (dummy) values'
  gtens = Zero
  maxes = Zero
  do i=1,3
    maxes(i,i) = One
  end do
  return
end if

! initialization:
do I=1,3
  do J=1,3
    AC_TENS(I,J) = cZero
    A_TENS_TERM(I,J) = Zero
  end do
end do

do ic1=1,3
  do ic2=1,3
    Ac_tens(ic1,ic2) = trace(dim,moment(ic1,1:dim,1:dim),moment(ic2,1:dim,1:dim))
  end do
end do

do ic1=1,3
  do ic2=1,3
    A_TENS_TERM(ic1,ic2) = real(Ac_tens(ic1,ic2)+Ac_tens(ic2,ic1))*Half
  end do
end do

factor = Zero
factor = Twelve/real(dim**3-dim,kind=wp)
do ic1=1,3
  do ic2=1,3
    A_TENS_TERM(ic1,ic2) = factor*A_TENS_TERM(ic1,ic2)
  end do
end do

if (iprint > 2) then
  write(u6,'(/)')
  write(u6,'(5X,A)') 'A_TENS_TERM(ic1,ic2):'
  write(u6,*)
  do ic1=1,3
    write(u6,'(5X,3(2F14.7,3x))') (A_TENS_TERM(ic1,ic2),ic2=1,3)
  end do
end if

! Diagonalization of A_tens - g tensors

do I=1,3
  main(I) = Zero
  w(I) = Zero
  do J=1,3
    z(I,J) = Zero
  end do
end do
info = 0

call DIAG_R2(A_TENS_TERM(1:3,1:3),3,info,W(1:3),Z(1:3,1:3))

if (INFO /= 0) go to 199
if ((w(1) < Zero) .and. (w(2) < Zero) .and. (w(3) < Zero)) then
  write(u6,'(2x,A)') 'ALL EIGENVALUES OF THE A-TENSOR ARE NEGATIVE'
  write(u6,'(2X,A)') 'THIS IS A VERY UNUSUAL SITUATION. PLEASE CHECK MANUALLY '
  write(u6,'(2x,A)') 'THE FOLLOWING PART OF THE PSEUDoSPIN SECTION'
  write(u6,'(2x,A)') 'MUST BE DISREGARDED. THE RESULTS ARE NOT TRUSTABLE.'
  go to 199
end if

if (iprint > 2) then
  write(u6,*)
  write(u6,'(4x,A)') 'A_TENS_TERM TENSOR:'
  write(u6,'(65a)') ('-',i=1,56),'|'
  write(u6,'(4x,A,4x,A,13x,A,5x,a,3x,a)') 'MAIN VALUES','|','MAIN MAGNETIC AXES','|','x , y , z  -- initial Cartesian axes'
  write(u6,'(57a,3x,a)') ('-',i=1,19),'|',('-',i=1,36),'|','Xm, Ym, Zm -- main magnetic axes'
  write(u6,'(19x,a,4x,a,5x,a,9x,a,9x,a,5x,a)') '|','|','x','y','z','|'
  write(u6,'(65a)') ('-',i=1,19),'|',('-',i=1,4),'|',('-',i=1,31),'|'
  write(u6,'(A,F12.9,A,3F10.6,1x,A)') ' gX = ',w(1),' | Xm |',(z(j,1),j=1,3),'|'
  write(u6,'(A,F12.9,A,3F10.6,1x,A)') ' gY = ',w(2),' | Ym |',(z(j,2),j=1,3),'|'
  write(u6,'(A,F12.9,A,3F10.6,1x,A)') ' gZ = ',w(3),' | Zm |',(z(j,3),j=1,3),'|'
  write(u6,'(65a)') ('-',i=1,56),'|'
end if

do I=1,3
  if (W(I) < Zero) W(I) = 1.0e-24_wp
  MAIN(i) = sqrt(W(i))
end do

! Check the sign of the coordinate system. If CS is Left-handed,
! Then change it to RIGHT-handed
Det_gtens = Zero
do I=1,3
  do J=1,3
    ZR(I,J) = Zero
    ZR(I,J) = Z(I,J)
  end do
end do

Det_gtens = FindDetR(ZR,3)

if (Det_gtens < Zero) then
  do i=1,3
    Z(i,1) = -Z(i,1)
  end do
  if (iprint > 2) then
    write(u6,'(a)') 'The coordinate system is LEFT-handed.'
    write(u6,'(a)') 'It has been changed to RIGHT-handed'
  end if
end if

diff12 = Zero
diff23 = Zero
diff12 = MAIN(2)-MAIN(1)
diff23 = MAIN(3)-MAIN(2)

if (iprint > 2) then
  write(u6,'(5x,a,3F19.15)') 'diff12 = ',diff12
  write(u6,'(5x,a,3F19.15)') 'diff23 = ',diff23
end if

do i=1,3
  gtens(i) = Zero
  do j=1,3
    maxes(i,j) = Zero
  end do
end do
! set the main Z axis:
if (diff12 > diff23) then
  gtens(3) = MAIN(1)
  gtens(2) = MAIN(2)
  gtens(1) = MAIN(3)

  if (Z(3,1) >= Zero) then
    do i=1,3
      maxes(i,3) = Z(i,1)
      maxes(i,1) = Z(i,3)
    end do
  else if (Z(3,1) < Zero) then
    do i=1,3
      maxes(i,3) = -Z(i,1)
      maxes(i,1) = Z(i,3)
    end do
  end if

  maxes(1,2) = maxes(2,3)*maxes(3,1)-maxes(2,1)*maxes(3,3)
  maxes(2,2) = maxes(1,1)*maxes(3,3)-maxes(1,3)*maxes(3,1)
  maxes(3,2) = maxes(1,3)*maxes(2,1)-maxes(1,1)*maxes(2,3)

else if (diff23 > diff12) then
  gtens(3) = MAIN(3)
  gtens(2) = MAIN(2)
  gtens(1) = MAIN(1)

  if (Z(3,3) >= Zero) then
    do i=1,3
      maxes(i,3) = Z(i,3)
      maxes(i,1) = Z(i,1)
    end do
  else if (Z(3,3) < Zero) then
    do i=1,3
      maxes(i,3) = -Z(i,3)
      maxes(i,1) = Z(i,1)
    end do
  end if

  maxes(1,2) = maxes(2,3)*maxes(3,1)-maxes(2,1)*maxes(3,3)
  maxes(2,2) = maxes(1,1)*maxes(3,3)-maxes(1,3)*maxes(3,1)
  maxes(3,2) = maxes(1,3)*maxes(2,1)-maxes(1,1)*maxes(2,3)
else !( diff23 == diff12)
  ! this special case is isotropic:
  ! therefore assign main axes as close as to be to the cartesian xyz
  do i=1,3
    gtens(i) = MAIN(i)
    maxes(i,i) = One
  end do

end if
if (iprint > 2) then
  write(u6,*)
  write(u6,'(20X,A)') 'A-TENSOR:'
  write(u6,*)
  write(u6,'(10X,A,10X,3(F11.5,2X))') '|  xx    xy    xz  |',(A_TENS_TERM(1,i),i=1,3)
  write(u6,'(10X,A,10X,3(F11.5,2X))') '|  yx    yy    yz  |',(A_TENS_TERM(2,i),i=1,3)
  write(u6,'(10X,A,10X,3(F11.5,2X))') '|  zx    zy    zz  |',(A_TENS_TERM(3,i),i=1,3)
end if

if (iprint >= 2) then
  write(u6,*)
  write(u6,'(4x,A)') 'g TENSOR:'
  write(u6,'(65a)') ('-',i=1,56),'|'
  write(u6,'(4x,A,4x,A,13x,A,5x,a,3x,a)') 'MAIN VALUES','|','MAIN MAGNETIC AXES','|','x , y , z  -- initial Cartesian axes'
  write(u6,'(26a,3x,a)') ('-',i=1,19),'|',('-',i=1,4),'|','----- x ------- y ------- z ---|','Xm, Ym, Zm -- main magnetic axes'
  !write(u6,'(A,F12.9,A,3F10.6,1x,A)') ' gX = ',gtens(1),' | Xm |',(maxes(j,1),j=1,3),'|'
  !write(u6,'(A,F12.9,A,3F10.6,1x,A)') ' gY = ',gtens(2),' | Ym |',(maxes(j,2),j=1,3),'|'
  !write(u6,'(A,F12.9,A,3F10.6,1x,A)') ' gZ = ',gtens(3),' | Zm |',(maxes(j,3),j=1,3),'|'
  write(u6,'(A,F18.14,A,3F18.14,1x,A)') ' gX = ',gtens(1),' | Xm |',(maxes(j,1),j=1,3),'|'
  write(u6,'(A,F18.14,A,3F18.14,1x,A)') ' gY = ',gtens(2),' | Ym |',(maxes(j,2),j=1,3),'|'
  write(u6,'(A,F18.14,A,3F18.14,1x,A)') ' gZ = ',gtens(3),' | Zm |',(maxes(j,3),j=1,3),'|'
  write(u6,'(65a)') ('-',i=1,56),'|'
end if

199 continue

return

end subroutine atens
