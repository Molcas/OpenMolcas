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

subroutine atens(moment,d,gtens,maxes,iprint)
!-----------------------------------------------------------------------
! d      -- size of the magnetic moment
!           d = muliplicity of the pseudospin ( 2*S+1, where S is the pseudospin);
! moment -- matrix of size (3,d,d) of the moment (magnetic, spin or angular)
! gtens  -- array of size (3) keeping the main values of the A tensor ( sqrt(main_values) )
! maxes  -- array of size (3,3) keeping the main axes of the A tensor Writen in
!           the right coordinate system (Determinant = +1)
! iprint -- the print level of the Subroutine
!           iprint = 1 => no output
!           iprint = 2 => standard
!           iprint = 3 => print for debug
!-----------------------------------------------------------------------

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Twelve, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: d, iprint
complex(kind=wp), intent(in) :: moment(3,d,d)
real(kind=wp), intent(out) :: gtens(3)
real(kind=wp), intent(out) :: maxes(3,3)
integer(kind=iwp) :: i, ic1, ic2, info, j
real(kind=wp) :: A_TENS_TERM(3,3), Det_gtens, diff12, diff23, dnorm, MAIN(3), W(3), Z(3,3), ZR(3,3)
complex(kind=wp) :: AC_TENS(3,3)
complex(kind=wp), allocatable :: tmp1(:,:), tmp2(:,:)
real(kind=wp), external :: dznrm2, FindDetR
complex(kind=wp), external :: trace

dnorm = Zero
dnorm = dznrm2(3*d*d,moment,1)

if (dnorm <= tiny(dnorm)) then
  write(u6,'(A)') 'Norm of the magnetic moment is zero.'
  write(u6,'(A)') 'Returning the default (dummy) values'
  gtens(:) = Zero
  call unitmat(maxes,3)
  return
end if

! initialization:

call mma_allocate(tmp1,d,d,label='tmp1')
call mma_allocate(tmp2,d,d,label='tmp1')
do ic1=1,3
  tmp1(:,:) = moment(ic1,:,:)
  do ic2=1,3
    tmp2(:,:) = moment(ic2,:,:)
    Ac_tens(ic1,ic2) = trace(d,tmp1,tmp2)
  end do
end do
call mma_deallocate(tmp1)
call mma_deallocate(tmp2)

do ic1=1,3
  do ic2=1,3
    A_TENS_TERM(ic1,ic2) = real(Ac_tens(ic1,ic2)+Ac_tens(ic2,ic1))*Half
  end do
end do

A_TENS_TERM(:,:) = Twelve/real(d**3-d,kind=wp)*A_TENS_TERM(:,:)

if (iprint > 2) then
  write(u6,'(/)')
  write(u6,'(5X,A)') 'A_TENS_TERM(ic1,ic2):'
  write(u6,*)
  do ic1=1,3
    write(u6,'(5X,3(2F14.7,3x))') (A_TENS_TERM(ic1,ic2),ic2=1,3)
  end do
end if

! Diagonalization of A_tens - g tensors

call DIAG_R2(A_TENS_TERM,3,info,W,Z)

if (INFO /= 0) return
if ((w(1) < Zero) .and. (w(2) < Zero) .and. (w(3) < Zero)) then
  write(u6,'(2x,A)') 'ALL EIGENVALUES OF THE A-TENSOR ARE NEGATIVE'
  write(u6,'(2X,A)') 'THIS IS A VERY UNUSUAL SITUATION. PLEASE CHECK MANUALLY '
  write(u6,'(2x,A)') 'THE FOLLOWING PART OF THE PSEUDOSPIN SECTION'
  write(u6,'(2x,A)') 'MUST BE DISREGARDED. THE RESULTS ARE NOT TRUSTABLE.'
  return
end if

if (iprint > 2) then
  write(u6,*)
  write(u6,'(4x,A)') 'A_TENS_TERM TENSOR:'
  write(u6,'(2a)') repeat('-',56),'|'
  write(u6,'(4x,A,4x,A,13x,A,5x,a,3x,a)') 'MAIN VALUES','|','MAIN MAGNETIC AXES','|','x , y , z  -- initial Cartesian axes'
  write(u6,'(4a,3x,a)') repeat('-',19),'|',repeat('-',36),'|','Xm, Ym, Zm -- main magnetic axes'
  write(u6,'(19x,a,4x,a,5x,a,9x,a,9x,a,5x,a)') '|','|','x','y','z','|'
  write(u6,'(6a)') repeat('-',19),'|',repeat('-',4),'|',repeat('-',31),'|'
  write(u6,'(A,F12.9,A,3F10.6,1x,A)') ' gX = ',w(1),' | Xm |',(z(j,1),j=1,3),'|'
  write(u6,'(A,F12.9,A,3F10.6,1x,A)') ' gY = ',w(2),' | Ym |',(z(j,2),j=1,3),'|'
  write(u6,'(A,F12.9,A,3F10.6,1x,A)') ' gZ = ',w(3),' | Zm |',(z(j,3),j=1,3),'|'
  write(u6,'(2a)') repeat('-',56),'|'
end if

do I=1,3
  if (W(I) < Zero) W(I) = 1.0e-24_wp
end do
MAIN(:) = sqrt(W(:))

! Check the sign of the coordinate system. If CS is Left-handed,
! Then change it to RIGHT-handed
ZR(:,:) = Z(:,:)

Det_gtens = FindDetR(ZR,3)

if (Det_gtens < Zero) then
  Z(:,1) = -Z(:,1)
  if (iprint > 2) then
    write(u6,'(a)') 'The coordinate system is LEFT-handed.'
    write(u6,'(a)') 'It has been changed to RIGHT-handed'
  end if
end if

diff12 = MAIN(2)-MAIN(1)
diff23 = MAIN(3)-MAIN(2)

if (iprint > 2) then
  write(u6,'(5x,a,3F19.15)') 'diff12 = ',diff12
  write(u6,'(5x,a,3F19.15)') 'diff23 = ',diff23
end if

maxes(:,:) = Zero
! set the main Z axis:
if (diff12 > diff23) then
  gtens(:) = MAIN(3:1:-1)

  if (Z(3,1) >= Zero) then
    maxes(:,3) = Z(:,1)
    maxes(:,1) = Z(:,3)
  else if (Z(3,1) < Zero) then
    maxes(:,3) = -Z(:,1)
    maxes(:,1) = Z(:,3)
  end if

  maxes(1,2) = maxes(2,3)*maxes(3,1)-maxes(2,1)*maxes(3,3)
  maxes(2,2) = maxes(1,1)*maxes(3,3)-maxes(1,3)*maxes(3,1)
  maxes(3,2) = maxes(1,3)*maxes(2,1)-maxes(1,1)*maxes(2,3)

else if (diff23 > diff12) then
  gtens(:) = MAIN(:)

  if (Z(3,3) >= Zero) then
    maxes(:,3) = Z(:,3)
    maxes(:,1) = Z(:,1)
  else if (Z(3,3) < Zero) then
    maxes(:,3) = -Z(:,3)
    maxes(:,1) = Z(:,1)
  end if

  maxes(1,2) = maxes(2,3)*maxes(3,1)-maxes(2,1)*maxes(3,3)
  maxes(2,2) = maxes(1,1)*maxes(3,3)-maxes(1,3)*maxes(3,1)
  maxes(3,2) = maxes(1,3)*maxes(2,1)-maxes(1,1)*maxes(2,3)

else !( diff23 == diff12)
  ! this special case is isotropic:
  ! therefore assign main axes as close as to be to the cartesian xyz
  gtens(:) = MAIN(:)
  call unitmat(maxes,3)

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
  write(u6,'(2a)') repeat('-',56),'|'
  write(u6,'(4x,A,4x,A,13x,A,5x,a,3x,a)') 'MAIN VALUES','|','MAIN MAGNETIC AXES','|','x , y , z  -- initial Cartesian axes'
  write(u6,'(4a,3x,a)') repeat('-',19),'|',repeat('-',4),'|','----- x ------- y ------- z ---|','Xm, Ym, Zm -- main magnetic axes'
  !write(u6,'(A,F12.9,A,3F10.6,1x,A)') ' gX = ',gtens(1),' | Xm |',(maxes(j,1),j=1,3),'|'
  !write(u6,'(A,F12.9,A,3F10.6,1x,A)') ' gY = ',gtens(2),' | Ym |',(maxes(j,2),j=1,3),'|'
  !write(u6,'(A,F12.9,A,3F10.6,1x,A)') ' gZ = ',gtens(3),' | Zm |',(maxes(j,3),j=1,3),'|'
  write(u6,'(A,F18.14,A,3F18.14,1x,A)') ' gX = ',gtens(1),' | Xm |',(maxes(j,1),j=1,3),'|'
  write(u6,'(A,F18.14,A,3F18.14,1x,A)') ' gY = ',gtens(2),' | Ym |',(maxes(j,2),j=1,3),'|'
  write(u6,'(A,F18.14,A,3F18.14,1x,A)') ' gZ = ',gtens(3),' | Zm |',(maxes(j,3),j=1,3),'|'
  write(u6,'(2a)') repeat('-',56),'|'
end if

return

end subroutine atens
