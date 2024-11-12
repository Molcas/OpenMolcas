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

subroutine DMATRIX(E,F,Z,IPRINT)
! THIS ROUTINE COMPUTES THE USUAL D-TENSOR ON THE BASIS OF COEFFICIENTS
! OF THE STEVENS OPERATORS OF ORDER 2 (ES AND FS) AND DIAGONALIZE IT
! TO OBTAIN THE MAIN ANISOTROPY AXES

use Constants, only: Zero, One, Two, Half, OneHalf
use Definitions, only: wp, iwp, u6

implicit none
complex(kind=wp), intent(in) :: E(0:2), F(0:2)
real(kind=wp), intent(in) :: Z(3,3)
integer(kind=iwp), intent(in) :: iprint
integer(kind=iwp) :: INFO, J
real(kind=wp) :: CF, D_factor, daxes(3,3), diff12, diff23, DMATR(3,3), dtens(3), E_factor, SMAT(3,3), Unity(3,3), WD(3), ZD(3,3), &
                 ZD2(3,3)
complex(kind=wp) :: DMAT(3,3)

CF = sqrt(OneHalf)

DMAT(1,1) = CF*E(2)-E(0)
DMAT(2,2) = -CF*E(2)-E(0)
DMAT(3,3) = Two*E(0)
DMAT(1,2) = CF*F(2)
DMAT(1,3) = CF*E(1)
DMAT(2,3) = CF*F(1)
DMAT(2,1) = DMAT(1,2)
DMAT(3,1) = DMAT(1,3)
DMAT(3,2) = DMAT(2,3)

DMATR(:,:) = real(DMAT(:,:))

call DIAG_R2(DMATR,3,INFO,WD,ZD)

! calculate the rotation matrix:
call unitmat(Unity,3)

call DGEMM_('N','N',3,3,3,One,ZD,3,Unity,3,Zero,SMAT,3)

! Set the dtens and daxes with respect to the gtens and maxes.
!ccccccccccccccccccccccccccccccc
dtens(:) = Zero
if ((abs(SMAT(1,1)) > abs(SMAT(1,2))) .and. (abs(SMAT(1,1)) > abs(SMAT(1,3)))) then
  ! the WD(1) and ZD(i,1) correspond to gtens(1) and maxes(i,1)
  dtens(1) = WD(1)
  if (SMAT(1,1) > Zero) then
    if (IPRINT > 2) write(u6,'(a)') 'SMAT(1,1) is larger than SMAT(1,2) and SMAT(1,3) and is positive'
    daxes(:,1) = ZD(:,1)
  else
    if (IPRINT > 2) write(u6,'(a)') 'SMAT(1,1) is larger than SMAT(1,2) and SMAT(1,3) and is negative'
    daxes(:,1) = -ZD(:,1)
  end if

else if ((abs(SMAT(1,2)) > abs(SMAT(1,1))) .and. (abs(SMAT(1,2)) > abs(SMAT(1,3)))) then
  ! the WD(1) and ZD(i,1) correspond to gtens(2) and maxes(i,2)
  dtens(1) = WD(2)
  if (SMAT(1,2) > Zero) then
    if (IPRINT > 2) write(u6,'(a)') 'SMAT(1,2) is larger than SMAT(1,1) and SMAT(1,3) and is positive'
    daxes(:,1) = ZD(:,2)
  else
    if (IPRINT > 2) write(u6,'(a)') 'SMAT(1,2) is larger than SMAT(1,1) and SMAT(1,3) and is negative'
    daxes(:,1) = -ZD(:,2)
  end if

else if ((abs(SMAT(1,3)) > abs(SMAT(1,1))) .and. (abs(SMAT(1,3)) > abs(SMAT(1,2)))) then
  ! the WD(1) and ZD(i,1) correspond to gtens(3) and maxes(i,3)
  dtens(1) = WD(3)
  if (SMAT(1,3) > Zero) then
    if (IPRINT > 2) write(u6,'(a)') 'SMAT(1,3) is larger than SMAT(1,1) and SMAT(1,2) and is positive'
    daxes(:,1) = ZD(:,3)
  else
    if (IPRINT > 2) write(u6,'(a)') 'SMAT(1,3) is larger than SMAT(1,1) and SMAT(1,2) and is negative'
    daxes(:,1) = -ZD(:,3)
  end if
end if
!ccccccccccccccccccccccccccccccc

if ((abs(SMAT(2,1)) > abs(SMAT(2,2))) .and. (abs(SMAT(2,1)) > abs(SMAT(2,3)))) then
  ! the WD(2) and ZD(i,2) correspond to gtens(1) and maxes(i,1)
  dtens(2) = WD(1)
  if (SMAT(2,1) > Zero) then
    if (IPRINT > 2) write(u6,'(a)') 'SMAT(2,1) is larger than SMAT(2,2) and SMAT(2,3) and is positive'
    daxes(:,2) = ZD(:,1)
  else
    if (IPRINT > 2) write(u6,'(a)') 'SMAT(2,1) is larger than SMAT(2,2) and SMAT(2,3) and is negative'
    daxes(:,2) = -ZD(:,1)
  end if

else if ((abs(SMAT(2,2)) > abs(SMAT(2,1))) .and. (abs(SMAT(2,2)) > abs(SMAT(2,3)))) then
  ! the WD(2) and ZD(i,2) correspond to gtens(2) and maxes(i,2)
  dtens(2) = WD(2)
  if (SMAT(2,2) > Zero) then
    if (IPRINT > 2) write(u6,'(a)') 'SMAT(2,2) is larger than SMAT(2,1) and SMAT(2,3) and is positive'
    daxes(:,2) = ZD(:,2)
  else
    if (IPRINT > 2) write(u6,'(a)') 'SMAT(2,2) is larger than SMAT(2,1) and SMAT(2,3) and is negative'
    daxes(:,2) = -ZD(:,2)
  end if

else if ((abs(SMAT(2,3)) > abs(SMAT(2,1))) .and. (abs(SMAT(2,3)) > abs(SMAT(2,2)))) then
  ! the WD(2) and ZD(i,2) correspond to gtens(3) and maxes(i,3)
  dtens(2) = WD(3)
  if (SMAT(2,3) > Zero) then
    if (IPRINT > 2) write(u6,'(a)') 'SMAT(2,3) is larger than SMAT(2,1) and SMAT(2,2) and is positive'
    daxes(:,2) = ZD(:,3)
  else
    if (IPRINT > 2) write(u6,'(a)') 'SMAT(2,3) is larger than SMAT(2,1) and SMAT(2,2) and is negative'
    daxes(:,2) = -ZD(:,3)
  end if
end if

!ccccccccccccccccccccccccccccccc
if ((abs(SMAT(3,1)) > abs(SMAT(3,2))) .and. (abs(SMAT(3,1)) > abs(SMAT(3,3)))) then
  ! the WD(3) and ZD(i,3) correspond to gtens(1) and maxes(i,1)
  dtens(3) = WD(1)
  if (SMAT(3,1) > Zero) then
    if (IPRINT > 2) write(u6,'(a)') 'SMAT(3,1) is larger than SMAT(3,2) and SMAT(3,3) and is positive'
    daxes(:,3) = ZD(:,1)
  else
    if (IPRINT > 2) write(u6,'(a)') 'SMAT(3,1) is larger than SMAT(3,2) and SMAT(3,3) and is negative'
    daxes(:,3) = -ZD(:,1)
  end if

else if ((abs(SMAT(3,2)) > abs(SMAT(3,1))) .and. (abs(SMAT(3,2)) > abs(SMAT(3,3)))) then
  ! the WD(3) and ZD(i,3) correspond to gtens(2) and maxes(i,2)
  dtens(3) = WD(2)
  if (SMAT(3,2) > Zero) then
    if (IPRINT > 2) write(u6,'(a)') 'SMAT(3,2) is larger than SMAT(3,1) and SMAT(3,3) and is positive'
    daxes(:,3) = ZD(:,2)
  else
    if (IPRINT > 2) write(u6,'(a)') 'SMAT(3,2) is larger than SMAT(3,1) and SMAT(3,3) and is negative'
    daxes(:,3) = -ZD(:,2)
  end if

else if ((abs(SMAT(3,3)) > abs(SMAT(3,1))) .and. (abs(SMAT(3,3)) > abs(SMAT(3,2)))) then
  ! the WD(3) and ZD(i,3) correspond to gtens(3) and maxes(i,3)
  dtens(3) = WD(3)
  if (SMAT(3,3) > Zero) then
    if (IPRINT > 2) write(u6,'(a)') 'SMAT(3,3) is larger than SMAT(3,1) and SMAT(3,2) and is positive'
    daxes(:,3) = ZD(:,3)
  else
    if (IPRINT > 2) write(u6,'(a)') 'SMAT(3,3) is larger than SMAT(3,1) and SMAT(3,2) and is negative'
    daxes(:,3) = -ZD(:,3)
  end if
end if

call DGEMM_('N','N',3,3,3,One,Z,3,daxes,3,Zero,ZD2,3)

diff12 = abs(dtens(1)-dtens(2))
diff23 = abs(dtens(2)-dtens(3))

if (diff12 > diff23) then
  D_factor = OneHalf*dtens(1)
  E_factor = (dtens(2)-dtens(3))*Half
else
  D_factor = OneHalf*dtens(3)
  E_factor = (dtens(1)-dtens(2))*Half
end if

if (iprint > 2) then
  write(u6,'(20X,A)') 'D-TENSOR:'
  write(u6,*)
  write(u6,'(10X,A,10X,3(F9.5,2X))') '|  xx    xy    xz  |',(DMATR(1,J),J=1,3)
  write(u6,'(10X,A,10X,3(F9.5,2X))') '|  yx    yy    yz  |',(DMATR(2,J),J=1,3)
  write(u6,'(10X,A,10X,3(F9.5,2X))') '|  zx    zy    zz  |',(DMATR(3,J),J=1,3)
  write(u6,*)
end if

if (iprint >= 2) then
  write(u6,*)
  write(u6,'(A)') 'D TENSOR:'
  write(u6,'(2a)') repeat('-',84),'|'
  write(u6,'(A,4x,A,27x,A,21x,a,3x,a)') 'MAIN VALUES','|','MAIN ANISOTROPY AXES','|','x , y , z  -- initial Cartesian axes'
  write(u6,'(6a,3x,a)') repeat('-',15),'|',repeat('-',36),'|',repeat('-',31),'|','Xm, Ym, Zm -- main magnetic axes'
  write(u6,'(15x,a,4x,a,5x,a,8x,a,8x,a,4x,a,5x,a,9x,a,9x,a,5x,a,3x,a)') '|','|','Xm','Ym','Zm','|','x','y','z','|', &
                                                                        'Xa, Ya, Za -- main anisotropy axes'
  write(u6,'(8a)') repeat('-',15),'|',repeat('-',4),'|',repeat('-',31),'|',repeat('-',31),'|'
  write(u6,'(A,F9.3,A,3F10.6,1x,A,3F10.6,1x,A)') ' Dx =',dtens(1),' | Xa |',(daxes(j,1),j=1,3),'|',(ZD2(j,1),j=1,3),'|'
  write(u6,'(A,F9.3,A,3F10.6,1x,A,3F10.6,1x,A)') ' Dy =',dtens(2),' | Ya |',(daxes(j,2),j=1,3),'|',(ZD2(j,2),j=1,3),'|'
  write(u6,'(A,F9.3,A,3F10.6,1x,A,3F10.6,1x,A)') ' Dz =',dtens(3),' | Za |',(daxes(j,3),j=1,3),'|',(ZD2(j,3),j=1,3),'|'
  write(u6,'(2a)') repeat('-',84),'|'
  write(u6,*)
  write(u6,'(A)') '2-nd order ZFS Hamiltonian:'
  write(u6,*)
  write(u6,'(A)') 'H_zfs^{2}= D * [S_{Za}^2 - S*(S+1)/3] + E * [S_{Xa}^2 - S_{Ya}^2]'
  write(u6,*)
  write(u6,'(A)') 'Anisotropy parameters: D = 3/2 * Dz;  E = (Dx-Dy)/2;'
  write(u6,'(a,F9.4)') 'D = ',D_factor
  write(u6,'(a,F9.4)') 'E = ',E_factor
end if

return

end subroutine DMATRIX
