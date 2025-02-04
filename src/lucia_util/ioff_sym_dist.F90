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

integer function IOFF_SYM_DIST(ISYM,NGASL,IOFF,MAXVAL,MINVAL)
! A ts block of string is given and the individual
! symmetrydisrtributions has been obtained ( for example
! by TS_SYM_PNT)
!
! Obtain offset for symmetrycombination defined by ISYM

use Definitions, only: u6

implicit real*8(A-H,O-Z)
integer ISYM(*), IOFF(*), maxval(*), minval(*)

! Address in IOFF is
!     1
!     +  (ISM1-MINVAL(1))
!     +  (ISM2-MINVAL(2))*(MAXVAL(1)-MINVAL(1)+1)
!     +  (ISM3-MINVAL(3))*(MAXVAL(1)-MINVAL(1)+1)*(MAXVAL(2)-MINVAL(2)+1)
!     +  ...
!     +  (ISM L-1-MINVAL(L-1))*Prod(i=1,L-2)(MAXVAL(i)-MINVAL(i)+1)

NTEST = 0
if (NTEST >= 100) then
  write(u6,*) ' Isym, minval, ioff:'
  call iwrtma(isym,1,ngasl,1,ngasl)
  call iwrtma(minval,1,ngasl,1,ngasl)
  call iwrtma(ioff,1,ngasl,1,ngasl)
end if
! Offset for this symmetry distribution in IOFFI
I = 1
IMULT = 1
do IGAS=1,NGASL-1
  I = I+(ISYM(IGAS)-minval(IGAS))*IMULT
  IMULT = IMULT*(maxval(IGAS)-minval(IGAS)+1)
  !write(u6,*) ' igas,i,imult ',igas,i,imult
end do
! The following IF block is needed for avoinging going outside the IOFF bounds.
! This is possible for certain GAS setups. Test 897 helped in finding this issue.
if (I <= 0) then
  IOFF_SYM_DIST = 0
else
  IOFF_SYM_DIST = IOFF(I)
end if

if (NTEST >= 100) then
  write(u6,*) ' Info from IOFF_SYM_DIST'
  write(u6,*) ' ======================='
  write(u6,*)
  write(u6,*) ' Address and offset ',I,IOFF_SYM_DIST
  write(u6,*) ' Symmetry distribution : ',(ISYM(J),J=1,NGASL)
end if

end function IOFF_SYM_DIST
