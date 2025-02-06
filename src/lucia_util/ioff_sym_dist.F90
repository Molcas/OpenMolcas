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

function IOFF_SYM_DIST(ISYM,NGASL,IOFF,MXVAL,MNVAL)
! A ts block of string is given and the individual
! symmetrydisrtributions has been obtained ( for example
! by TS_SYM_PNT)
!
! Obtain offset for symmetrycombination defined by ISYM

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: IOFF_SYM_DIST
integer(kind=iwp) :: ISYM(*), NGASL, IOFF(*), MXVAL(*), MNVAL(*)
integer(kind=iwp) :: I, IGAS, IMULT, J, NTEST

! Address in IOFF is
!     1
!     +  (ISM1-MNVAL(1))
!     +  (ISM2-MNVAL(2))*(MXVAL(1)-MNVAL(1)+1)
!     +  (ISM3-MNVAL(3))*(MXVAL(1)-MNVAL(1)+1)*(MXVAL(2)-MNVAL(2)+1)
!     +  ...
!     +  (ISM L-1-MNVAL(L-1))*Prod(i=1,L-2)(MXVAL(i)-MNVAL(i)+1)

NTEST = 0
if (NTEST >= 100) then
  write(u6,*) ' Isym, mnval, ioff:'
  call iwrtma(isym,1,ngasl,1,ngasl)
  call iwrtma(MNVAL,1,ngasl,1,ngasl)
  call iwrtma(ioff,1,ngasl,1,ngasl)
end if
! Offset for this symmetry distribution in IOFFI
I = 1
IMULT = 1
do IGAS=1,NGASL-1
  I = I+(ISYM(IGAS)-MNVAL(IGAS))*IMULT
  IMULT = IMULT*(MXVAL(IGAS)-MNVAL(IGAS)+1)
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
