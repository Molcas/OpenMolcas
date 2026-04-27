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

function OVERLAP_RASSI(IFSBTAB1,IFSBTAB2,PSI1,PSI2)
! Purpose: Compute the overlap of two wave functions
! The FS blocks of the two wave functions:

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: OVERLAP_RASSI
integer(kind=iwp), intent(in) :: IFSBTAB1(*), IFSBTAB2(*)
real(kind=wp), intent(in) :: PSI1(*), PSI2(*)
integer(kind=iwp) :: IBLKPOS1, IBLKPOS2, IFSB1, IFSB2, ISSTARR(50), KHSH2, KPOS1, KPOS2, KSTARR1, KSTARR2, NASPRT1, NASPRT2, &
                     NBLKSIZ1, NBLKSIZ2, NDETS1, NDETS2, NFSB1, NFSB2, NHSH2
real(kind=wp), external :: ddot_

NFSB1 = IFSBTAB1(3)
NASPRT1 = IFSBTAB1(4)
NDETS1 = IFSBTAB1(5)
KSTARR1 = 8
NFSB2 = IFSBTAB2(3)
NASPRT2 = IFSBTAB2(4)
NDETS2 = IFSBTAB2(5)
NHSH2 = IFSBTAB2(6)
KHSH2 = IFSBTAB2(7)
KSTARR2 = 8
OVERLAP_RASSI = Zero
if (NFSB1 == 0) return
if (NFSB2 == 0) return
if (NASPRT1 /= NASPRT2) then
  write(u6,*) ' OVERLAP Error: The two wave function structures'
  write(u6,*) ' have different nr of subpartitions!'
  call ABEND()
end if
if (NDETS1 == 0) return
if (NDETS2 == 0) return

! Loop over FS blocks of the PSI1 wave function
do IFSB1=1,NFSB1
  KPOS1 = KSTARR1+(NASPRT1+2)*(IFSB1-1)
  ISSTARR(1:NASPRT1) = IFSBTAB1(KPOS1:KPOS1+NASPRT1-1)
  NBLKSIZ1 = IFSBTAB1(KPOS1+NASPRT1)
  IBLKPOS1 = IFSBTAB1(KPOS1+NASPRT1+1)
  ! Find this block in the PSI2 structure.
  call HSHGET(ISSTARR,NASPRT2,NASPRT2+2,IFSBTAB2(KSTARR2),NHSH2,IFSBTAB2(KHSH2),IFSB2)
  if (IFSB2 == 0) cycle
  KPOS2 = KSTARR2+(NASPRT2+2)*(IFSB2-1)
  NBLKSIZ2 = IFSBTAB2(KPOS2+NASPRT2)
  if (NBLKSIZ1 /= NBLKSIZ2) then
    write(u6,*) ' OVERLAP Error: The same FS block has not'
    write(u6,*) ' the same size in PSI1 and PSI2.'
    call ABEND()
  end if
  IBLKPOS2 = IFSBTAB2(KPOS2+NASPRT2+1)
  OVERLAP_RASSI = OVERLAP_RASSI+DDOT_(NBLKSIZ1,PSI1(IBLKPOS1),1,PSI2(IBLKPOS2),1)
end do

end function OVERLAP_RASSI
