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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine COUL(ISYP,ISYQ,ISYI,ISYJ,II,IJ,ERI,SCR)

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ISYP, ISYQ, ISYI, ISYJ, II, IJ
real(kind=wp), intent(_OUT_) :: ERI(*), SCR(*)
integer(kind=iwp) :: I3, I34, I4, IDISK, IS12, IS34, ISY1, ISY2, ISY3, ISY4, NBUF, NDIM2M, NO1, NO2
logical(kind=iwp) :: TRANSP, TRIANG
#include "intgrl.fh"

! Return a matrix ERI(IP,IQ) of two-electron integrals (pq,ij).
! IP=1..NORB(ISYP), IQ=1..NORB(ISYQ),1<=II<=NOSH(ISYI),
! 1<=IJ<=NOSH(ISYJ). Normal Fortran storage order.

NDIM2M = (NSYMZ*(NSYMZ+1))/2
TRANSP = .false.
TRIANG = .false.
if (ISYP >= ISYQ) then
  ISY1 = ISYP
  ISY2 = ISYQ
  TRANSP = .true.
  !TYPE = 'TRANSPOS'
  if (ISYP == ISYQ) then
    TRIANG = .true.
  end if
  !if (ISYP == ISYQ) TYPE = 'TRIANGUL'
else
  ISY1 = ISYQ
  ISY2 = ISYP
  !TYPE = 'NORMAL  '
end if
if (ISYI >= ISYJ) then
  ISY3 = ISYI
  ISY4 = ISYJ
  I3 = II
  I4 = IJ
  if ((ISYI == ISYJ) .and. (II < IJ)) then
    I3 = IJ
    I4 = II
  end if
else
  ISY3 = ISYJ
  ISY4 = ISYI
  I3 = IJ
  I4 = II
end if

! Disk address to symmetry block:
IS12 = (ISY1*(ISY1-1))/2+ISY2
IS34 = (ISY3*(ISY3-1))/2+ISY4
IDISK = IAD2M(1,IS12+NDIM2M*(IS34-1))

! Record number, in symmetry block:
I34 = I4+NOSHZ(ISY4)*(I3-1)
if (ISY3 == ISY4) I34 = (I3*(I3-1))/2+I4

! Buffer sizes:
NO1 = NORBZ(ISY1)
NO2 = NORBZ(ISY2)
NBUF = NO1*NO2
if (TRIANG) NBUF = (NBUF+NO1)/2
if (NBUF == 0) return

! Address update for earlier records, then read:
IDISK = IDISK+NBUF*(I34-1)

if (.not. TRANSP) then
  call dDAFILE(LUINTMZ,2,ERI,NBUF,IDISK)
else
  call dDAFILE(LUINTMZ,2,SCR,NBUF,IDISK)

  if (.not. TRIANG) then
    call TRNSPS(NO2,NO1,SCR,ERI)
  else
    call SQUARE(SCR,ERI,NO1,1,NO1)
  end if
end if

return

end subroutine COUL
