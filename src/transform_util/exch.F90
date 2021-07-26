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

subroutine EXCH(ISYP,ISYI,ISYQ,ISYJ,II,IJ,ERI,SCR)

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ISYP, ISYI, ISYQ, ISYJ, II, IJ
real(kind=wp), intent(_OUT_) :: ERI(*), SCR(*)
integer(kind=iwp) :: I, I3, I4, IDISK, IREC, IS12, IS34, ISY1, ISY2, ISY3, ISY4, NBUF, NDIM2M, NO1, NO2
logical(kind=iwp) :: TRANSP
#include "intgrl.fh"

NDIM2M = (NSYMZ*(NSYMZ+1))/2
if (ISYP >= ISYQ) then
  ISY1 = ISYP
  ISY2 = ISYQ
  ISY3 = ISYI
  ISY4 = ISYJ
  I3 = II
  I4 = IJ
  TRANSP = .true.
else
  ISY1 = ISYQ
  ISY2 = ISYP
  ISY3 = ISYJ
  ISY4 = ISYI
  I3 = IJ
  I4 = II
  TRANSP = .false.
end if

! We are going for an EXCH1 or EXCH2 integral of symmetry
! block ISY1,ISY3,ISY2,ISY4.
! Disk address to symmetry block, and record nr:

IS12 = (ISY1*(ISY1-1))/2+ISY2
if (ISY3 > ISY4) then
  IS34 = (ISY3*(ISY3-1))/2+ISY4
  IDISK = IAD2M(2,IS12+NDIM2M*(IS34-1))
  IREC = I4+NOSHZ(ISY4)*(I3-1)
else if (ISY3 == ISY4) then
  IS34 = (ISY3*(ISY3-1))/2+ISY4
  if (I3 >= I4) then
    IDISK = IAD2M(2,IS12+NDIM2M*(IS34-1))
    IREC = (I3*(I3-1))/2+I4
  else
    IDISK = IAD2M(2,IS12+NDIM2M*(IS34-1))
    IREC = (I4*(I4-1))/2+I3
    TRANSP = .not. TRANSP
  end if
else
  IS34 = (ISY4*(ISY4-1))/2+ISY3
  IDISK = IAD2M(3,IS12+NDIM2M*(IS34-1))
  IREC = I3+NOSHZ(ISY3)*(I4-1)
end if

! Buffer size:
NO1 = NORBZ(ISY1)
NO2 = NORBZ(ISY2)
NBUF = NO1*NO2
if (NBUF == 0) return

! Address update for earlier records, then read:
! PAM07 * Eliminate unsafe IPOSFILE call
!IDISK = IDISK+iPosFile(NBUF)*(IREC-1)
! Replace with dummy i/o operations: Is this efficience issue??
do I=1,IREC-1
  call dDAFILE(LUINTMZ,0,SCR,NBUF,IDISK)
end do

if (TRANSP) then
  call dDAFILE(LUINTMZ,2,SCR,NBUF,IDISK)
  call TRNSPS(NO2,NO1,SCR,ERI)
else
  call dDAFILE(LUINTMZ,2,ERI,NBUF,IDISK)
end if

return

end subroutine EXCH
