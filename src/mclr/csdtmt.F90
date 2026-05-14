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
! Copyright (C) 1984,1989-1993, Jeppe Olsen                            *
!***********************************************************************

subroutine CSDTMT(IDFTP,ICFTP,DTOC,PSSIGN)
! Construct list of prototype combinations in IDFTP
! Construct list of prototype CSF'S, in ICFTP
! Construct matrix expanding prototype CSF's in terms of
! prototype combinations in DTOC

use MCLR_Data, only: MINOP, MS2P, MULTSP, NCPCNT, NDPCNT, NTYP
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(inout) :: IDFTP(*), ICFTP(*)
real(kind=wp), intent(_OUT_) :: DTOC(*)
real(kind=wp), intent(in) :: PSSIGN
integer(kind=iwp) :: ICDCBS, ICSBS, IDTBS, IOPEN, ITP, MS2, MULTS, NNDET
integer(kind=iwp), allocatable :: SCR7(:)

MULTS = MULTSP
MS2 = MS2P

! Set up determinants and upper determinants

IDTBS = 0 ! dummy initialize
ICSBS = 0 ! dummy initialize
do ITP=1,NTYP
  IOPEN = MINOP+ITP-1
  if (ITP == 1) then
    IDTBS = 1
    ICSBS = 1
  else
    IDTBS = IDTBS+(IOPEN-1)*NDPCNT(ITP-1)
    ICSBS = ICSBS+(IOPEN-1)*NCPCNT(ITP-1)
  end if

  if (IOPEN /= 0) then
    call mma_allocate(SCR7,IOPEN+1,Label='SCR7')
    ! Proto type determinants and upper determinants
    if (MS2+1 == MULTS) then
      call SPNCOM_MCLR(scr7,IOPEN,MS2,NNDET,IDFTP(IDTBS),ICFTP(ICSBS),2,PSSIGN)
    else
      call SPNCOM_MCLR(scr7,IOPEN,MS2,NNDET,IDFTP(IDTBS),ICFTP(ICSBS),1,PSSIGN)
      call SPNCOM_MCLR(scr7,IOPEN,MULTS-1,NNDET,IDFTP(IDTBS),ICFTP(ICSBS),3,PSSIGN)
    end if
    call mma_deallocate(SCR7)
  end if
end do
! Matrix expressing csf's in terms of combinations
IDTBS = 1
ICSBS = 1
ICDCBS = 1
do ITP=1,NTYP
  IOPEN = MINOP+ITP-1
  if (ITP > 1) then
    IDTBS = IDTBS+(IOPEN-1)*NDPCNT(ITP-1)
    ICSBS = ICSBS+(IOPEN-1)*NCPCNT(ITP-1)
    ICDCBS = ICDCBS+NDPCNT(ITP-1)*NCPCNT(ITP-1)
  end if
  if (NDPCNT(ITP)*NCPCNT(ITP) == 0) cycle
  if (IOPEN == 0) then
    DTOC(ICDCBS) = One
  else
    call CSFDET_MCLR(IOPEN,IDFTP(IDTBS),NDPCNT(ITP),ICFTP(ICSBS),NCPCNT(ITP),DTOC(ICDCBS),PSSIGN)
  end if
end do

end subroutine CSDTMT
