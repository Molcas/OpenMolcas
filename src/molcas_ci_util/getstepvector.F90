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

subroutine GETSTEPVECTOR(NOW,IOW,MV,IDWN,IUP,ICS)

use Definitions, only: iwp

implicit none
#include "rasdim.fh"
#include "general.fh"
#include "output_ras.fh"
#include "gugx.fh"
#include "WrkSpc.fh"
integer(kind=iwp), intent(in) :: NOW(2,NSYM,NMIDV), IOW(2,NSYM,NMIDV)
integer(kind=iwp), intent(inout) :: MV, IDWN, IUP
integer(kind=iwp), intent(out) :: ICS(NLEV)
integer(kind=iwp) :: IC1, ICDPOS, ICDWN, ICUP, ICUPOS, IDW0, IUW0, LEV, NDWN, NNN

! RECONSTRUCT THE CASE LIST

NICASE = NWALK*NIPWLK
!NSCR = 3*(NLEV+1)
!call mma_allocate(SCR,NSCR,label='SCR1')
!call mma_allocate(ICASE,NICASE,label='CASE')
!call MKCLIST(NSM,IWORK(LDOWN),IWORK(LNOW),IWORK(LIOW),ICASE,LSCR)
!call mma_deallocate(SCR)

! ENTER THE MAIN LOOP IS OVER BLOCKS OF THE ARRAY CI
! WITH SPECIFIED MIDVERTEX MV, AND UPPERWALK SYMMETRY ISYUP.
!
!do MV=1,NMIDV
!  do ISYUP=1,NSYM
NUP = NOW(1,1,MV)
NDWN = NOW(2,1,MV)
IUW0 = LICASE-NIPWLK+IOW(1,1,MV)
IDW0 = LICASE-NIPWLK+IOW(2,1,MV)
!    IDWNSV = 0
!    do IDWN=1,NDWN
!      do IUP=1,NUP
! determine the stepvector
!        if (IDWNSV /= IDWN) then
ICDPOS = IDW0+IDWN*NIPWLK
ICDWN = IWORK(ICDPOS)
! unpack lower walk
NNN = 0
do LEV=1,MIDLEV
  NNN = NNN+1
  if (NNN == 16) then
    NNN = 1
    ICDPOS = ICDPOS+1
    ICDWN = IWORK(ICDPOS)
  end if
  IC1 = ICDWN/4
  ICS(LEV) = ICDWN-4*IC1
  ICDWN = IC1
end do
!          IDWNSV = IDWN
!        end if
ICUPOS = IUW0+NIPWLK*IUP
ICUP = IWORK(ICUPOS)
! unpack upper walk
NNN = 0
do LEV=MIDLEV+1,NLEV
  NNN = NNN+1
  if (NNN == 16) then
    NNN = 1
    ICUPOS = ICUPOS+1
    ICUP = IWORK(ICUPOS)
  end if
  IC1 = ICUP/4
  ICS(LEV) = ICUP-4*IC1
  ICUP = IC1
end do
!      end do
!    end do
!  end do
!end do

! compute the next set of indices
if (IUP == NUP) then
  if (IDWN == NDWN) then
    if (MV == NMIDV) then
      MV = 0
    else
      MV = MV+1
    end if
    IDWN = 1
  else
    IDWN = IDWN+1
  end if
  IUP = 1
else
  IUP = IUP+1
end if

!call mma_deallocate(ICASE)
return

end subroutine GETSTEPVECTOR
