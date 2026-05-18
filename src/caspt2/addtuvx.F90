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

subroutine ADDTUVX(NP,NI,NQ,NK,NASHT,iOffP,iOffI,iOffQ,iOffK,TUVX,nTUVX,PIQK,nPIQK,NUMERR)

use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#else
use Constants, only: One
#endif

implicit none
integer(kind=iwp), intent(in) :: NP, NI, NQ, NK, NASHT, iOffP, iOffI, iOffQ, iOffK, nTUVX, nPIQK
real(kind=wp), intent(inout) :: TUVX(nTUVX)
real(kind=wp), intent(in) :: PIQK(nPIQK)
integer(kind=iwp), intent(inout) :: NUMERR
integer(kind=iwp) :: iU, iUVX1, iUVX2, iV, iVX1, iVX2, iX, iX1, iX2
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iPIQK, iT, iTUVX
#include "warnings.h"
#else
#include "macros.fh"
unused_var(NUMERR)
#endif

! Add into correct positions in TUVX:

do iX=0,NK-1
  iX1 = NASHT*(iX+iOffK)
  iX2 = NQ*iX
  do iV=0,NQ-1
    iVX1 = NASHT*(iV+iOffQ+iX1)
    iVX2 = NI*(iV+iX2)
    do iU=0,NI-1
      iUVX1 = NASHT*(iU+iOffI+iVX1)
      iUVX2 = NP*(iU+iVX2)
#     ifdef _DEBUGPRINT_
      do iT=1,NP
        iTUVX = iT+iOffP+iUVX1
        iPIQK = iT+iUVX2
        ! Temporary test statements -- remove after debug!
        if ((ITUVX < 1) .or. (ITUVX > NTUVX+1)) then
          ITUVX = NTUVX+1
          NUMERR = NUMERR+1
          if (NUMERR > 100) then
            write(u6,*) ' THIS IS TOO MUCH -- STOP.'
            call QUIT(_RC_INTERNAL_ERROR_)
          end if
        end if
        ! End of temporary test statements
        TUVX(iTUVX) = TUVX(iTUVX)+PIQK(iPIQK)
      end do
#     else
      call DaXpY_(nP,One,PIQK(1+iUVX2),1,TUVX(1+iOffP+iUVX1),1)
#     endif
    end do
  end do
end do

end subroutine ADDTUVX
