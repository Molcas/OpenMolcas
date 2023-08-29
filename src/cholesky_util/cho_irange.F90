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

function CHO_IRANGE(I,BIN,NBIN,LOWER)
!
! Purpose: return range to which I belongs.
!          It is assumed that BIN is sorted
!          in ascending order so that:
!
!            I > BIN(NBIN)   ==> CHO_IRANGE = NBIN
!          Else
!            I > BIN(NBIN-1) ==> CHO_IRANGE = NBIN-1
!          Else
!            I > BIN(NBIN-2) ==> CHO_IRANGE = NBIN-2
!            ..
!            ..
!            ..
!          Else
!            I < BIN(1)      ==> CHO_IRANGE = 0
!
!          If NBIN <= 0      ==> CHO_IRANGE = -1
!
!          In case of degeneracies, the smallest bin is
!          returned if flag LOWER is set, i.e.,
!
!          I > BIN(IBIN) = BIN(IBIN-1)
!                        = ...
!                        = BIN(IBIN-N) ==> CHO_IRANGE = IBIN-N
!
!          else the largest bin (IBIN in the case above).

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: CHO_IRANGE
integer(kind=iwp), intent(in) :: I, BIN(*), NBIN
logical(kind=iwp), intent(in) :: LOWER
integer(kind=iwp) :: IBIN, JBIN, JBIN1

CHO_IRANGE = 0

if (NBIN < 1) then
  CHO_IRANGE = -1
else
  if (LOWER) then
    do IBIN=NBIN,1,-1
      if (I > BIN(IBIN)) then
        CHO_IRANGE = IBIN
        JBIN1 = IBIN-1
        do JBIN=JBIN1,1,-1
          if (BIN(JBIN) /= BIN(IBIN)) exit
          CHO_IRANGE = JBIN
        end do
        exit
      end if
    end do
  else
    CHO_IRANGE = 1
    do IBIN=NBIN,2,-1
      if (I > BIN(IBIN)) then
        CHO_IRANGE = IBIN
        exit
      end if
    end do
  end if
end if

end function CHO_IRANGE
