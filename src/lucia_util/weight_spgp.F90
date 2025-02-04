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

subroutine WEIGHT_SPGP(Z,NORBTP,NELFTP,NORBFTP,ISCR,NTEST)
! construct vertex weights for given supergroup
!
! Reverse lexical ordering is used

use Definitions, only: u6

implicit real*8(A-H,O-Z)
! Input
integer NELFTP(NORBTP), NORBFTP(NORBTP)
! Ouput
integer Z(*)
! Scratch length : 2 * NORB + (NEL+1)*(NORB+1)
integer ISCR(*)
integer, external :: IELSUM

NORB = IELSUM(NORBFTP,NORBTP)
NEL = IELSUM(NELFTP,NORBTP)

if (NTEST >= 100) then
  write(u6,*) ' Subroutine WEIGHT_SPGP in action'
  write(u6,*) ' ================================'
  write(u6,*) 'NELFTP'
  call IWRTMA(NELFTP,1,NORBTP,1,NORBTP)
end if

KLFREE = 1
KLMAX = KLFREE
KLFREE = KLFREE+NORB

KLMIN = KLFREE
KLFREE = KLFREE+NORB

KW = KLFREE
KLFREE = KW+(NEL+1)*(NORB+1)
! Max and min arrays for strings
call MXMNOC_SPGP(ISCR(KLMIN),ISCR(KLMAX),NORBTP,NORBFTP,NELFTP,NTEST)
! Arc weights
call GRAPW(ISCR(KW),Z,ISCR(KLMIN),ISCR(KLMAX),NORB,NEL,NTEST)

end subroutine WEIGHT_SPGP
