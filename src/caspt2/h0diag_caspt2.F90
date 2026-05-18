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

subroutine H0DIAG_CASPT2(ISYCI,DIAG,nDiag,NOW,IOW,nMidV)
! PURPOSE: FORM AN ARRAY OF DIAGONAL HAMILTONIAN MATRIX ELEMENTS
! FOR THE SPECIFIED TOTAL SYMMETRY ISYCI

use Symmetry_Info, only: Mul
use sguga, only: CIS
use caspt2_module, only: nSym
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ISYCI, nDiag, nMidV, NOW(2,NSYM,NMIDV), IOW(2,NSYM,NMIDV)
real(kind=wp), intent(out) :: DIAG(nDiag)
integer(kind=iwp) :: ICS, IEMU, ISYDWN, ISYUP, JCS, MV, NC, NDWN, NUP

DIAG(:) = Zero
IEMU = 1
do MV=1,NMIDV
  do ISYUP=1,NSYM
    NUP = NOW(1,ISYUP,MV)
    if (NUP == 0) cycle
    ISYDWN = Mul(ISYUP,ISYCI)
    NDWN = NOW(2,ISYDWN,MV)
    if (NDWN == 0) cycle
    ICS = 1+IOW(1,ISYUP,MV)
    JCS = 1+IOW(2,ISYDWN,MV)
    NC = NUP*NDWN
    call DIELMV(CIS%ICASE(ICS:),size(CIS%ICASE(ICS:)),CIS%ICASE(JCS:),size(CIS%ICASE(JCS:)),NUP,NDWN,DIAG(IEMU))
    IEMU = IEMU+NC
  end do
end do

end subroutine H0DIAG_CASPT2
