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

subroutine CSDIAG_MCLR(CSFDIA,DETDIA,NCNFTP,NTYP,ICTSDT,NDTFTP,NCSFTP)
! obtain average CI diagonal elements and store in
! CSFDIA as CSF diagonal

use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: CSFDIA(*)
real(kind=wp), intent(in) :: DETDIA(*)
integer(kind=iwp), intent(in) :: NTYP, NCNFTP(NTYP), ICTSDT(*), NDTFTP(NTYP), NCSFTP(NTYP)
integer(kind=iwp) :: ICNF, ICSF, ICSOFF, IDET, IDTOFF, ITYP, JCNABS, JCNF, JDET
real(kind=wp) :: EAVER

ICSOFF = 0
IDTOFF = 0
JCNABS = 0
do ITYP=1,NTYP
  IDET = NDTFTP(ITYP)
  ICSF = NCSFTP(ITYP)
  ICNF = NCNFTP(ITYP)
  do JCNF=1,ICNF
    JCNABS = JCNABS+1
    EAVER = Zero
    do JDET=1,IDET
      EAVER = EAVER+DETDIA(abs(ICTSDT(IDTOFF+JDET)))
    end do
    if (IDET /= 0) EAVER = EAVER/real(IDET,kind=wp)
    CSFDIA(ICSOFF+1:ICSOFF+ICSF) = EAVER
    ICSOFF = ICSOFF+ICSF
    IDTOFF = IDTOFF+IDET
  end do
end do

end subroutine CSDIAG_MCLR
