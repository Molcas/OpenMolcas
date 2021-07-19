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

subroutine CSDIAG_CI_UTIL(CSFDIA,DETDIA,NCNFTP,NTYP,ICTSDT,NDTFTP,NCSFTP,IPRINT)
! PURPOSE: OBTAIN AVERAGE CI DIAGONAL ELEMENTS AND STORE IN CSFDIA
!
! CALLING PARAMETERS.
! CSFDIA  : CI DIAGONAL IN SCF BASIS
! DETDIA  : CI DIAGONAL IN DETERMINENT BASIS
! NCNFTP  :
! NTYP    :
! ICTSDT  :
! NDTFTP  :
! NCSFTP  :

use Constants, only: Zero
use Definitions, only: wp, iwp,u6

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: CSFDIA(*) !IFG
real(kind=wp), intent(in) :: DETDIA(*) !IFG
integer(kind=iwp), intent(in) :: NTYP, NCNFTP(NTYP), ICTSDT(*), NDTFTP(NTYP), NCSFTP(NTYP), IPRINT !IFG
real(kind=wp) :: EAVER
integer(kind=iwp) :: ICNF, ICSF, ICSOFF, IDET, IDTOFF, ITYP, JCNABS, JCNF, JDET, NCSTOT, NDTTOT

ICSOFF = 1
IDTOFF = 1
JCNABS = 0
do ITYP=1,NTYP
  IDET = NDTFTP(ITYP)
  ICSF = NCSFTP(ITYP)
  ICNF = NCNFTP(ITYP)
  do JCNF=1,ICNF
    JCNABS = JCNABS+1
    EAVER = Zero
    do JDET=1,IDET
      EAVER = EAVER+DETDIA(abs(ICTSDT(IDTOFF-1+JDET)))
    end do
    if (IDET /= 0) EAVER = EAVER/real(IDET,kind=wp)
    call DCOPY_(ICSF,[EAVER],0,CSFDIA(ICSOFF),1)
    ICSOFF = ICSOFF+ICSF
    IDTOFF = IDTOFF+IDET
  end do
end do

if (IPRINT >= 40) then
  NCSTOT = ICSOFF-1
  NDTTOT = IDTOFF-1
  write(u6,*)
  write(u6,*) ' CIDIAGONAL IN DET BASIS'
  call WRTMAT(DETDIA,1,NDTTOT,1,NDTTOT)
  write(u6,*)
  write(u6,*) ' CIDIAGONAL IN CSF BASIS'
  call WRTMAT(CSFDIA,1,NCSTOT,1,NCSTOT)
end if

return

end subroutine CSDIAG_CI_UTIL
