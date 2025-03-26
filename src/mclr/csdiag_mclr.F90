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

subroutine CSDIAG_MCLR(CSFDIA,DETDIA,NCNFTP,NTYP,ICTSDT,NDTFTP,NCSFTP,IFLAG,NCNFCN,ICNFOK)
! obtain average CI diagonal elements and store in
! CSFDIA as CSF diagonal

implicit real*8(A-H,O-Z)
dimension CSFDIA(*), DETDIA(*)
dimension NCNFTP(NTYP), NDTFTP(NTYP), NCSFTP(NTYP)
dimension ICTSDT(*), ICNFOK(*)

ICSOFF = 1
IDTOFF = 1
JCNABS = 0
do ITYP=1,NTYP
  IDET = NDTFTP(ITYP)
  ICSF = NCSFTP(ITYP)
  ICNF = NCNFTP(ITYP)
  do JCNF=1,ICNF
    JCNABS = JCNABS+1
    EAVER = 0.0d0
    do JDET=1,IDET
      EAVER = EAVER+DETDIA(abs(ICTSDT(IDTOFF-1+JDET)))
    end do
    if (IDET /= 0) EAVER = EAVER/dble(IDET)
    CSFDIA(ICSOFF:ICSOFF+ICSF-1) = EAVER
    ICSOFF = ICSOFF+ICSF
    IDTOFF = IDTOFF+IDET
  end do
end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(IFLAG)
  call Unused_integer(NCNFCN)
  call Unused_integer_array(ICNFOK)
end if

end subroutine CSDIAG_MCLR
