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
! Copyright (C) 1996, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine SPGRPCON(IOFSPGRP,NSPGRP,NGAS,MXPNGAS,IELFSPGRP,ISPGRPCON)
! Find connection matrix for string types
!
! ISPGRPCON(ISPGP,JSPGRP) = 0 => spgrps are identical
!                         = 1 => spgrps are connected by single excitation
!                         = 2 => spgrps are connected by double excitation
!                        >= 3 => spgrps are connected by triple or higher excitation
!
! Jeppe Olsen, September 1996

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: wp, u6
#endif

implicit none
integer(kind=iwp), intent(in) :: IOFSPGRP, NSPGRP, NGAS, MXPNGAS, IELFSPGRP(MXPNGAS,*)
integer(kind=iwp), intent(out) :: ISPGRPCON(NSPGRP,NSPGRP)
integer(kind=iwp) :: IDIF, ISPGRP, ISPGRPA, JSPGRP, JSPGRPA, NEXC
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: NEXC1, NEXC2
#endif

do ISPGRP=1,NSPGRP
  ISPGRPA = IOFSPGRP-1+ISPGRP
  do JSPGRP=1,ISPGRP
    JSPGRPA = IOFSPGRP-1+JSPGRP
    IDIF = sum(abs(IELFSPGRP(1:NGAS,ISPGRPA)-IELFSPGRP(1:NGAS,JSPGRPA)))
    NEXC = IDIF/2
    ISPGRPCON(ISPGRP,JSPGRP) = NEXC
    ISPGRPCON(JSPGRP,ISPGRP) = NEXC
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) '===================='
write(u6,*) 'output from SPGRPCON'
write(u6,*) '===================='
write(u6,*)
NEXC1 = 0
NEXC2 = 0
do ISPGRP=1,NSPGRP
  do JSPGRP=1,NSPGRP
    if (ISPGRPCON(ISPGRP,JSPGRP) == 1) then
      NEXC1 = NEXC1+1
    else if (ISPGRPCON(ISPGRP,JSPGRP) == 2) then
      NEXC2 = NEXC2+1
    end if
  end do
end do

write(u6,*) ' single excitation interactions',NEXC1,'( ',real(NEXC1,kind=wp)*100.0_wp/real(NSPGRP,kind=wp)**2,' % )'
write(u6,*) ' double excitation interactions',NEXC2,'( ',real(NEXC2,kind=wp)*100.0_wp/real(NSPGRP,kind=wp)**2,' % )'

write(u6,*) ' Supergroup connection matrix'
call IWRTMA(ISPGRPCON,NSPGRP,NSPGRP,NSPGRP,NSPGRP)
#endif

end subroutine SPGRPCON
