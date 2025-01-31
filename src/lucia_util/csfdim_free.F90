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
! Copyright (C) 1994, Jeppe Olsen                                      *
!               2024, Giovanni Li Manni                                *
!***********************************************************************

subroutine CSFDIM_FREE(ISYM)
! Free resources allocated by CSFDIM_GAS

use stdalloc, only: mma_deallocate
use GLBBAS, only: DFTP, CFTP, DTOC, SDREO_I, CONF_OCC, CONF_REO, Z_PTDT, REO_PTDT, SDREO
use lucia_data, only: MINOP, MAXOP

implicit none
integer ISYM
#include "warnings.h"
integer IOPEN, ITYP

do IOPEN=MINOP,MAXOP
  ITYP = IOPEN+1

  call mma_deallocate(Z_PTDT(ITYP)%I)
  call mma_deallocate(REO_PTDT(ITYP)%I)
end do
deallocate(Z_PTDT)
deallocate(REO_PTDT)

!LDET = NSD_PER_SYM(ISYM)
!LCONF = 0
!LLCONF = 0
!do IOPEN=0,MAXOP
!  ITYP = IOPEN+1
!  !FIXME: NELEC is undefined
!  ICL = (NELEC-IOPEN)/2
!  LLCONF = LLCONF+NCONF_PER_OPEN(ITYP,ISYM)*(IOPEN+ICL)
!end do
!LCONF = max(LCONF,LLCONF)

call mma_deallocate(DFTP)
call mma_deallocate(CFTP)
call mma_deallocate(DTOC)

call mma_deallocate(CONF_OCC(ISYM)%I)
call mma_deallocate(CONF_REO(ISYM)%I)

call mma_deallocate(SDREO_I(ISYM)%I)
nullify(SDREO)

end subroutine CSFDIM_FREE
