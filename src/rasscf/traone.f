************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE TRAONE_fciqmc(PAO,PMO,TEMP,CMO)
*
*     Transformation program: one-electron section
*
*     Objective: transformes a one-electron matrix PAO in AO-basis
*                to a molecular orbital matrix PMO.
*
*     Subroutine calls: none
*
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
#include "fciqmc_global.fh"
#include "WrkSpc.fh"
*
      Real*8 CMO(*)
      DIMENSION PAO(*),PMO(*),TEMP(*)
*
      Call qEnter('TraOne')
*
      ICMO=1
      IAO =1
      IMO =1
      DO 100 ISYM=1,NSYM
       ICMO=ICMO+NBAS(ISYM)*NFRO(ISYM)
       IOFF=1+NBAS(ISYM)*NBAS(ISYM)
       if(NORB(ISYM).ne.0) then
         CALL SQUARE(PAO(IAO),TEMP(1),1,NBAS(ISYM),NBAS(ISYM))
         CALL DGEMM_('T','N',NORB(ISYM),NBAS(ISYM),
     *                NBAS(ISYM),1.0d0,CMO(ICMO),
     *                 NBAS(ISYM),TEMP,NBAS(ISYM),
     *                 0.0d0,TEMP(IOFF),NORB(ISYM))
         CALL MXMT(TEMP(IOFF),    1,NORB(ISYM),
     *             CMO(ICMO),     1,NBAS(ISYM),
     *             PMO(IMO),
     *             NORB(ISYM),NBAS(ISYM))
       end if
       ICMO=ICMO+NBAS(ISYM)*(NORB(ISYM)+NDEL(ISYM))
       IAO =IAO +NBAS(ISYM)*(NBAS(ISYM)+1)/2
       IMO =IMO +NORB(ISYM)*(NORB(ISYM)+1)/2
100   CONTINUE
*
      Call qExit('TraOne')
*
      RETURN
      END
