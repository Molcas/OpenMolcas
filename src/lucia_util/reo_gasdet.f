************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2001, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE REO_GASDET(IBLOCK,NBLOCK,ISYM,IREO)
      use stdalloc, only: mma_allocate, mma_deallocate
      use GLBBAS
*
* Create reorder array for determinants : configuration order => Ab order
*
*
* Jeppe Olsen, November 2001, from GASANA
*
*
#include "implicit.fh"
*
* =====
*.Input
* =====
*
#include "mxpdim.fh"
#include "orbinp.fh"
#include "strbas.fh"
#include "cicisp.fh"
#include "cstate.fh"
#include "strinp.fh"
#include "stinf.fh"
#include "csm.fh"
#include "WrkSpc.fh"
#include "cgas.fh"
#include "gasstr.fh"
#include "cprnt.fh"
#include "spinfo_lucia.fh"
*
      DIMENSION IBLOCK(8,NBLOCK)
*
*. Output
      INTEGER IREO(*)
      Integer, Allocatable:: LASTR(:), LBSTR(:)
      Integer, Allocatable:: ZSCR(:), Z(:)

c      write(6,*)'nconf_per_open in reo_gasdet'
c      call iwrtma(nconf_per_open,1,4,1,4)
*
** Specifications of internal space
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRDIA)
* Type of alpha and beta strings
      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NEL= NAEL+NBEL
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)

*
*
**. Info on block structure of space
*
*
*Space for alpha and beta strings
      Call mma_allocate(LASTR,MXNSTR*NAEL,Label='LASTR')
      Call mma_allocate(LBSTR,MXNSTR*NBEL,Label='LBSTR')
*. Space for constructing arc weights for configurations
      CALL mma_allocate(ZSCR,(NOCOB+1)*(NEL+1),Label='ZSCR')
      CALL mma_allocate(Z,NOCOB*NEL*2,Label='Z')
      CALL GETMEM('OCMIN ','ALLO','INTE',KLOCMIN,NOCOB)
      CALL GETMEM('OCMAX ','ALLO','INTE',KLOCMAX,NOCOB)
*. Occupation and projections of a given determinant
C??? Jesper      CALL MEMMAN(KLDET_OC,NOCOB,'ADDL  ',1,'CONF_O')
      CALL GETMEM('CONF_O','ALLO','INTE',KLDET_OC,NAEL+NBEL)
      CALL GETMEM('CONF_M','ALLO','INTE',KLDET_MS,NAEL+NBEL)
C??? Jesper      CALL MEMMAN(KLDET_MS,NOCOB,'ADDL  ',1,'CONF_M')
      CALL GETMEM('CONF_M','ALLO','INTE',KLDET_VC,NOCOB)

*
*     KIB_OCCLS removed from the argument list, since it's not defined
*
*     / Jesper Wisborg Krogh, 2005-06-22
      CALL REO_GASDET_S(    IREO,
     &                  IWORK(KNSTSO(IATP)),
     &                 IWORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,MXPNGAS,IOCTPA,
     &                    IOCTPB,  NBLOCK,  IBLOCK,    NAEL,    NBEL,
     &                 LASTR,LBSTR,NSMST,NELFSPGP,
*
     &                  NMXOCCLS,    NGAS,IWORK(KIOCLS),NTOOB,  NOBPT,
     &                  IWORK(KDFTP),
     &                  IB_CONF_REO,
     &                  iwork(kiconf_reo(isym)),
     &                  nconf_tot,
*
     &                  ib_conf_reo,
     &                     maxop,
     &                  nconf_per_open(1,isym),
     &                  IB_SD_FOR_OPEN,
     &                  ZSCR,
*
     &                  Z(:),
     &                  IWORK(KLOCMIN),
     &                  IWORK(KLOCMAX),
     &                  IWORK(KLDET_OC),
     &                  IWORK(KLDET_MS),
*
     &                  IWORK(KLDET_VC),iWORK,KZ_PTDT,KREO_PTDT, MINOP,
     &                  IBCONF_ALL_SYM_FOR_OCCLS,PSSIGN,NPDTCNF)
*
      Call mma_deallocate(LASTR)
      Call mma_deallocate(LBSTR)
      Call mma_deallocate(ZSCR)
      Call mma_deallocate(Z)
      CALL GETMEM('OCMIN ','FREE','INTE',KLOCMIN,NOCOB)
      CALL GETMEM('OCMAX ','FREE','INTE',KLOCMAX,NOCOB)
      CALL GETMEM('CONF_O','FREE','INTE',KLDET_OC,NAEL+NBEL)
      CALL GETMEM('CONF_M','FREE','INTE',KLDET_MS,NAEL+NBEL)
      CALL GETMEM('CONF_M','FREE','INTE',KLDET_VC,NOCOB)
*
      RETURN
      END
