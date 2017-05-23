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
      SUBROUTINE CNFORD_GAS( IOCCLS, NOCCLS,   ISYM, PSSIGN, IPRCSF,
     &                      ICONF_OCC,ICONF_REO,ICTSDT,IBLOCK,NBLOCK)
*
* Generate configurations in ICONF
*
* Generate determinants in configuration order and obtain
* sign array for switching between the two formats.
*
*
* It is assumed that CSFDIM has been called
*
* Jeppe Olsen Dec. 2001 from CNFORD
*
*
#include "implicit.fh"
#include "mxpdim.fh"
#include "spinfo_lucia.fh"
#include "glbbas.fh"
#include "WrkSpc.fh"
#include "orbinp.fh"
#include "cgas.fh"
*. Specific input
       INTEGER IOCCLS(NGAS,NOCCLS)
       INTEGER IBLOCK(8,NBLOCK)
*. Output
      INTEGER ICONF_OCC(*),ICONF_REO(*)
      DIMENSION ICTSDT(*)
*
      NTEST = 0
      NELEC = IELSUM(IOCCLS(1,1),NGAS)
C
      IDUM=1
      CALL GETMEM('ZSCR  ','ALLO','INTE',KLZSCR,(NOCOB+1)*(NELEC+1))
      CALL GETMEM('Z     ','ALLO','INTE',KLZ,NOCOB*NELEC*2)
      CALL GETMEM('OCMIN ','ALLO','INTE',KLOCMIN,NOCOB)
      CALL GETMEM('OCMAX ','ALLO','INTE',KLOCMAX,NOCOB)
*. Zero configuration reorder array using NCONF_ALL_SYM
      IZERO = 0
c     CALL ISETVC(WORK(KICONF_REO(ISYM)),IZERO,NCONF_ALL_SYM)
      CALL ISETVC(iWORK(KICONF_REO(ISYM)),IZERO,NCONF_tot)
*
* Generate configurations for all occupation classes
*
      IB_OCCLS = 1
      DO JOCCLS = 1, NOCCLS
*. Save offset to current occupation class
*
*        KIB_OCCLS seems no longer to be used. Therefore the call to
*        ITOR is commented out to avoid having unitialized arrays
*        floating around and in case of ITOR there's even written to
*        this undefined piece of memory.
*
*        / Jesper Wisborg Krogh, 2005-06-22
*
C_REMOVED         CALL ITOR(WORK(KIB_OCCLS(ISYM)),1,IB_OCCLS,JOCCLS)
*.Max and min arrays for strings
      CALL MXMNOC_OCCLS(iWORK(KLOCMIN),iWORK(KLOCMAX),
     &                  NGAS,NOBPT,IOCCLS(1,JOCCLS),MINOP,NTEST)
*. the arcweights
         CALL CONF_GRAPH(iWORK(KLOCMIN),iWORK(KLOCMAX),
     &                   NOCOB, NELEC,iWORK(KLZ),NCONF_P,iWORK(KLZSCR))
*
         IF(JOCCLS.EQ.1) THEN
            INITIALIZE_CONF_COUNTERS = 1
         ELSE
            INITIALIZE_CONF_COUNTERS = 0
         END IF
         IDOREO = 1
*. Lexical addressing for configurations of this type
         IB_OCCLS = IBCONF_ALL_SYM_FOR_OCCLS(JOCCLS)
*
         CALL GEN_CONF_FOR_OCCLS(IOCCLS(1,JOCCLS),
     &                           IB_OCCLS,
     &                           INITIALIZE_CONF_COUNTERS,
     &                             NGAS,  ISYM, MINOP, MAXOP, NSMST,
     &                                0,
*
     &                            NOCOB,
     &                            NOBPT,
     &                           NCONF_PER_OPEN(1,ISYM),
     &                           NCONF_OCCLS,
     &                           IB_CONF_REO,
*
     &                           IB_CONF_OCC,
     &                           iWORK(KICONF_OCC(ISYM)),
     &                           IDOREO,
     &                           iWORK(KLZ),
     &                           NCONF_ALL_SYM,
*
     &                           iwork(kiconf_reo(isym)),
     &                           nconf_tot)
c      write(6,*)'nconf_per_open after 2nd call of gen_conf_for_occls'
c      call iwrtma(nconf_per_open,1,4,1,4)
c     write(6,*)'nconf_tot in gen_conf_for_occls',nconf_tot
*
C     GEN_CONF_FOR_OCCLS(
C    &    IOCCLS,IB_OCCLS,INITIALIZE_CONF_COUNTERS,
C    &    NGAS,ISYM,MINOP,MAXOP,NSMST,IONLY_NCONF,NTORB,NOBPT,
C    &
C    NCONF_OP,IBCONF_REO,IBCONF_OCC,ICONF,IDOREO,IZ_CONF,IREO,NCONF_ALL_SYM)
CError  IB_OCCLS = IB_OCCLS + NCONF_ALL_SYM
      END DO
*
* Reorder Determinants from configuration order to ab-order
*
C          REO_GASDET(INBLOCK,NBLOCK,ISYM,IREO,SREO )
      CALL REO_GASDET(IBLOCK,NBLOCK,ISYM,ICTSDT)

C     CALL CNTOST(ICONF,SGNCTS,ICTSDT,IWORK,IDCNF(1),IREFSM,
C    &            IREFML,IREFG,XNDXCI,NORB,NEL,ORBSYM,
C    &            ICNFOK,IGENSG,ISGNA,ISGNB,IPRCSF)
*
      CALL GETMEM('ZSCR  ','FREE','INTE',KLZSCR,(NOCOB+1)*(NELEC+1))
      CALL GETMEM('Z     ','FREE','INTE',KLZ,NOCOB*NELEC*2)
      CALL GETMEM('OCMIN ','FREE','INTE',KLOCMIN,NOCOB)
      CALL GETMEM('OCMAX ','FREE','INTE',KLOCMAX,NOCOB)
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_real(PSSIGN)
        CALL Unused_integer(IPRCSF)
        CALL Unused_integer_array(ICONF_OCC)
        CALL Unused_integer_array(ICONF_REO)
      END IF
      END
