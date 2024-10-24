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
*               2024, Giovanni Li Manni                                *
************************************************************************
* G. Li Manni, June 2024: Scale-up capability for single SD ROHF type calculations
      SUBROUTINE GEN_CONF_FOR_OCCLS(     IOCCLS,
     &                              IB_OCCLS,
     &                              INITIALIZE_CONF_COUNTERS,
     &                                NGAS,  ISYM, MINOP, MAXOP, NSMST,
     &                              IONLY_NCONF,
*
     &                               NTORB,
     &                               NOBPT,
     &                              NCONF_OP,
     &                               NCONF,
     &                              IBCONF_REO,
*
     &                              IBCONF_OCC,
     &                               ICONF,
     &                              IDOREO,
     &                              IZ_CONF,
     &                              NCONF_ALL_SYM,
*
     &                                ireo,
     &                              nconf_tot)
*
* IONLY_NCONF = 1 :
*
* Generate number of configurations of occclass IOCCLS and sym ISYM
*
* IONLY_NCONF = 0 :
*
* Generate number and actual configurations of occclass IOCCLS
* and sym ISYM
*
*
* Jeppe Olsen, Nov. 2001
*
      Implicit REAL*8 (A-H,O-Z)
#include "mxpdim.fh"
*
*.. Input
*
*. Number of electrons per gas space
      INTEGER IOCCLS(NGAS)
*. Number of orbitals per gasspace
      INTEGER NOBPT(NGAS)
*. Arc weights for configurations
      INTEGER IZ_CONF(*)
*. Offset for reordering array and occupation array
      INTEGER IBCONF_REO(*), IBCONF_OCC(*)
*
*.. Output
*
*. Number of configurations per number of open shells, all symmetries
      INTEGER NCONF_OP(MAXOP+1)
*. And the actual configurations
      INTEGER ICONF(*)
*. Reorder array : Lex number => Actual number
c     INTEGER IREO(*)
      integer ireo(nconf_tot)
*. Local scratch
      INTEGER JCONF(2*MXPORB)

      INTEGER IDUM_ARR(1)
*
      NTEST = 00
*. Total number of electrons
       NEL = IELSUM(IOCCLS,NGAS)
       IF(INITIALIZE_CONF_COUNTERS.EQ.1) THEN
         IZERO = 0
         CALL ISETVC(NCONF_OP,IZERO,MAXOP+1)
         NCONF_ALL_SYM = 0
       END IF
*. Loop over configurations
       INI = 1
       NCONF = 0
       ISUM = 0
       CALL ISETVC(JCONF,IZERO,2*MXPORB)

 1000  CONTINUE
*  Generate an array of integers from 1 to NEL.
*  It is the only CSF that matter for HS calculations.
*  Skip any loop below... It is an overwhelmingly long loop.
         If(NEL.eq.MINOP.and.NEL.eq.NOBPT(2)) then
          do i = 1, NEL
            JCONF(i) = i
          end do
          NONEW = 0
         Else
           CALL NEXT_CONF_FOR_OCCLS (JCONF,IOCCLS,NGAS,NOBPT,INI,NONEW)
         end if
         ISUM = ISUM + 1
         INI = 0
         IF(NONEW.EQ.0) THEN
*. Check symmetry and number of open orbitals for this space
           ISYM_CONF = ISYMST(JCONF,NEL)
           NOPEN     = NOP_FOR_CONF(JCONF,NEL)
           NOCOB =  NOPEN + (NEL-NOPEN)/2
           IF(NOPEN.GE.MINOP .OR. IONLY_NCONF .NE. 0)
     &           NCONF_ALL_SYM = NCONF_ALL_SYM + 1
           IF(ISYM_CONF.EQ.ISYM.AND.NOPEN.GE.MINOP) THEN
*. A new configuration to be included, reform and save in packed form
             NCONF = NCONF + 1
             NCONF_OP(NOPEN+1) = NCONF_OP(NOPEN+1) + 1
             IF(IONLY_NCONF .EQ. 0 ) THEN
*. Lexical number of this configuration
               IB_OCC = IBCONF_OCC(NOPEN+1)
     &                + (NCONF_OP(NOPEN+1)-1)*NOCOB
               CALL REFORM_CONF_OCC(JCONF,ICONF(IB_OCC),NEL,NOCOB,1)
               IF(IDOREO.NE.0) THEN
c.. Giovanni and Dongxia 2011.1.31
                 ilexnum = ilex_for_conf_new(iconf(ib_occ),nocob,
     &                ntorb,nel,iz_conf,0,idum_arr,idum,idum)
                 JREO = IBCONF_REO(NOPEN+1) -1 + NCONF_OP(NOPEN+1)
c.. Giovanni and Dongxia 2011
                 ireo(jreo)=ib_occls-1+ilexnum
               END IF
             END IF
           END IF
           If(NEL.ne.MINOP.or.NEL.ne.NOBPT(2)) GOTO 1000
         END IF
*        ^ End if nonew = 0
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*)  ' ====================================== '
        WRITE(6,*)  ' Results from configuration generator : '
        WRITE(6,*)  ' ====================================== '
        WRITE(6,*)
        WRITE(6,*) ' Occupation class in action : '
        CALL IWRTMA(IOCCLS,1,NGAS,1,NGAS)
        WRITE(6,*) ' Number of configurations of correct symmetry ',
     &       NCONF
        WRITE(6,*) ' Number of configurations of all symmetries   ',
     &       NCONF_ALL_SYM
        WRITE(6,*)
     &  ' Number of configurations for various number of open orbs'
        CALL IWRTMA(NCONF_OP,1,MAXOP+1,1,MAXOP+1)
        IF(IONLY_NCONF.EQ.0) THEN
          WRITE(6,*)
     &   ' Updated list of configurations (may not be the final...)'
          CALL WRT_CONF_LIST(ICONF,NCONF_OP,MAXOP,NCONF,NEL)
          WRITE(6,*)
     &   ' Updated reordering of conf, Lex=>Act (may not be the final'
          CALL IWRTMA(IREO,1,NCONF_tot,1,NCONF_tot)
        END IF
      END IF
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(NSMST)
      END
