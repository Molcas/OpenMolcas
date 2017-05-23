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
* Copyright (C) 1991, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE ORBINF(IPRNT)
*
* Obtain information about orbitals from shell information
*
* =====
* input
* =====
* Shell and symmetry information in /LUCINP/
*
* ======
* Output
* ======
* Orbital information in /ORBINP/
*
* Jeppe Olsen, Winter of 1991
*
COLD  INTEGER CITYP
#include "mxpdim.fh"
#include "lucinp.fh"
#include "cgas.fh"
*
#include "orbinp.fh"
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)
************************************************
*                                              *
* Part 1 : From shell format to orbital format *
*                                              *
************************************************
      CALL OSPIR(    NOSPIR,    IOSPIR,    PNTGRP,    NIRREP,    MXPIRR,
     &               MXPOBS,     IPRNT)
*
* 2 : Shell information to orbital information for each group of orbital
*
*
* ===============
*     GAS case
* ===============
*
        DO IGAS = 1, NGAS
*. Shell => orbitals for each GAS space
          CALL SHTOOB(NGSSH(1,IGAS),NIRREP,  MXPOBS,   NSMOB,  NOSPIR,
     &                   IOSPIR,NGSOB(1,IGAS),NGSOBT(IGAS))
        END DO
*
*  ========================================================
*. Number of inactive, active, occupied , deleted orbitals
*  ========================================================
*
*
* current inactive and deleted orbitals are not identified so
       IGSINA = 0
       IGSDEL = 0
*
       CALL ISETVC(NTOOBS,0,NSMOB)
       CALL ISETVC(NOCOBS,0,NSMOB)
       CALL ISETVC(NACOBS,0,NSMOB)
*
       NTOOB = 0
       NACOB = 0
       NOCOB = 0
       DO IGAS = 1, NGAS
*. Inactive orbitals
         IF(IGAS.EQ.IGSINA) THEN
           CALL ICOPVE(NGSOB(1,IGAS),NINOBS,NSMOB)
           NINOB = NGSOBT(IGAS)
         END IF
*. Deleted orbitals
         IF(IGAS.EQ.IGSDEL) THEN
           CALL ICOPVE(NGSOB(1,IGAS),NDEOBS,NSMOB)
           NDEOB = NGSOBT(IGAS)
         END IF
*. Add to total number of orbitals
         CALL IVCSUM(   NTOOBS,   NTOOBS,NGSOB(1,IGAS),    1,    1,
     &                   NSMOB)
         NTOOB = NTOOB + NGSOBT(IGAS)
*. Add to occupied orbitals
         IF(IGAS.NE.IGSDEL) THEN
           CALL IVCSUM(  NOCOBS,  NOCOBS,NGSOB(1,IGAS),  1,  1,
     &                    NSMOB)
           NOCOB = NOCOB + NGSOBT(IGAS)
         END IF
*. Add to active orbitals
         IF(IGAS.NE.IGSINA.AND.IGAS.NE.IGSDEL) THEN
           CALL IVCSUM(  NACOBS,  NACOBS,NGSOB(1,IGAS),  1,  1,
     &                    NSMOB)
           NACOB = NACOB + NGSOBT(IGAS)
         END IF
       END DO
* ===============================================
*. Well, report back
* ===============================================
        IF(NTEST.GT.0) THEN
          WRITE(6,*)
          WRITE(6,*) ' Number of orbitals per symmetry :'
          WRITE(6,*) ' ================================='
          WRITE(6,*)
          WRITE(6,'(1H ,A,10I4,A)')
     &    '            Symmetry  ',(I,I = 1,NSMOB)
          WRITE(6,'(1H ,A,2X,10A)')
     &    '           ========== ',('====',I = 1,NSMOB)
          DO IGAS = 1, NGAS
            WRITE(6,'(1H      ,A,I3,7X,A,10I4,8X,I3)')
     &      '   GAS',IGAS,'      ',(NGSOB(I,IGAS),I=1,NSMOB),
     &      NGSOBT(IGAS)
          END DO
*
          WRITE(6,*) ' Total number of orbitals ', NTOOB
          WRITE(6,*) ' Total number of occupied orbitals ', NOCOB
        END IF
*. Offsets for orbitals of given symmetry
        ITOOBS(1) = 1
        DO  ISMOB = 2, NSMOB
          ITOOBS(ISMOB) = ITOOBS(ISMOB-1)+NTOOBS(ISMOB-1)
        END DO
*
        IF(NTEST.GT.0) THEN
          WRITE(6,*) ' Offsets for orbital of given symmetry '
          CALL IWRTMA(ITOOBS,1,NSMOB,1,NSMOB)
        END IF

********************************************
*                                          *
* Part 2 : Reordering arrays for orbitals  *
*                                          *
********************************************
        CALL ORBORD_GAS(   NSMOB,  MXPOBS, MXPNGAS,    NGAS,   NGSOB,
     &                    NGSOBT,  NOCOBS,  NTOOBS,   NTOOB,  IREOST,
     &                    IREOTS,  ISMFTO,  ITPFSO,    IBSO,  NOBPTS,
     &                    IOBPTS,  ISMFSO, ITPFTO,    NOBPT,   IPRNT)
*
*. Largest number of orbitals of given sym and type
      MXTSOB = 0
      MXTOB = 0
      DO IOBTP = 1, NGAS
        LTOB = 0
        DO IOBSM = 1, NSMOB
         MXTSOB = MAX(MXTSOB,NOBPTS(IOBTP,IOBSM))
         LTOB = LTOB + NOBPTS(IOBTP,IOBSM)
        END DO
        MXTOB= MAX(LTOB,MXTOB)
      END DO
      IF (NTEST .GT. 0)
     &   WRITE(6,*) ' MXTSOB,MXTOB from ORBINF = ', MXTSOB,MXTOB
*
      RETURN
      END
