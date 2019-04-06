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
      SUBROUTINE ORBINF_MCLR(NIRREP,NSMOB,NRAS1,NRAS2,NRAS3,
     &                  MXR4tp,IPRNT)
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
      DIMENSION NRAS1(NIRREP),NRAS2(NIRREP),NRAS3(NIRREP)
      DIMENSION ITFSO(1)
*
#include "detdim.fh"
#include "orbinp_mclr.fh"
*
*
      NTEST = 0000
      NTEST = MAX(NTEST,IPRNT)
************************************************
*                                              *
* Part 1 : From shell format to orbital format *
*                                              *
************************************************
      DO 10 IPIRR = 1, MXPIRR
        NOSPIR(IPIRR) = 1
        IOSPIR(1,IPIRR) = IPIRR
10    CONTINUE
*
* 2 : Shell information to orbital information for each group of orbital
*
*. RAS1
      Call iCopy(3*MXPOBS,[0],0,NRSOBS,1)
      NORB1 = 0
      DO 20 IRREP = 1, NIRREP
        DO 30 ISM = 1, NOSPIR(IRREP)
          IISM = IOSPIR(ISM,IRREP)
          NRSOBS(IISM,1) = NRSOBS(IISM,1) + NRAS1(IRREP)
          NORB1 = NORB1 + NRAS1(IRREP)
30      CONTINUE
20    CONTINUE
*. RAS2

      NORB2 = 0
      DO 40 IRREP = 1, NIRREP
        DO 50 ISM = 1, NOSPIR(IRREP)
          IISM = IOSPIR(ISM,IRREP)
          NRSOBS(IISM,2) = NRSOBS(IISM,2) + NRAS2(IRREP)
          NORB2 = NORB2 + NRAS2(IRREP)
50      CONTINUE
40    CONTINUE

*. RAS3
      NORB3 = 0
      DO 60 IRREP = 1, NIRREP
        DO 70 ISM = 1, NOSPIR(IRREP)
          IISM = IOSPIR(ISM,IRREP)
          NRSOBS(IISM,3) = NRSOBS(IISM,3) + NRAS3(IRREP)
          NORB3 = NORB3 + NRAS3(IRREP)
70      CONTINUE
60    CONTINUE
*. Inactive, RAS0, RAS4, deleted
      NORB4=0
      NORB0=0
      NINOB=0
      NDEOB=0
      DO 80 IRREP = 1, NIRREP
        DO 90 ISM = 1, NOSPIR(IRREP)
          IISM = IOSPIR(ISM,IRREP)
          NINOBS(IISM) = 0
          NDEOBS(IISM) = 0
          NR0OBS(1,IISM) = 0
          DO 100 IPR4T = 1 , MXPR4T
             NR4OBS(IISM,IPR4T)=0
100       CONTINUE
90      CONTINUE
80    CONTINUE
*. Active, occupied  and total number of orbitals
      NACOB = NORB1+NORB2+NORB3
      NOCOB = NACOB+NINOB+NORB0+NORB4
      NTOOB = NOCOB + NDEOB
      DO 110 ISMOB = 1, NSMOB
         NACOBS(ISMOB)=NRSOBS(ISMOB,1)
     &                +NRSOBS(ISMOB,2)
     &                +NRSOBS(ISMOB,3)
         NOCOBS(ISMOB)=NACOBS(ISMOB)
         NTOOBS(ISMOB)=NACOBS(ISMOB)
110   CONTINUE
*
      IF(NTEST.GT.0) THEN
        WRITE(6,*)
        WRITE(6,*) ' ORBINF speaking'
        WRITE(6,*) ' ==============='
        WRITE(6,*) ' Number of orbitals per symmetry + total '
        WRITE(6,'(1H ,A,10I4,8X,I3)')
     *  '     Ras1             ',(NRSOBS(I,1),I=1,NSMOB),NORB1
        WRITE(6,'(1H ,A,10I4,8X,I3)')
     *  '     Ras2             ',(NRSOBS(I,2),I=1,NSMOB),NORB2
        WRITE(6,'(1H ,A,10I4,8X,I3)')
     *  '     Ras3             ',(NRSOBS(I,3),I=1,NSMOB),NORB3
        WRITE(6,'(1H ,A,10I4,8X,I3)')
     *  '     Active           ',(NACOBS(I),I=1,NSMOB),NACOB
        WRITE(6,'(1H ,A,10I4,8X,I3)')
     *  '     Total            ',(NTOOBS(I),I=1,NSMOB),NTOOB
      END IF
*. Offsets for orbitals of given symmetry
      ITOOBS(1) = 1
      DO 500 ISMOB = 2, NSMOB
        ITOOBS(ISMOB) = ITOOBS(ISMOB-1)+NTOOBS(ISMOB-1)
500   CONTINUE
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
      CALL ORBORD(NSMOB,MXPOBS,MXR4TP,NDEOBS,NINOBS,NR0OBS,NACOBS,
     *     NRSOBS,NR4OBS,NOCOBS,NTOOBS,IREOST,IREOTS,ISMFTO,ITFSO,IPRNT,
     *     IBSO,NTSOB,IBTSOB,ITSOB,NOBPTS,IOBPTS,MXPR4T,
     *     ISMFSO,ITPFTO,NOBPT)
*
*
      RETURN
      END
