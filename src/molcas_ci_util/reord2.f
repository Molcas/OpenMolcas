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
* Copyright (C) 1990,1996, Markus P. Fuelscher                         *
*               1990, Jeppe Olsen                                      *
************************************************************************
      Subroutine Reord2(NORB,NEL,IREFSM,IMODE,
     &                  ICONF,ISPIN,CIOLD,CINEW,kcnf)
************************************************************************
*                                                                      *
*     Rearrange CI-vectors                                             *
*     iMode=0 --> from SGA to split graph GUGA order                   *
*     iMode=1 --> from split graph GUGA to SGA order                   *
*                                                                      *
*     calling arguments:                                               *
*     nOrb    : integer                                                *
*               total number of active orbitals                        *
*     nEl     : integer                                                *
*               total number of active electrons                       *
*     iRefSm  : integer                                                *
*               state symmetry                                         *
*     iMode   : integer                                                *
*               switch selecting reordering mode (see above)           *
*     iConf   : array of integer                                       *
*               string information                                     *
*     iSpin   : array of integer                                       *
*               spin coupling information                              *
*     nSm     : array of integer                                       *
*               symmetry per active orbital                            *
*     CIold   : array of real*8                                        *
*               incoming CI vector                                     *
*     CInew   : array of real*8                                        *
*               outgoing CI vector                                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher and J. Olsen                                      *
*     University of Lund, Sweden, 1990                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*     - updated for integral direct and reaction field calculations    *
*       M.P. Fuelscher, University of Lund, Sweden, 1996               *
*                                                                      *
************************************************************************

      Implicit Real*8 (A-H,O-Z)


#include "rasdim.fh"
#include "strnum.fh"
#include "ciinfo.fh"
#include "spinfo.fh"
#include "gugx.fh"
#include "WrkSpc.fh"
#include "output_ras.fh"
      PARAMETER (ROUTINE='REORD2  ')
C
      DIMENSION ICONF(*),ISPIN(*)
      DIMENSION CIOLD(*),CINEW(*)
      DIMENSION KCNF(NEL)
C
      DIMENSION IWALK(mxAct)

      IPRLEV=IPRLOC(3)
C
C     LOOP OVER CONFIGURATIONS TYPES
C
      ICSFJP = 0
      ICNBS0 = 0
      IPBAS = 0
      DO 1000 ITYP = 1, NTYP
        IOPEN = ITYP + MINOP - 1
        ICL = (NEL - IOPEN) / 2
C      BASE ADRESS FOR CONFIGURATION OF THIS TYPE
        IF( ITYP .EQ. 1 ) THEN
          ICNBS0 = 1
        ELSE
          ICNBS0 = ICNBS0 + NCNFTP(ITYP-1,IREFSM)*(NEL+IOPEN-1)/2
        END IF
C      BASE ADRESS FOR PROTOTYPE SPIN COUPLINGS
        IF( ITYP .EQ. 1 ) THEN
          IPBAS = 1
        ELSE
          IPBAS = IPBAS + NCSFTP(ITYP-1)*(IOPEN-1)
        END IF
C
C     LOOP OVER NUMBER OF CONFIGURATIONS OF TYPE ITYP AND PROTOTYPE
C     SPIN COUPLINGS
C
        DO 900  IC = 1, NCNFTP(ITYP,IREFSM)
          ICNBS = ICNBS0 + (IC-1)*(IOPEN+ICL)
          DO 800 IICSF = 1,NCSFTP(ITYP)
            ICSFJP = ICSFJP + 1
            ICSBAS = IPBAS + (IICSF-1)*IOPEN
*. Obtain configuration in standard RASSCF form
            IIBOP = 1
            IIBCL = 1
            JOCC  = ICL + IOPEN
            DO KOCC = 0, JOCC-1
              KORB = ICONF(ICNBS+KOCC)
              IF(KORB.LT.0) THEN
*. Doubly occupied orbitals
                KCNF(IIBCL) = ABS(KORB)
                IIBCL = IIBCL + 1
              ELSE
*. Singly occupied orbital
                KCNF(ICL+IIBOP) = KORB
                IIBOP = IIBOP + 1
              END IF
            END DO
*
C     COMPUTE STEP VECTOR
            CALL STEPVEC(KCNF(1),KCNF(ICL+1),ICL,IOPEN,
     &                   ISPIN(ICSBAS),NORB,IWALK)
C     GET SPLIT GRAPH ORDERING NUMBER
            ISG=ISGNUM(IWORK(LDOWN),IWORK(LUP),
     &                 IWORK(LDAW),IWORK(LRAW),
     &                 IWORK(LUSGN),IWORK(LLSGN),
     &                 IWALK)
C     GET PHASE PHASE FACTOR
            IP=IPHASE(IWORK(LDRT),IWORK(LUP),IWALK)
            IF( IMODE.EQ.0 ) THEN
              CINEW(ISG)=CIOLD(ICSFJP)
              IF( IP.LT.0 ) CINEW(ISG)=-CIOLD(ICSFJP)
            ELSE
              CINEW(ICSFJP)=CIOLD(ISG)
              IF( IP.LT.0 ) CINEW(ICSFJP)=-CIOLD(ISG)
            ENDIF
800       CONTINUE
900     CONTINUE
1000  CONTINUE
C
      IF( IPRLEV.GE.DEBUG ) THEN
        LPRINT=MIN(200,ICSFJP)
        WRITE (6,*)
        WRITE (6,*) ' OLD CI-VECTOR IN SUBROUTINE REORD',
     &              ' (MAX. 200 ELEMENTS)'
        WRITE (6,'(10F12.8)') (CIOLD(I),I=1,LPRINT)
        WRITE (6,*) ' NEW CI-VECTOR IN SUBROUTINE REORD',
     &              ' (MAX. 200 ELEMENTS)'
        WRITE (6,'(10F12.8)') (CINEW(I),I=1,LPRINT)
        WRITE (6,*)
      ENDIF
C
      RETURN
      END
