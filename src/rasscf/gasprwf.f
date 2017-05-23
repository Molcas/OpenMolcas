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
* Copyright (C) 1996, Markus P. Fuelscher                              *
************************************************************************
      Subroutine gasprwf(isel,NORB,NEL,IREFSM,
     &                  ICONF,ISPIN,CICOEF,kcnf)
************************************************************************
*                                                                      *
*     PURPOSE: PRINT THE WAVEFUNCTION FOR GAS                          *
*                                                                      *
************************************************************************
*                                                                      *
*     calling arguments:                                               *
*     nOrb    : integer                                                *
*               total number of active orbitals                        *
*     nEl     : integer                                                *
*               total number of active electrons                       *
*     iRefSm  : integer                                                *
*               state symmetry                                         *
*     iConf   : array of integer                                       *
*               string information                                     *
*     iSpin   : array of integer                                       *
*               spin coupling information                              *
*     nSm     : array of integer                                       *
*               symmetry per active orbital                            *
*     CiCoef   : array of real*8                                       *
*               incoming CI vector                                     *
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
#include "rasscf.fh"
#include "strnum.fh"
#include "ciinfo.fh"
#include "spinfo.fh"
#include "gugx.fh"
#include "gas.fh"
#include "WrkSpc.fh"
#include "output_ras.fh"
      PARAMETER (ROUTINE='GASPRWF ')
C
      DIMENSION ICONF(*),ISPIN(*)
      DIMENSION CICOEF(*)
      DIMENSION KCNF(NEL)
C
      DIMENSION IWALK(mxAct)
      integer   iSel(*)
      character*120 Line
C
C     print headline
C
      Line(1:16)='      Conf/sym  '
      iOff=16
      iSym=nSm(1)
      Do iorb=1,norb
         if ( nsm(iorb).ne.isym) ioff=ioff+1
         write(line(ioff+iorb:),'(I1)') nsm(iorb)
         if (nsm(iorb).ne.isym) isym=nsm(iorb)
      End do
      iOff=iOff+norb+3
      Line(iOff:iOff+15)='   Coeff Weight'
      Write(LF,'(A)') Line(1:iOff+15)
      Line=' '
C
C     Loop over configuration types
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
C     SKIP IT OR PRINT IT?
            COEF=CICOEF(ICSFJP)
            IF(ABS(COEF).LT.PRWTHR) GOTO 800
C     PRINT IT
              iSel(icsfjp) = 1
              Write(Line(1:),'(I8)') icsfjp
              iOff=10
              iSym=nSm(1)
              Do iorb=1,norb
                 If ( nSm(iorb).ne.iSym ) iOff=iOff+1
                 If ( iwalk(iorb).eq.3 ) then
                    Write(Line(iOff+iorb:),'(A1)') '2'
                 Else If (iwalk(iorb).eq.2) then
                    Write(Line(iOff+iorb:),'(A1)') 'd'
                 Else If (iwalk(iorb).eq.1) then
                    Write(Line(iOff+iorb:),'(A1)') 'u'
                 Else If (iwalk(iorb).eq.0) then
                    Write(Line(iOff+iorb:),'(A1)') '0'
                 End If
                 If ( nSm(iorb).ne.iSym ) iSym=nSm(iorb)
              End Do
              iOff=iOff+norb+3
              Write(Line(iOff:),'(2F8.5)') COEF,COEF**2
              Write(LF,'(6X,A)') Line(1:iOff+15)
              Line=' '

800       CONTINUE
900     CONTINUE
1000  CONTINUE

      goto 2000
2000  continue
      RETURN
      END
