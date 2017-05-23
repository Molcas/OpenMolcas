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
* Copyright (C) 1989,2003, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE CNHCN(ICNL,ITPL,ICNR,ITPR,CNHCNM,SCR,
     *                 NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,
     *                 NORB,TUVX,IPREXH,ExFac,IREOTS)
*
* Obtain Hamilton matrix over CSF's of configurations ICNL,ICNR
*
* Jeppe Olsen , Summer of '89
*               IREOTS added August 2003
*
      IMPLICIT REAL*8(A-H,O-Z)
C
#include "spinfo.fh"
#include "ciinfo.fh"
C
      DIMENSION ICNL(*),ICNR(*),SCR(*)
      DIMENSION IPRODT(*),DTOC(*),ONEBOD(*)
      DIMENSION CNHCNM(*)
      DIMENSION TUVX(*)
      INTEGER IREOTS(*)
C
      NTEST = 0
*
** 1 : Obtain determints corresponding to configurations
*
      NDETL = NDTFTP(ITPL)
      NDETR = NDTFTP(ITPR)
*
      NCSFL = NCSFTP(ITPL)
      NCSFR = NCSFTP(ITPR)
*
      KLFREE = 1
*
      KLDTLA = KLFREE
      KLFREE = KLFREE + NDETL*NAEL
*
      KLDTLB = KLFREE
      KLFREE = KLFREE + NDETL*NBEL
*
      KLISL  = KLFREE
      KLFREE = KLFREE + NDETL
*
      KLDTRA = KLFREE
      KLFREE = KLFREE + NDETR*NAEL
*
      KLDTRB = KLFREE
      KLFREE = KLFREE + NDETR*NBEL
*
      KLISR  = KLFREE
      KLFREE = KLFREE + NDETR
*
      CALL CNFSTR(ICNL,ITPL,SCR(KLDTLA),SCR(KLDTLB),
     *            NORB,NAEL,NBEL,NDETL,IPRODT,SCR(KLFREE),
     *            SCR(KLISL),IPREXH)
      CALL CNFSTR(ICNR,ITPR,SCR(KLDTRA),SCR(KLDTRB),
     *            NORB,NAEL,NBEL,NDETR,IPRODT,SCR(KLFREE),
     *            SCR(KLISR),IPREXH)
*
** Hamiltonin matrix over determinants of ths configuration
*
      KLDHD = KLFREE
      KLFREE = KLFREE + NDETL*NDETR
      CALL COMBINATIONS(ICOMBI,PSIGN)
      CALL DIHDJ_MOLCAS(SCR(KLDTLA),SCR(KLDTLB),NDETL,
     *           SCR(KLDTRA),SCR(KLDTRB),NDETR,
     *           NAEL,NBEL,SCR(KLFREE),NORB,ONEBOD,
     *           SCR(KLDHD),0,0,ECORE,ICOMBI,PSIGN,0,TUVX,ExFac,IREOTS)
*
** Transform matrix to CSF basis
*
** : sign changes
      CALL DGMM2_MOLCAS (SCR(KLDHD),SCR(KLDHD),SCR(KLISL),1,NDETL,NDETR)
      CALL DGMM2_MOLCAS (SCR(KLDHD),SCR(KLDHD),SCR(KLISR),2,NDETL,NDETR)
      IPL = 1
      DO 200 JTYP = 1, ITPL - 1
        JCSF = NCSFTP(JTYP)
        JDET = NDTFTP(JTYP)
        IPL = IPL + JCSF*JDET
200   CONTINUE
      IF(ITPR.EQ. ITPL) THEN
        IPR = IPL
      ELSE
        IPR = 1
        DO 300 JTYP = 1, ITPR - 1
          JCSF = NCSFTP(JTYP)
          JDET = NDTFTP(JTYP)
          IPR = IPR + JCSF*JDET
300     CONTINUE
      END IF
*
      KLCHD = KLFREE
      KLFREE = KLFREE + NCSFL*NDETR
      CALL MATML4(SCR(KLCHD),DTOC(IPL),SCR(KLDHD),
     *            NCSFL,NDETR,NDETL,NCSFL,NDETL,NDETR,1)
      CALL MATML4(CNHCNM,SCR(KLCHD),DTOC(IPR),
     *            NCSFL,NCSFR,NCSFL,NDETR,NDETR,NCSFR,0)
*
      IF( NTEST.GE.20 ) THEN
        WRITE(6,*)
     *  ' CSF-Hamiltonian matrix between two configurations'
        CALL WRTMAT(CNHCNM,NCSFL,NCSFR,NCSFL,NCSFR)
      END IF
*
      RETURN
      END
