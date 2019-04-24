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
      CALL CNHCN_INTERNAL(SCR)
*
*     This is to allow type punning without an explicit interface
      CONTAINS
      SUBROUTINE CNHCN_INTERNAL(SCR)
      USE ISO_C_BINDING
      REAL*8, TARGET :: SCR(*)
      INTEGER, POINTER :: iSCRla(:),iSCRlb(:),iSCRra(:),iSCRrb(:),
     &                    iSCRf(:)
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
      CALL C_F_POINTER(C_LOC(SCR(KLDTLA)),iSCRla,[1])
      CALL C_F_POINTER(C_LOC(SCR(KLDTLB)),iSCRlb,[1])
      CALL C_F_POINTER(C_LOC(SCR(KLFREE)),iSCRf,[1])
      CALL CNFSTR(ICNL,ITPL,iSCRla,iSCRlb,
     *            NORB,NAEL,NBEL,NDETL,IPRODT,iSCRf,
     *            SCR(KLISL),IPREXH)
      NULLIFY(iSCRla,iSCRlb,iSCRf)
      CALL C_F_POINTER(C_LOC(SCR(KLDTRA)),iSCRra,[1])
      CALL C_F_POINTER(C_LOC(SCR(KLDTRB)),iSCRrb,[1])
      CALL C_F_POINTER(C_LOC(SCR(KLFREE)),iSCRf,[1])
      CALL CNFSTR(ICNR,ITPR,iSCRra,iSCRrb,
     *            NORB,NAEL,NBEL,NDETR,IPRODT,iSCRf,
     *            SCR(KLISR),IPREXH)
      NULLIFY(iSCRra,iSCRrb,iSCRf)
*
** Hamiltonin matrix over determinants of ths configuration
*
      KLDHD = KLFREE
      KLFREE = KLFREE + NDETL*NDETR
      CALL COMBINATIONS(ICOMBI,PSIGN)
      CALL C_F_POINTER(C_LOC(SCR(KLDTLA)),iSCRla,[1])
      CALL C_F_POINTER(C_LOC(SCR(KLDTLB)),iSCRlb,[1])
      CALL C_F_POINTER(C_LOC(SCR(KLDTRA)),iSCRra,[1])
      CALL C_F_POINTER(C_LOC(SCR(KLDTRB)),iSCRrb,[1])
      CALL C_F_POINTER(C_LOC(SCR(KLFREE)),iSCRf,[1])
      CALL DIHDJ_MOLCAS(iSCRla,iSCRlb,NDETL,
     *           iSCRra,iSCRrb,NDETR,
     *           NAEL,NBEL,iSCRf,NORB,ONEBOD,
     *           SCR(KLDHD),0,0,ECORE,ICOMBI,PSIGN,0,TUVX,ExFac,IREOTS)
      NULLIFY(iSCRla,iSCRlb,iSCRra,iSCRrb,iSCRf)
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
      END SUBROUTINE CNHCN_INTERNAL
*
      END
