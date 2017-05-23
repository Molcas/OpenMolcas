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
* Copyright (C) 1989,1993, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE CNHCN2(ICNL,ITPL,ICNR,ITPR,CNHCNM,SCR,NEL,
     &                 NAEL,NBEL,INTSPC,
     &                 NINOC,ECORE,IPRODT,DTOC,
     &                 NORB,ICOMBI,PSSIGN,NTERMS,NDIF0,NDIF1,NDIF2,
     &                 NTEST)
*
* Obtain Hamilton matrix over CSFs of configurations ICNL,ICNR
*
* Jeppe Olsen , Summer of 89
*
*. Modified for LUCIA, september 1993
*
      IMPLICIT real*8(A-H,O-Z)
*. Specific input
      DIMENSION ICNL(*),ICNR(*)
*. General input
      DIMENSION IPRODT(*),DTOC(*)
*. Scratch
      DIMENSION SCR(*)
*. Output
      DIMENSION CNHCNM(*)
*. Interface to LUCIA common blocks in order to access strings
#include "detdim.fh"
#include "WrkSpc.fh"
#include "strbas_mclr.fh"
*
#include "spinfo_mclr.fh"
#include "cicisp_mclr.fh"


*
* Length of SCR : 6 * NDET + NDET**2 + NDET*NCSF +
*                 MAX((NDET*NEL + 2*NEL),4*NORB + 2*NEL)
* where the last line arises from the local memory in CNFSTR and DIHDJ,
* respectively
*
      IAGRP = IASTFI(INTSPC)
      IBGRP = IBSTFI(INTSPC)
*
** 1 : Obtain determinants corresponding to configurations
*
      IOPL = ITPL - 1 + MINOP
      IOPR = ITPR - 1 + MINOP
*
      ICLL = (NEL-IOPL)/2
      ICLR = (NEL-IOPR)/2
*
      NDETL = NDPCNT(ITPL)
      NDETR = NDPCNT(ITPR)
*
      NCSFL = NCPCNT(ITPL)
      NCSFR = NCPCNT(ITPR)
*
*
      KLFREE = 1
* :  6* Ndet for holding string numbers and signs
      KLDTLA = KLFREE
      KLFREE = KLFREE + NDETL
*
      KLDTLB = KLFREE
      KLFREE = KLFREE + NDETL
*
      KLISL  = KLFREE
      KLFREE = KLFREE + NDETL
*
      KLDTRA = KLFREE
      KLFREE = KLFREE + NDETR
*
      KLDTRB = KLFREE
      KLFREE = KLFREE + NDETR
*
      KLISR  = KLFREE
      KLFREE = KLFREE + NDETR
*. : for transforming to CSFs
      KLDHD = KLFREE
      KLFREE = KLFREE + NDETL*NDETR
*
      KLCHD = KLFREE
      KLFREE = KLFREE + NCSFL*NDETR
*. : Scratch space for routines called ( DIHDJ, CNFSTR )
      KLROU = KLFREE
      LDIHDJ = 4 * NORB + 2*NEL
      LCNFST = MAX(NDETL,NDETR)*NEL + 2 * NEL
      KLFREE = KLROU + MAX(LDIHDJ,LCNFST)
*
*. Prescreen for CNFs that do not interact
*
      CALL CMP2CN(ICNL,ICLL,IOPL,ICNR,ICLR,IOPR,SCR,NORB,NDIFF,
     &            NTEST)
*
      IF(NDIFF .LE. 2 ) THEN
*. Strings for determinants of these configurations
        CALL CNFSTR_MCLR(ICNL,ITPL,SCR(KLDTLA),SCR(KLDTLB),
     &              NORB,NAEL,NBEL,NDETL,IPRODT,IAGRP,IBGRP,
     &              SCR(KLROU ),SCR(KLISL),NTEST)
        CALL CNFSTR_MCLR(ICNR,ITPR,SCR(KLDTRA),SCR(KLDTRB),
     &              NORB,NAEL,NBEL,NDETR,IPRODT,IAGRP,IBGRP,
     &              SCR(KLROU ),SCR(KLISR),NTEST)
*
** Hamiltonian matrix over determinants of the configurations
*
        NTESTP = MAX(0,NTEST-5)
        isym=0       ! eaw
        ecorep=0.0d0 ! eaw
        CALL DIHDJ2_MCLR(SCR(KLDTLA),SCR(KLDTLB),NDETL,
     &              SCR(KLDTRA),SCR(KLDTRB),NDETR,
     &              NAEL,NBEL,SCR(KLROU+2*NEL),LWORK,NORB,
     &              SCR(KLDHD),ISYM,0,ECOREP,ICOMBI,PSSIGN,
     &              iWORK(KOCSTR(IAGRP)),iWORK(KOCSTR(IBGRP)),
     &              iWORK(KOCSTR(IAGRP)),iWORK(KOCSTR(IBGRP)),
     &              0,IDUMMY,IDUMMY,IDUMMY,IDUMMY,SCR(KLROU),
     &              SCR(KLROU+NEL),NDIF0,NDIF1,NDIF2,NTEST)
*
** Transform matrix to CSF basis
*
** : sign changes
        CALL DGMM2 (SCR(KLDHD),SCR(KLDHD),SCR(KLISL),1,NDETL,NDETR)
        CALL DGMM2 (SCR(KLDHD),SCR(KLDHD),SCR(KLISR),2,NDETL,NDETR)
        IPL = 1
        DO 200 JTYP = 1, ITPL - 1
          JCSF = NCPCNT(JTYP)
          JDET = NDPCNT(JTYP)
          IPL = IPL + JCSF*JDET
  200   CONTINUE
        IF(ITPR.EQ. ITPL) THEN
          IPR = IPL
        ELSE
          IPR = 1
          DO 300 JTYP = 1, ITPR - 1
          JCSF = NCPCNT(JTYP)
          JDET = NDPCNT(JTYP)
          IPR = IPR + JCSF*JDET
  300     CONTINUE
        END IF
*
        CALL MATML4(SCR(KLCHD),DTOC(IPL),SCR(KLDHD),
     &              NCSFL,NDETR,NDETL,NCSFL,NDETL,NDETR,1)
        CALL MATML4(CNHCNM,SCR(KLCHD),DTOC(IPR),
     &              NCSFL,NCSFR,NCSFL,NDETR,NDETR,NCSFR,0)
*
      ELSE IF(NDIFF.GT. 2) THEN
        CALL SETVEC(CNHCNM,0.0D0,NCSFL*NCSFR)
      END IF
*
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer(NINOC)
        CALL Unused_real(ECORE)
        CALL Unused_integer(NTERMS)
      END IF
      END
