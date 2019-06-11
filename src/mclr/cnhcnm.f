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
* Copyright (C) 1990, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE CNHCNM(HSUB,ISYM,ILCNF,NLCNF,IRCNF,NRCNF,NLCSF,NRCSF,
     &                  SCR,ICONF,NEL,
     &                  IREFSM,NAEL,NBEL,NINOB,NACOB,ECORE,
     &                  IPRODT,DTOC,INTSPC,ICOMBI,PSSIGN,
     &                  NTEST)
*
* Calculate  Hamiltonian block defined by configuration
* lists ILCNF,IRCNF
* If ISYM.Ne. 0 only the lower half of the matrix is constructed
*
* Jeppe Olsen April 1990
* ========================
*
* No modifications
* ================
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Specific input
      DIMENSION ILCNF(*),IRCNF(*)
*. General input
      DIMENSION ICONF(*),IPRODT(*),DTOC(*)
*.Output
      DIMENSION HSUB(*)
*.Scratch
      DIMENSION SCR(*)
*. Length of scratch : 2 * NEL + MXCSFC                   (used in CNHCNM)
*                    + 6*MXDTFC+MXDTFC**2+MXDTFC+MXCSFC   (used in CNHCN2)
*                    + MAX(MXDTFC*NEL+2*NEL,4*NORB+2*NEL) (used in DIHDJ,CNFSTR)
*.Standard common block
#include "detdim.fh"
#include "spinfo_mclr.fh"
*
      CALL CNHCNM_INTERNAL(SCR)
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(NRCSF)
*
*     This is to allow type punning without an explicit interface
      CONTAINS
      SUBROUTINE CNHCNM_INTERNAL(SCR)
      USE ISO_C_BINDING
      REAL*8, TARGET :: SCR(*)
      INTEGER, POINTER :: iSCRl(:),iSCRr(:)
*
      NTERMS = 0
      NDIF0 = 0
      NDIF1 = 0
      NDIF2 = 0
      NMAX = 0
*. Largest configuration block possible
      MXCSFC = 0
      DO 10 ITYP = 1, NTYP
        MXCSFC = MAX(MXCSFC,NCPCNT(ITYP) )
   10 CONTINUE
*
*
      KLFREE = 1
*
      KLCONF = KLFREE
      KLFREE = KLFREE + NEL
*
      KRCONF = KLFREE
      KLFREE = KLFREE + NEL
*
      KLPHPS = KLFREE
      KLFREE = KLFREE + MXCSFC ** 2
*
*. LHR
      IILB = 1
      CALL C_F_POINTER(C_LOC(SCR(KLCONF)),iSCRl,[1])
      CALL C_F_POINTER(C_LOC(SCR(KRCONF)),iSCRr,[1])
      DO 200 ICNL = 1, NLCNF
        CALL GETCNF(iSCRl,ILTYP,ILCNF(ICNL),ICONF,IREFSM,NEL,NTEST)
        NCSFL = NCPCNT(ILTYP)
        IIRB = 1
        IF(ISYM.EQ.0) THEN
          MXR = NRCNF
        ELSE
          MXR = ICNL
        END IF
        DO 190 ICNR = 1, MXR
          CALL GETCNF(iSCRr,IRTYP,IRCNF(ICNR),ICONF,IREFSM,NEL,NTEST)
          NCSFR = NCPCNT(IRTYP)
          CALL CNHCN2(iSCRl,ILTYP,iSCRr,IRTYP,SCR(KLPHPS),
     &               SCR(KLFREE),NEL,NAEL,NBEL,INTSPC,
     &               NINOB,ECORE,
     &               IPRODT,DTOC,NACOB,ICOMBI,PSSIGN,
     &               NTERMS,MDIF0,MDIF1,MDIF2,NTEST)
          NDIF0 = NDIF0 + MDIF0
          NDIF1 = NDIF1 + MDIF1
          NDIF2 = NDIF2 + MDIF2
*

* Copy to HSUB matrix
          IF(ISYM.NE.0) THEN
*. Copy to lower half format
            DO 160 IIL = 1, NCSFL
              IF(IILB.EQ.IIRB) THEN
                IIRMAX = IIL
              ELSE
                IIRMAX = NCSFR
              END IF
              DO 150 IIR = 1, IIRMAX
                IIRACT = IIRB - 1 + IIR
                IILACT = IILB - 1 + IIL
                ILRO = IILACT*(IILACT-1)/2 + IIRACT
                ILRI = (IIR-1)*NCSFL + IIL
                HSUB(ILRO) = SCR(KLPHPS-1+ILRI)
  150         CONTINUE
  160       CONTINUE
          ELSE
*. Pack to full format
            DO 260 IIL = 1, NCSFL
              DO 250 IIR = 1, NCSFR
                IIRACT = IIRB - 1 + IIR
                IILACT = IILB - 1 + IIL
                ILRO = (IIRACT-1)*NLCSF + IILACT
                ILRI = (IIR-1)*NCSFL + IIL
                HSUB(ILRO) = SCR(KLPHPS-1+ILRI)
  250         CONTINUE
  260       CONTINUE
          END IF
          IIRB = IIRB + NCSFR
  190   CONTINUE
      IILB = IILB + NCSFL
 200  CONTINUE
      NULLIFY(iSCRl,iSCRr)
*
*
      RETURN
      END SUBROUTINE CNHCNM_INTERNAL
*
      END
