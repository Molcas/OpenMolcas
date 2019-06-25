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
* Copyright (C) 2019, Stefano Battaglia                                *
************************************************************************
      SUBROUTINE H0WEIGHTS

      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "intgrl.fh"
#include "eqsolv.fh"
#include "warnings.fh"

      CALL QENTER('H0WEIGHTS')

      IF (IPRGLB.GE.VERBOSE) THEN
        WRITE(6,*)' Entered H0WEIGHTS.'
      END IF

* Initialize array of weights with all zeros
      CALL DCOPY_(NSTATE**2,[0.0D0],0,WORK(LFWGT),1)

      DO ISTATE=1,NSTATE

* Compute weights for constructing H0
        IF (IFFDW) THEN
          EBETA = REFENE(ISTATE)
* Compute normalization factor
          DO JSTATE=1,NSTATE
            EALPHA = REFENE(JSTATE)
            FAC = 0.0D0
            DO KSTATE=1,NSTATE
              EGAMMA = REFENE(KSTATE)
              FAC = FAC + EXP(-ZETAF*(EALPHA - EGAMMA)**2)
            END DO
            IJ = (ISTATE-1) + NSTATE*(JSTATE-1)
            WORK(LFWGT+IJ) = EXP(-ZETAF*(EALPHA - EBETA)**2)/FAC
          END DO
* If it is an XMS-CASPT2 calculation, all the weights are 1/NSTATE
        ELSE IF (IFXMS) THEN
          CALL DCOPY_(NSTATE**2,1.0D0/NSTATE,0,WORK(LFWGT),1)
        ELSE
* It is a normal MS-CASPT2 and the weight vectors are the standard e_1, e_2, ...
          WORK(LFWGT + (NSTATE*(ISTATE-1)) + (ISTATE-1)) = 1.0D0
        END IF

      END DO

        WRITE(6,*)
        WRITE(6,*)' Weights used for H0:'
        DO I=1,NSTATE
          WRITE(6,'(1x,10f8.4)')(WORK(LFWGT + (I-1) + NSTATE*(J-1)),
     &    J=1,NSTATE)
        END DO
        WRITE(6,*)

      CALL QEXIT('H0WEIGHTS')
      RETURN
      END


********************************************************************


      SUBROUTINE H1WEIGHTS

      IMPLICIT REAL*8 (A-H,O-Z)

      LOGICAL IF_TRNSF
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "intgrl.fh"
#include "eqsolv.fh"
#include "warnings.fh"

      CALL QENTER('H1WEIGHTS')

      IF (IPRGLB.GE.VERBOSE) THEN
        WRITE(6,*)' Entered H1WEIGHTS.'
      END IF

* Initialize array of weights with all zeros
      CALL DCOPY_(NSTATE**2,[0.0D0],0,WORK(LVWGT),1)

      if (IFVDW) THEN
* GET ORIGINAL CASSCF CMO COEFFICIENTS.
      CALL GETMEM('LCMO','ALLO','REAL',LCMO,NCMO)
      IDISK=IAD1M(1)
      CALL DDAFILE(LUONEM,2,WORK(LCMO),NCMO,IDISK)
* Also (for temporary back-compatibility with older code) save as
*  'current' CMO data on LUONEM:
      IAD1M(2)=IDISK
      CALL DDAFILE(LUONEM,1,WORK(LCMO),NCMO,IDISK)
      IEOF1M=IDISK
      END IF

***************************************

      DO ISTATE=1,NSTATE

        WRITE(6,*)
        WRITE(6,*)' ISTATE = ',ISTATE
* Compute weights for constructing V
        IF (IFVDW) THEN
* build sa-Fock matrix (we do this at every cycle but in reality is always the same...)
        CALL DCOPY_(NDREF,WORK(LDMIX+(ISTATE-1)*NDREF),1,WORK(LDREF),1)
* This builds FIMO, FAMO, FIFA, ...
          If (IfChol) then
* INTCTL2 uses TraCho2 and FMatCho to get matrices in MO basis.
            IF_TRNSF=.FALSE.
            CALL INTCTL2(IF_TRNSF)
          Else
* INTCTL1 uses TRAONE and FOCK_RPT2, to get the matrices in MO basis.
            CALL INTCTL1(WORK(LCMO))
            CALL DCOPY_(NCMO,WORK(LCMO),1,WORK(LCMOPT2),1)
          End If

          ! EBETA = REFENE(ISTATE)
          DO JSTATE=1,NSTATE
          WRITE(6,*)
          WRITE(6,*)' JSTATE = ',JSTATE
            ! EALPHA = REFENE(JSTATE)
            FIJ = 0.0D0
            EIJ = EI - EJ
            CALL FOPAB(WORK(LFIFA),ISTATE,JSTATE,FIJ)
            WRITE(6,*)
            WRITE(6,*)' FIJ = ',FIJ
* Compute normalization factor FAC (i.e. the denominator)
            FAC = 0.0D0
            DO KSTATE=1,NSTATE
              ! EGAMMA = REFENE(KSTATE)
              FKJ = 0.0D0
              EKJ = EK - EJ
              CALL FOPAB(WORK(LFIFA),KSTATE,JSTATE,FKJ)
              FAC = FAC + EXP(-ZETAV*(EJK/FKJ)**2)
            END DO
            IJ = (ISTATE-1) + NSTATE*(JSTATE-1)
            WORK(LVWGT+IJ) = EXP(-ZETAV*(EIJ/FIJ)**2)/FAC
          END DO
* It is a normal (X)MS-CASPT2 and the weight vectors are the standard e_1, e_2, ...
        ELSE
          WORK(LVWGT + (ISTATE-1) + NSTATE*(ISTATE-1)) = 1.0D0
        END IF

      END DO

      IF (IFVDW) THEN
      CALL GETMEM('LCMO','FREE','REAL',LCMO,NCMO)
      END IF

        WRITE(6,*)
        WRITE(6,*)' Weights used for V:'
        DO I=1,NSTATE
          WRITE(6,'(1x,10f8.4)')(WORK(LVWGT + (I-1) + NSTATE*(J-1)),
     &    J=1,NSTATE)
        END DO
        WRITE(6,*)

      CALL QEXIT('H1WEIGHTS')
      RETURN
      END
