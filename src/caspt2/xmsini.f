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
* Copyright (C) 2012, Per Ake Malmqvist                                *
*               2019, Stefano Battaglia                                *
************************************************************************
      SUBROUTINE XMSINI(HEFF,H0)
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
      REAL*8 HEFF(NSTATE,NSTATE)
      REAL*8 H0(NSTATE,NSTATE)
      LOGICAL IF_TRNSF

      CALL QENTER('XMSINI')

* Load CASSCF MO coefficients
      CALL GETMEM('LCMO','ALLO','REAL',LCMO,NCMO)
      IDISK=IAD1M(1)
      CALL DDAFILE(LUONEM,2,WORK(LCMO),NCMO,IDISK)

* Allocate memory for CI array state averaged 1-RDM
      CALL GETMEM('LCI','ALLO','REAL',LCI,NCONF)
      CALL GETMEM('LDAVE','ALLO','REAL',LDAVE,NDREF)
      CALL DCOPY_(NDREF,[0.0D0],0,WORK(LDAVE),1)

* Set the weight for the density averaging
      WGT = 1.0D0/DBLE(NSTATE)

* Loop over all states to compute state-average 1-RDM with equal weights
      DO ISTATE=1,NSTATE

        IF (ISCF.NE.0) THEN
* Special case for a single Slater determinant
          WORK(LCI)=1.0D0
        ELSE
* Get the CI array
          ID=IDCIEX
          DO I=1,ISTATE-1
            CALL DDAFILE(LUCIEX,0,WORK(LCI),NCONF,ID)
          END DO
          CALL DDAFILE(LUCIEX,2,WORK(LCI),NCONF,ID)
        END IF

        IF (IPRGLB.GE.DEBUG) THEN
          WRITE(6,*)
          WRITE(6,*)' CI array of CASSCF state nr. ',MSTATE(ISTATE)
          CALL PRWF_CP2(LSYM,NCONF,WORK(LCI),CITHR)
        END IF

* Compute 1-particle active density matrix GAMMA1
        CALL POLY1(WORK(LCI))
* Restructure GAMMA1 as DREF array
        CALL GETDREF(WORK(LDREF))

        IFTEST = 0
        IF ( IFTEST.NE.0 ) THEN
          WRITE(6,*)' DREF for state nr. ',MSTATE(ISTATE)
          DO I=1,NASHT
            WRITE(6,'(1x,14f10.6)')(WORK(LDREF+(I*(I-1))/2+J-1),J=1,I)
          END DO
          WRITE(6,*)
        END IF

* Average the density
        CALL DAXPY_(NDREF,WGT,WORK(LDREF),1,WORK(LDAVE),1)

      END DO


      IFTEST = 0
      IF ( IFTEST.NE.0 ) THEN
        WRITE(6,*)' State averaged 1-RDM'
        DO I=1,NASHT
          WRITE(6,'(1x,14f10.6)')(WORK(LDAVE+(I*(I-1))/2+J-1),J=1,I)
        END DO
        WRITE(6,*)
      END IF

* Copy the state averaged 1-RDM into LDREF and release memory
* for both the CI array and SA 1-RDM
      CALL DCOPY_(NDREF,WORK(LDAVE),1,WORK(LDREF),1)
      CALL GETMEM('LCI','FREE','REAL',LCI,NCONF)
      CALL GETMEM('LDAVE','FREE','REAL',LDAVE,NDREF)


* Build the state-averaged Fock matrix in MO basis
      If (IfChol) then
* INTCTL2 uses TraCho2 and FMatCho to get matrices in MO basis.
        IF_TRNSF=.FALSE.
        CALL INTCTL2(IF_TRNSF)
      Else
* INTCTL1 uses TRAONE and FOCK_RPT2, to get the matrices in MO basis.
        CALL INTCTL1(WORK(LCMO))
        ! CALL DCOPY_(NCMO,WORK(LCMO),1,WORK(LCMOPT2),1)
      End If


* Initialize model space Fock matrix (i.e. the Fock matrix expressed
* in the basis of the selected CASSCF states)
      CALL GETMEM('FOPXMS','ALLO','REAL',LFOPXMS,NSTATE**2)
      CALL DCOPY_(NSTATE**2,[0.0D0],0,WORK(LFOPXMS),1)

* Loop again over all states in the calculation to build FOPXMS
      DO J2=1,NSTATE
        JSTATE=J2

* Compute elements <BRA|FOP|KET> where FOP is the (state averaged)
* Fock operator (i.e. the FIFA array)
        IKET=JSTATE
* Loop over bra functions:
        DO J1=1,NSTATE
          IBRA=J1
* Compute matrix element and put it into FOPXMS and H0:
          FOPEL = 0.0D0
          CALL FOPAB(WORK(LFIFA),IBRA,IKET,FOPEL)
          WORK(LFOPXMS + (J1-1) + (NSTATE*(J2-1))) = FOPEL
          H0(J1,J2) = FOPEL
        END DO

      END DO
* End of loop over states


      IF (IPRGLB.GE.VERBOSE.AND.(NSTATE.GT.1)) THEN
        WRITE(6,*)
        WRITE(6,*)' FOPXMS:'
        DO ISTA=1,NSTATE,5
          IEND=MIN(ISTA+4,NSTATE)
          WRITE(6,*)
          WRITE(6,'(1x,5I16)')(I,I=ISTA,IEND)
          DO I=ISTA,NSTATE
            II0=(I*(I-1))/2
            WRITE(6,'(1x,I3,2X,5F16.8)')
     &            I,(WORK(LFOPXMS+I-1+NSTATE*(J-1)),J=ISTA,IEND)
          END DO
        END DO

        WRITE(6,*)
        WRITE(6,*)' H0 (should be equal to FOPXMS):'
        DO J1=1,NSTATE
          WRITE(6,'(1x,5F16.8)')(H0(J1,J2),J2=1,NSTATE)
        END DO
      END IF

* Transform the CI arrays to make the FOP matrix diagonal.
      CALL GETMEM('EVEC','ALLO','REAL',LEVEC,NSTATE**2)
      CALL DIAFOP(NSTATE,WORK(LFOPXMS),WORK(LEVEC))

* Transform FOPXMS and HEFF in the basis that diagonalizes FOPXMS
      CALL TRANSMAT(H0,WORK(LEVEC),NSTATE)
      CALL TRANSMAT(HEFF,WORK(LEVEC),NSTATE)

      IF (IPRGLB.GE.VERBOSE) THEN
        WRITE(6,*)
        WRITE(6,'(6X,A)')' H0 eigenvectors:'
        DO ISTA=1,NSTATE,5
          IEND=MIN(ISTA+4,NSTATE)
          DO J1=0,NSTATE-1
            WRITE(6,'(6x,5F16.8)')(WORK(LEVEC+J1+NSTATE*(I-1)),
     &                             I=ISTA,IEND)
          END DO
          WRITE(6,*)
        END DO
      END IF

      IF (IPRGLB.GE.VERBOSE) THEN
        WRITE(6,*)
        WRITE(6,*)' Heff[1] in H0 basis:'
        DO J1=1,NSTATE
          WRITE(6,'(1x,5F16.8)')(HEFF(J1,J2),J2=1,NSTATE)
        END DO
        WRITE(6,*)
        WRITE(6,*)' H0 in H0 basis (should be digonal):'
        DO J1=1,NSTATE
          WRITE(6,'(1x,5F16.8)')(H0(J1,J2),J2=1,NSTATE)
        END DO
      END IF

* And then, transform the CI arrays. Assume we can put all the
* original ones in memory, but put the resulting vectors one by
* one in a buffer.
      WRITE(6,*)
      WRITE(6,*)' Mixing the CASSCF states according to XMS...'
      WRITE(6,*)
      CALL GETMEM('CIREF','ALLO','REAL',LCIREF,NSTATE*NCONF)
      DO J=1,NSTATE
        IK=JSTATE_OFF+J
        ID=IDCIEX
        DO I=1,IK-1
          CALL DDAFILE(LUCIEX,0,WORK(LCIREF),NCONF,ID)
        END DO
        CALL DDAFILE(LUCIEX,2,WORK(LCIREF+NCONF*(J-1)),NCONF,ID)
      END DO
      CALL GETMEM('CIXMS','ALLO','REAL',LCIXMS,NCONF)
      DO J=1,NSTATE
        CALL DGEMM_('N','N',NCONF,1,NSTATE,
     &              1.0D0,WORK(LCIREF),NCONF,
     &              WORK(LEVEC+NSTATE*(J-1)),NSTATE,
     &              0.0D0,WORK(LCIXMS),NCONF)
        IK=JSTATE_OFF+J
        ID=IDCIEX
        DO I=1,IK-1
          CALL DDAFILE(LUCIEX,0,WORK(LCIXMS),NCONF,ID)
        END DO
        CALL DDAFILE(LUCIEX,1,WORK(LCIXMS),NCONF,ID)

        IF(IPRGLB.GE.VERBOSE) THEN
          WRITE(6,'(1x,a,i3)')
     &    ' The CI coefficients for the MIXED state nr. ',MSTATE(IK)
          CALL PRWF_CP2(LSYM,NCONF,WORK(LCIXMS),CITHR)
        END IF
      END DO

* Release all memory
      CALL GETMEM('CIREF','FREE','REAL',LCIREF,NSTATE*NCONF)
      CALL GETMEM('CIXMS','FREE','REAL',LCIXMS,NCONF)
      CALL GETMEM('EVEC','FREE','REAL',LEVEC,NSTATE**2)
      CALL GETMEM('FOPXMS','FREE','REAL',LFOPXMS,NSTATE**2)
      CALL GETMEM('LCMO','FREE','REAL',LCMO,NCMO)

      CALL QEXIT('XMSINI')
      RETURN
      END

