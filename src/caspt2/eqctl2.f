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
* Copyright (C) 2005, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 2005  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE EQCTL2(ICONV)
      IMPLICIT REAL*8 (A-H,O-Z)
C On return, the following data sets will be defined and stored
C on LUSOLV.
C At position IVEC=IRHS, the RHS array, in SR representation.
C At position IVEC=IVECX, the solution array, in SR representation.
C At position IVEC=IVECR, the residual array, in SR representation.
C At position IVEC=IVECC, the solution array, in contravariant rep.
C At position IVEC=IVECC2, the solution array, in covariant repr.
C At position IVEC=IVECW, the RHS array, in contravariant repr.
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "chocaspt2.fh"


      IF (IPRGLB.GE.VERBOSE) THEN
        WRITE(6,'(1X,A)')
        WRITE(6,'(1X,A)') 'Computing the S/B matrices'
        WRITE(6,'(1X,A)') '--------------------------'
      END IF

C Compute S and (possibly) B matrices and write them on LUSBT
C this uses previously stored data on LUSOLV!!
CPAM98      IF(SMATRIX.NE.'NO      ')CALL SBMAT
      CALL GASync
      CALL TIMING(CPU0,CPU,TIO0,TIO)
* PAM14 Necessary to reset NINDEP to conservative estimate. It gets its
* final value after SBDIAG call, but must have its original value when
* calling MKSMAT and MKBMAT.
      DO ICASE=1,13
       DO ISYM=1,NSYM
        NINDEP(ISYM,ICASE)=NASUP(ISYM,ICASE)
        IF(NISUP(ISYM,ICASE).EQ.0) NINDEP(ISYM,ICASE)=0
       END DO
      END DO
      IF(SMATRIX.NE.'NO      ') THEN
        CALL MKSMAT
        CALL MKBMAT
      END IF
C Modify B matrices, if necessary:
      IF(HZERO.EQ.'CUSTOM') THEN
        CALL NEWB
      END IF
      CALL GASync
      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUSBM=CPU1-CPU0
      TIOSBM=TIO1-TIO0

C Linear dependence removal, ON transformation of S matrix,
C and spectral resolution of H0:
      CALL GASync
      CALL TIMING(CPU0,CPU,TIO0,TIO)
      IF(SDECOM.NE.'NO      ') THEN
        CALL SBDIAG
      END IF
      CALL GASync
      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUEIG=CPU1-CPU0
      TIOEIG=TIO1-TIO0
C The transformation matrices have now been computed and
C written at disk addresses given by IDTMAT().
C The B matrices were destroyed here. In their place,
C at the IDBMAT() addresses, we find two sets of diagonal
C energy values. Each consists of: first the energies of
C the active combination, which are eigenvalues of the
C B matrix using overlap S; then the energies of the
C non-active combinations. Usually, only the first set of
C energies is used. The second is used, in a manner specific
C to future modifications, for additional energy parameters
C of more sophisticated H0 or resolvent definitions.
C Modif PAM 961022: If BMATRIX='YES     ', diagonal energies
C are put at IDBMAT(). If BTRANS='YES     ', these are the
C diagonal values of B after transformation to ON basis.
C If furthermore BSPECT='YES     ', the ON basis consists
C of eigenvectors. This is the usual CASPT2 situation.
C However, if only BMATRIX is 'YES     ', then the values
C are the diagonal values of B divided by diagonal values
C of S.

      CALL GASync
      CALL TIMING(CPU0,CPU,TIO0,TIO)
C Non-active part of diagonal elements of H0 are computed
C and written to LUSBT:
      CALL NADIAG
C Modify diagonal elements, if requested:
      IF(HZERO.EQ.'CUSTOM') THEN
        CALL NEWDIA
      END IF
C A second set of energy parameters may now have been
C computed and written to LUSBT.
      CALL GASync
      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUNAD=CPU1-CPU0
      TIONAD=TIO1-TIO0

      IF (IPRGLB.GE.VERBOSE) THEN
        WRITE(6,'(1X,A)')
        WRITE(6,'(1X,A)') 'Computing the right-hand side (RHS) elements'
        WRITE(6,'(1X,A)') '--------------------------------------------'
      END IF

      IRHS  =1
      IVECX =2
      IVECR =3
      IVECC =4
      IVECC2=5
      IVECW =6

      CALL GASync
      CALL TIMING(CPU0,CPU,TIO0,TIO)
C-SVC: initialize the RHS array offsets
      CALL RHS_INIT
C at this point LUSOLV is not used as a temporary disk anymore, so it's
C safe to initialize it (safe to write zeros to LUSOLV or delete it)

C Set up right-hand side matrix elements.
      IF(IfChol .AND. iALGO.eq.1) THEN
        IF (RHSDIRECT) THEN
          IF (NSYM.EQ.1) THEN
            CALL RHSOD_NOSYM(IVECW)
          ELSE
            CALL RHSOD(IVECW)
          END IF
        ELSE
          CALL RHS_ZERO(IVECW)
          CALL RHSALL2(IVECW)
        END IF
      ELSE
        CALL MKRHS(IVECW)
      END IF
      CALL GASync
      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPURHS=CPU1-CPU0
      TIORHS=TIO1-TIO0

      IF (IPRGLB.GE.INSANE) THEN
        WRITE(6,'("DEBUG> ")')
        WRITE(6,'("DEBUG> ",A)')
     &   'Norms of the RHS blocks:'
        CALL RHS_FPRINT('C',IVECW)
      END IF

C-SVC: start PCG routine, set timers.
      CALL GASync
      CALL TIMING(CPU0,CPU,TIO0,TIO)
      CPUSCA=0
      CPULCS=0
      CPUOVL=0
      CPUSGM=0
      CPUVEC=0

C Transform RHS of CASPT2 equations to eigenbasis for H0:
      CALL PTRTOSR(1,IVECW,IRHS)

      IF (IPRGLB.GE.INSANE) THEN
        WRITE(6,'("DEBUG> ")')
        WRITE(6,'("DEBUG> ",A)')
     &   'Norms of the RHS blocks (H0 eigenbasis):'
        CALL RHS_FPRINT('SR',IRHS)
      END IF

      CALL PCG(ICONV)

      IF (ICONV .NE. 0) GOTO 100
      CALL PTRTOC(0,IVECX,IVECC)
      CALL PTRTOC(1,IVECX,IVECC2)

C-SVC: end of PCG routine, compute total time.
      CALL GASync
      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUPCG=CPU1-CPU0
      TIOPCG=TIO1-TIO0

C-SVC: collect and print information on coefficients/denominators
      IF(IPRGLB.GE.USUAL) THEN
        CALL H0SPCT
      END IF

      CALL GASync
      CALL TIMING(CPU0,CPU,TIO0,TIO)
      CALL GASync
      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUSER=CPU1-CPU0
      TIOSER=TIO1-TIO0

  100 CONTINUE
      RETURN
      END
