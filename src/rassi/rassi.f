************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE RASSI(IRETURN)

      !> module dependencies
#ifdef _DMRG_
      use qcmaquis_interface_cfg
      use qcmaquis_interface_environment, only:
     &    read_dmrg_info
#endif

      IMPLICIT REAL*8 (A-H,O-Z)
C Matrix elements over RAS wave functions.
C RAS state interaction.
#include "rasdim.fh"
#include "cntrl.fh"
#include "Files.fh"
#include "Morsel.fh"
#include "Struct.fh"
#include "SysDef.fh"
#include "WrkSpc.fh"
#include "rassi.fh"
#include "prgm.fh"
#include "rasdef.fh"
#include "jobin.fh"
#include "symmul.fh"
#include "rassiwfn.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='RASSI')
      Logical Fake_CMO2
      COMMON / CHO_JOBS / Fake_CMO2
#ifdef _DMRG_
      logical dmrg_interface_exists
#endif

      IRETURN=20

      Call StatusLine('RASSI:','Starting calculation')

      CALL GETPRINTLEVEL

      CALL QENTER(ROUTINE)

C Greetings. Default settings. Initialize data sets.
      CALL INIT_RASSI()

#ifdef _DMRG_
      !> initialize DMRG
      inquire(file="dmrg_interface.parameters",
     &        exist=dmrg_interface_exists)
      if(dmrg_interface_exists) call read_dmrg_info()
#endif

C Read and check keywords etc. from stdin. Print out.
      CALL INPCTL_RASSI()

#ifdef _DMRG_
!     !> disable Hamiltonian calculation for rassi-dmrg for the time
!     being. TODO: FIXME: requires general 2-particle RDMs - not
!     difficult but extra work...
!     if(dodmrg) ifham = .false.
#endif
CSVC: prepare HDF5 wavefunction file
      CALL CRE_RASSIWFN

C--------  RAS wave function section --------------------------
C First, in a double loop over the states, compute any matrix
C elements required for the input RASSCF state pairs.
C Needed generalized transition density matrices are computed by
C GTDMCTL. They are written on unit LUTDM.
C Needed matrix elements are computed by PROPER.
      NPROPSZ=NSTATE*NSTATE*NPROP
      CALL GETMEM('Prop','Allo','Real',LPROP,NPROPSZ)
      CALL DCOPY_(NPROPSZ,0.0D0,0,WORK(LPROP),1)
C Loop over jobiphs JOB1:
      DO JOB1=1,NJOB
        DO JOB2=1,JOB1

        Fake_CMO2 = JOB1.eq.JOB2  ! MOs1 = MOs2  ==> Fake_CMO2=.true.

C Compute generalized transition density matrices, as needed:
          CALL GTDMCTL(WORK(LPROP),JOB1,JOB2)
        END DO
      END DO

#ifdef _HDF5_
      CALL mh5_put_dset_array_real(wfn_overlap,
     &     OVLP(1:NSTATE,1:NSTATE),[NSTATE,NSTATE],[0,0])
#endif
      Call Put_dArray('State Overlaps',OVLP(1:NSTATE,1:NSTATE),
     &                NSTATE*NSTATE)

      IF(TRACK) CALL TRACK_STATE()
      IF(TRACK.OR.ONLY_OVERLAPS) THEN

C       Print the overlap matrix here, since MECTL is skipped
        IF(IPGLOB.GE.USUAL) THEN
          WRITE(6,*)
          WRITE(6,*)'     OVERLAP MATRIX FOR THE ORIGINAL STATES:'
          WRITE(6,*)
          DO ISTATE=1,NSTATE
            WRITE(6,'(5(1X,F15.8))')(OVLP(ISTATE,J),J=1,ISTATE)
          END DO
        END IF
        GOTO 100
      END IF

C Property matrix elements:
      Call StatusLine('RASSI:','Computing matrix elements.')
      CALL MECTL(WORK(LPROP))

C--------  SI wave function section --------------------------
C In a second section, if Hamiltonian elements were requested,
C then also a set of secular equations are solved. This gives
C a set of non-interacting, orthonormal wave functions expressed
C as linear combinations of the input wave functions. These
C results are written out, as well as the matrix elements  of
C the eigenstates, and if requested, their density matrices
C and perhaps GTDM's.

C Hamiltonian matrix elements, eigenvectors:
      IF(IFHAM) THEN
        Call StatusLine('RASSI:','Computing Hamiltonian.')
        CALL EIGCTL(WORK(LPROP))
      END IF

C Natural orbitals, if requested:
      IF(NATO) THEN
C CALCULATE AND WRITE OUT NATURAL ORBITALS.
        CALL GETMEM('DMAT  ','ALLO','REAL',LDMAT,NBSQ)
        CALL GETMEM('TDMZZ ','ALLO','REAL',LTDMZZ,NTDMZZ)
        CALL GETMEM('VNAT  ','ALLO','REAL',LVNAT,NBSQ)
        CALL GETMEM('OCC   ','ALLO','REAL',LOCC,NBST)
        CALL NATORB_RASSI(WORK(LDMAT),WORK(LTDMZZ),WORK(LVNAT),
     &                    WORK(LOCC))
        CALL NATSPIN_RASSI(WORK(LDMAT),WORK(LTDMZZ),WORK(LVNAT),
     &                    WORK(LOCC))
        CALL GETMEM('DMAT  ','FREE','REAL',LDMAT,NBSQ)
        CALL GETMEM('TDMZZ ','FREE','REAL',LTDMZZ,NTDMZZ)
        CALL GETMEM('VNAT  ','FREE','REAL',LVNAT,NBSQ)
        CALL GETMEM('OCC   ','FREE','REAL',LOCC,NBST)
      END IF
C Bi-natural orbitals, if requested:
      IF(BINA) THEN
        CALL BINAT()
      END IF
      IF(.NOT.IFHAM) GOTO 100

C-------- Spin-Orbit calculations   --------------------------
C In this third section, if spin-orbit coupling parameters were
C computed by the previous sections, then additionally the
C spin-orbit eigenfunctions and levels are computed.

C Nr of spin states and division of loops:
      NSS=0
      LOOPDIVIDE_TEMP = 0
      DO ISTATE=1,NSTATE
       JOB=JBNUM(ISTATE)
       MPLET=MLTPLT(JOB)
       NSS=NSS+MPLET
       IF(ISTATE.GT.LOOPDIVIDE) CYCLE
       LOOPDIVIDE_TEMP = LOOPDIVIDE_TEMP + MPLET
      END DO
      LOOPDIVIDE = LOOPDIVIDE_TEMP

      CALL GETMEM('UTOTR','ALLO','REAL',LUTOTR,NSS**2)
      CALL GETMEM('UTOTI','ALLO','REAL',LUTOTI,NSS**2)
      CALL GETMEM('SOENE','ALLO','REAL',LSOENE,NSS)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LUTOTR),1)
      CALL DCOPY_(NSS   ,1.0D0,0,WORK(LUTOTR),NSS+1)
      CALL DCOPY_(NSS**2,0.0D0,0,WORK(LUTOTI),1)
      CALL DCOPY_(NSS   ,0.0D0,0,WORK(LSOENE),1)

      IF(IFSO) THEN
        Call StatusLine('RASSI:','Computing SO Hamiltonian.')
        CALL SOEIG(WORK(LPROP),WORK(LUTOTR),WORK(LUTOTI),
     &             WORK(LSOENE),NSS)
      END IF

      CALL PRPROP(WORK(LPROP),WORK(LUTOTR),WORK(LUTOTI),
     &            WORK(LSOENE),NSS)


C Plot SO-Natural Orbitals if requested
C Will also handle mixing of states (sodiag.f)
      IF(SONATNSTATE.GT.0) THEN
        CALL DO_SONATORB(NSS, LUTOTR, LUTOTI)
      END IF


      CALL GETMEM('UTOTR','FREE','REAL',LUTOTR,NSS**2)
      CALL GETMEM('UTOTI','FREE','REAL',LUTOTI,NSS**2)
      CALL GETMEM('SOENE','FREE','REAL',LSOENE,NSS)

CIgorS 02/10-2007  Begin----------------------------------------------C
C   Trajectory Surface Hopping                                        C
C                                                                     C
C   Turns on the procedure if the Keyword HOP was specified.          C
C                                                                     C
      IF (HOP) then
        Call StatusLine('RASSI:','Trajectory Surface Hopping')
        CALL TSHinit()
      END IF
C                                                                     C
CIgorS End------------------------------------------------------------C
*                                                                      *
************************************************************************
*
* CEH April 2015 DQV diabatization scheme
* This passes the PROP matrix into the DQV diabatization subroutine
* The subroutine will compute a transformation matrix, which is used
* to compute diabats.
*
* The user has to compute x, y, z, xx, yy, zz, 1/r in this order.

      IF (DQVD) then
        Call StatusLine('RASSI:', 'DQV Diabatization')
        CALL DQVDiabat(WORK(LPROP))
      END IF
*                                                                      *
************************************************************************
*

 100  CONTINUE
      CALL GETMEM('Prop','Free','Real',LPROP,NPROPSZ)
      CALL GETMEM('NilPt','FREE','REAL',LNILPT,1)
      CALL GETMEM('INilPt','FREE','INTE',LINILPT,1)

#ifdef _DMRG_
      !> finalize DMRG
      if(dmrg_interface_exists.and.
     &  allocated(dmrg_external%dmrg_state_specific))
     &  deallocate(dmrg_external%dmrg_state_specific)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Close dafiles.
*
      Call DaClos(LuScr)
c jochen 02/15: sonatorb needs LUTDM
c     we'll make it conditional upon the keyword
      IF(SONATNSTATE.GT.0) THEN
        Call DaClos(LuTDM)
      end if
c ... jochen end
      Call DaClos(LuExc)
*                                                                      *
************************************************************************
*                                                                      *
*     PRINT I/O STATISTICS:
      IF ( IPGLOB.GE.USUAL ) CALL FASTIO('STATUS')

      Call StatusLine('RASSI:','Finished.')
      IRETURN=0
      CALL QEXIT(ROUTINE)
      RETURN
      END
