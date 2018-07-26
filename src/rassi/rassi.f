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
      use qcmaquis_interface_environment, only: finalize_dmrg
      use qcmaquis_info, only : qcmaquis_info_deinit
#endif
      use mspt2_eigenvectors, only : deinit_mspt2_eigenvectors

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

      IRETURN=20

      Call StatusLine('RASSI:','Starting calculation')

      CALL GETPRINTLEVEL

      CALL QENTER(ROUTINE)

C Greetings. Default settings. Initialize data sets.
      CALL INIT_RASSI()

C Read and check keywords etc. from stdin. Print out.
      CALL INPCTL_RASSI()

CSVC: prepare HDF5 wavefunction file
      CALL CRE_RASSIWFN

C--------  RAS wave function section --------------------------
C First, in a double loop over the states, compute any matrix
C elements required for the input RASSCF state pairs.
C Needed generalized transition density matrices are computed by
C GTDMCTL. They are written on unit LUTDM.
C Needed matrix elements are computed by PROPER.
      NSTATE2=NSTATE*NSTATE
      Call GetMem('OVLP','Allo','Real',LOVLP,NSTATE2)
      Call GetMem('DYSAMPS','Allo','Real',LDYSAMPS,NSTATE2)
      Call GetMem('EIGVEC','Allo','Real',LEIGVEC,NSTATE2)
      Call GetMem('ENERGY','Allo','Real',LENERGY,NSTATE)
      Call GetMem('ITOCM','Allo','Inte',liTocM,NSTATE*(NSTATE+1)/2)
      Call GetMem('IDDET1','Allo','Inte',lIDDET1,NSTATE)

      NPROPSZ=NSTATE*NSTATE*NPROP
      CALL GETMEM('Prop','Allo','Real',LPROP,NPROPSZ)
      CALL DCOPY_(NPROPSZ,0.0D0,0,WORK(LPROP),1)
C Loop over jobiphs JOB1:
      DO JOB1=1,NJOB
        DO JOB2=1,JOB1

        Fake_CMO2 = JOB1.eq.JOB2  ! MOs1 = MOs2  ==> Fake_CMO2=.true.

C Compute generalized transition density matrices, as needed:
          CALL GTDMCTL(WORK(LPROP),JOB1,JOB2,WORK(LOVLP),WORK(LDYSAMPS),
     &                Work(LHAM),iWork(lIDDET1))
        END DO
      END DO
      Call GetMem('IDDET1','Free','Inte',lIDDET1,NSTATE)

#ifdef _HDF5_
      CALL mh5_put_dset_array_real(wfn_overlap,
     &     WORK(LOVLP),[NSTATE,NSTATE],[0,0])
#endif
      Call Put_dArray('State Overlaps',Work(LOVLP),
     &                NSTATE*NSTATE)

      IF(TRACK) CALL TRACK_STATE(Work(LOVLP))
      IF(TRACK.OR.ONLY_OVERLAPS) THEN

C       Print the overlap matrix here, since MECTL is skipped
        IF(IPGLOB.GE.USUAL) THEN
          WRITE(6,*)
          WRITE(6,*)'     OVERLAP MATRIX FOR THE ORIGINAL STATES:'
          WRITE(6,*)
          DO ISTATE=0,NSTATE-1
            iadr=LOVLP+istate*nstate
            WRITE(6,'(5(1X,F15.8))')(Work(iadr+j),j=0,istate)
          END DO
        END IF
        GOTO 100
      END IF

C Property matrix elements:
      Call StatusLine('RASSI:','Computing matrix elements.')
      CALL MECTL(WORK(LPROP),WORK(LOVLP),WORK(LHAM),WORK(LESHFT))

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
        CALL EIGCTL(WORK(LPROP),WORK(LOVLP),WORK(LDYSAMPS),Work(LHAM),
     &              Work(LEIGVEC),WORK(LENERGY))
      END IF

C Natural orbitals, if requested:
      IF(NATO) THEN
C CALCULATE AND WRITE OUT NATURAL ORBITALS.
        CALL GETMEM('DMAT  ','ALLO','REAL',LDMAT,NBSQ)
        CALL GETMEM('TDMZZ ','ALLO','REAL',LTDMZZ,NTDMZZ)
        CALL GETMEM('VNAT  ','ALLO','REAL',LVNAT,NBSQ)
        CALL GETMEM('OCC   ','ALLO','REAL',LOCC,NBST)
        CALL NATORB_RASSI(WORK(LDMAT),WORK(LTDMZZ),WORK(LVNAT),
     &                    WORK(LOCC),WORK(LEIGVEC))
        CALL NATSPIN_RASSI(WORK(LDMAT),WORK(LTDMZZ),WORK(LVNAT),
     &                    WORK(LOCC),WORK(LEIGVEC))
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
       JOB=iWork(lJBNUM+ISTATE-1)
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
     &             WORK(LSOENE),NSS,WORK(LENERGY))
      END IF

! +++ J. Norell 19/7 - 2018
C Make the SO Dyson orbitals and amplitudes from the SF ones

      LSODYSAMPS=0
      CALL GETMEM('SODYSAMPS','ALLO','REAL',LSODYSAMPS,NSS*NSS)

      ! Number of basis functions
      NZ=0 ! (NBAS is already used...)
      DO ISY=1,NSYM
       NZ=NZ+NBASF(ISY)
      END DO
      LSFDYS=0
      CALL GETMEM('SFDYS','ALLO','REAL',LSFDYS,NZ*NSTATE*NSTATE)
      CALL DO_SODYSORB(NSS,LUTOTR,LUTOTI,WORK(LDYSAMPS),
     &                 WORK(LSFDYS),NZ,WORK(LSODYSAMPS))

      CALL GETMEM('SFDYS','FREE','REAL',LSFDYS,NZ*NSTATE*NSTATE)
! +++

      CALL PRPROP(WORK(LPROP),WORK(LUTOTR),WORK(LUTOTI),
     &            WORK(LSOENE),NSS,WORK(LOVLP),WORK(LSODYSAMPS),
     &            WORK(LENERGY),iWork(lJBNUM))

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
        CALL TSHinit(WORK(LENERGY))
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
        CALL DQVDiabat(WORK(LPROP),WORK(LHAM))
      END IF
*                                                                      *
************************************************************************
*

 100  CONTINUE
      Call GetMem('OVLP','Free','Real',LOVLP,NSTATE2)
      Call GetMem('DYSAMPS','Free','Real',LDYSAMPS,NSTATE2)
      Call GetMem('HAM','Free','Real',LHAM,NSTATE2)
      Call GetMem('EIGVEC','Free','Real',LEIGVEC,NSTATE2)
      Call GetMem('ENERGY','Free','Real',LENERGY,NSTATE)
      Call GetMem('ESHFT','Free','Real',LESHFT,NSTATE)
      Call GetMem('HDIAG','Free','Real',LHDIAG,NSTATE)
      Call GetMem('IDTDM','Free','Inte',lIDTDM,NSTATE2)
      Call GetMem('IDDYS','Free','Inte',LIDDYS,NSTATE2)
      Call GetMem('JBNUM','Free','Inte',LJBNUM,NSTATE)
      Call GetMem('LROOT','Free','Inte',LLROOT,NSTATE)
      Call GetMem('ITOCM','Free','Inte',liTocM,NSTATE*(NSTATE+1)/2)
      CALL GETMEM('Prop','Free','Real',LPROP,NPROPSZ)
      CALL GETMEM('NilPt','FREE','REAL',LNILPT,1)
      CALL GETMEM('INilPt','FREE','INTE',LINILPT,1)
      CALL GETMEM('SODYSAMPS','FREE','REAL',LSODYSAMP,NSS*NSS)

#ifdef _DMRG_
!     !> finalize MPS-SI interface
      if (doDMRG)then
        call finalize_dmrg()
        call qcmaquis_info_deinit
      end if
#endif
      !> free memory (if allocated at all - currently only for QD-NEVPT2 as ref wfn)
      call deinit_mspt2_eigenvectors()
*                                                                      *
************************************************************************
*                                                                      *
*     Close dafiles.
*
      Call DaClos(LuScr)
      Call DaClos(LuTDM)
      Call DaClos(LUDYS)
c jochen 02/15: sonatorb needs LUTDM
c     we'll make it conditional upon the keyword
      IF((SONATNSTATE.GT.0).OR.NATO) THEN
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
