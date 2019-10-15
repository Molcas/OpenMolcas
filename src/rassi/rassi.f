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
      use rassi_aux
      use kVectors
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
#include "stdalloc.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='RASSI')
      Logical Fake_CMO2,CLOSEONE
      COMMON / CHO_JOBS / Fake_CMO2
      INTEGER IRC,ISFREEUNIT
      EXTERNAL ISFREEUNIT
      Real*8, Allocatable:: PROP(:,:,:), EigVec(:,:), USOR(:,:),
     &                      USOI(:,:), OVLP(:,:), DYSAMPS(:,:),
     &                      ENERGY(:)
      Integer, Allocatable:: IDDET1(:)
*                                                                      *
************************************************************************
*                                                                      *
*     Prolouge
*
      IRETURN=20

      Call StatusLine('RASSI:','Starting calculation')

      CALL GETPRINTLEVEL

      CALL QENTER(ROUTINE)

C Greetings. Default settings. Initialize data sets.
      CALL INIT_RASSI()

      CLOSEONE=.FALSE.
      IRC=-1
      CALL OPNONE(IRC,0,'ONEINT',LUONE)
      IF (IRC.NE.0) Then
         WRITE (6,*) 'RASSI: Error opening file'
         CALL ABEND()
      END IF
      CLOSEONE=.TRUE.

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
      Call mma_allocate(OVLP,NSTATE,NSTATE,Label='OVLP')
      Call mma_allocate(DYSAMPS,NSTATE,NSTATE,Label='DYSAMPS')
      Call mma_allocate(EigVec,nState,nState,Label='EigVec')
      LEIGVEC=ip_of_work(EigVec(1,1))
      Call mma_allocate(ENERGY,nState,Label='Energy')
      Call mma_allocate(TocM,NSTATE*(NSTATE+1)/2,Label='TocM')

      NPROPSZ=NSTATE*NSTATE*NPROP
      Call mma_allocate(PROP,NSTATE,NSTATE,NPROP,LABEL='Prop')
      Prop(:,:,:)=0.0D0
*                                                                      *
************************************************************************
*                                                                      *
C Number of basis functions
      NZ=0                      ! (NBAS is already used...)
      DO ISY=1,NSYM
         NZ=NZ+NBASF(ISY)
      END DO
*
      IF (DYSO) THEN
       LSFDYS=0
       CALL GETMEM('SFDYS','ALLO','REAL',LSFDYS,NZ*NSTATE*NSTATE)
      END IF
*
C Loop over jobiphs JOB1:
      Call mma_allocate(IDDET1,nState,Label='IDDET1')
      IDISK=0  ! Initialize disk address for TDMs.
      DO JOB1=1,NJOB
        DO JOB2=1,JOB1

        Fake_CMO2 = JOB1.eq.JOB2  ! MOs1 = MOs2  ==> Fake_CMO2=.true.

C Compute generalized transition density matrices, as needed:
          CALL GTDMCTL(PROP,JOB1,JOB2,OVLP,DYSAMPS,
     &                 WORK(LSFDYS),NZ,Work(LHAM),iWork(lIDDET1),
     &                 IDISK)
        END DO
      END DO
      Call mma_deallocate(IDDET1)

#ifdef _HDF5_
      CALL mh5_put_dset_array_real(wfn_overlap,
     &     OVLP,[NSTATE,NSTATE],[0,0])
#endif
      Call Put_dArray('State Overlaps',OVLP,
     &                NSTATE*NSTATE)

      IF(TRACK) CALL TRACK_STATE(LOVLP)
      IF(TRACK.OR.ONLY_OVERLAPS) THEN

C       Print the overlap matrix here, since MECTL is skipped
        IF(IPGLOB.GE.USUAL) THEN
          WRITE(6,*)
          WRITE(6,*)'     OVERLAP MATRIX FOR THE ORIGINAL STATES:'
          WRITE(6,*)
          DO ISTATE=1,NSTATE
            WRITE(6,'(5(1X,F15.8))')(Ovlp(j,iState),j=1,istate)
          END DO
        END IF
        GOTO 100
      END IF

C Property matrix elements:
      Call StatusLine('RASSI:','Computing matrix elements.')
      CALL MECTL(PROP,OVLP,WORK(LHAM),WORK(LESHFT))
*                                                                      *
************************************************************************
*                                                                      *
*     Spin-free section
*
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
        CALL EIGCTL(PROP,OVLP,DYSAMPS,Work(LHAM),
     &              EIGVEC,ENERGY)
      END IF


! +++ J. Creutzberg, J. Norell - 2018
! Write the spin-free Dyson orbitals to .DysOrb and .molden
! files if requested
*----------------------------------------------------------------

      IF (DYSEXPORT) THEN

       CALL WRITEDYS(DYSAMPS,WORK(LSFDYS),NZ,ENERGY)

      END IF
! +++


*---------------------------------------------------------------------*
C Natural orbitals, if requested:
      IF(NATO) THEN
C CALCULATE AND WRITE OUT NATURAL ORBITALS.
        CALL GETMEM('DMAT  ','ALLO','REAL',LDMAT,NBSQ)
        CALL GETMEM('TDMZZ ','ALLO','REAL',LTDMZZ,NTDMZZ)
        CALL GETMEM('VNAT  ','ALLO','REAL',LVNAT,NBSQ)
        CALL GETMEM('OCC   ','ALLO','REAL',LOCC,NBST)
        CALL NATORB_RASSI(WORK(LDMAT),WORK(LTDMZZ),WORK(LVNAT),
     &                    WORK(LOCC),EIGVEC)
        CALL NATSPIN_RASSI(WORK(LDMAT),WORK(LTDMZZ),WORK(LVNAT),
     &                    WORK(LOCC),EIGVEC)
        CALL GETMEM('DMAT  ','FREE','REAL',LDMAT,NBSQ)
        CALL GETMEM('TDMZZ ','FREE','REAL',LTDMZZ,NTDMZZ)
        CALL GETMEM('VNAT  ','FREE','REAL',LVNAT,NBSQ)
        CALL GETMEM('OCC   ','FREE','REAL',LOCC,NBST)
      END IF
C Bi-natural orbitals, if requested:
      IF (BINA) CALL BINAT()
*
      IF(.NOT.IFHAM) GOTO 100

*                                                                      *
************************************************************************
*                                                                      *
*      Spin_Orbit section
*
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

      Call mma_allocate(USOR,NSS,NSS,Label='USOR')
      LUTOTR=ip_of_Work(USOR(1,1))
      USOR(:,:)=0.0D0
      For All (i=1:NSS) USOR(i,i)=1.0D0
      Call mma_allocate(USOI,NSS,NSS,Label='USOI')
      LUTOTI=ip_of_Work(USOI(1,1))
      USOI(:,:)=0.0D0
      CALL GETMEM('SOENE','ALLO','REAL',LSOENE,NSS)
      CALL DCOPY_(NSS   ,[0.0D0],0,WORK(LSOENE),1)

      IF(IFSO) THEN
        Call StatusLine('RASSI:','Computing SO Hamiltonian.')
        CALL SOEIG(PROP,USOR,USOI,WORK(LSOENE),NSS,ENERGY)
      END IF

! +++ J. Norell - 2018
C Make the SO Dyson orbitals and amplitudes from the SF ones

      IF (DYSO.AND.IFSO) THEN
       LSODYSAMPS=0
       CALL GETMEM('SODYSAMPS','ALLO','REAL',LSODYSAMPS,NSS*NSS)
       LSODYSAMPSR=0
       CALL GETMEM('SODYSAMPSR','ALLO','REAL',LSODYSAMPSR,NSS*NSS)
       LSODYSAMPSI=0
       CALL GETMEM('SODYSAMPSI','ALLO','REAL',LSODYSAMPSI,NSS*NSS)

       CALL SODYSORB(NSS,USOR,USOI,DYSAMPS,
     &     WORK(LSFDYS),NZ,WORK(LSODYSAMPS),
     &     WORK(LSODYSAMPSR),WORK(LSODYSAMPSI),WORK(LSOENE))
      END IF

      IF (DYSO) THEN
       CALL GETMEM('SFDYS','FREE','REAL',LSFDYS,NZ*NSTATE*NSTATE)
      END IF
! +++

      CALL PRPROP(PROP,USOR,USOI,
     &            WORK(LSOENE),NSS,OVLP,WORK(LSODYSAMPS),
     &            ENERGY,iWork(lJBNUM),EigVec)

C Plot SO-Natural Orbitals if requested
C Will also handle mixing of states (sodiag.f)
      IF(SONATNSTATE.GT.0) THEN
        CALL DO_SONATORB(NSS,USOR,USOI)
      END IF

      Call mma_deallocate(USOR)
      Call mma_deallocate(USOI)
      CALL GETMEM('SOENE','FREE','REAL',LSOENE,NSS)
      IF (DYSO.AND.IFSO) THEN
       CALL GETMEM('SODYSAMPS','FREE','REAL',LSODYSAMPS,NSS*NSS)
       CALL GETMEM('SODYSAMPSR','FREE','REAL',LSODYSAMPSR,NSS*NSS)
       CALL GETMEM('SODYSAMPSI','FREE','REAL',LSODYSAMPSI,NSS*NSS)
      END IF
*                                                                      *
************************************************************************
*                                                                      *
CIgorS 02/10-2007  Begin----------------------------------------------C
C   Trajectory Surface Hopping                                        C
C                                                                     C
C   Turns on the procedure if the Keyword HOP was specified.          C
C                                                                     C
      IF (HOP) then
        Call StatusLine('RASSI:','Trajectory Surface Hopping')
        CALL TSHinit(ENERGY)
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
        CALL DQVDiabat(PROP,WORK(LHAM))
      END IF
*                                                                      *
************************************************************************
*                                                                      *
*     EPILOUGE                                                         *
*                                                                      *
************************************************************************
*                                                                      *
 100  CONTINUE
      Call mma_deallocate(Ovlp)
      Call mma_deallocate(DYSAMPS)
      Call GetMem('HAM','Free','Real',LHAM,NSTATE2)
      Call mma_deallocate(EigVec)
      Call mma_deallocate(Energy)
      Call GetMem('ESHFT','Free','Real',LESHFT,NSTATE)
      Call GetMem('HDIAG','Free','Real',LHDIAG,NSTATE)
      Call mma_deallocate(jDisk_TDM)
      Call GetMem('JBNUM','Free','Inte',LJBNUM,NSTATE)
      Call GetMem('LROOT','Free','Inte',LLROOT,NSTATE)
      Call mma_deallocate(TocM)
      Call mma_deallocate(Prop)
      CALL GETMEM('NilPt','FREE','REAL',LNILPT,1)
      CALL GETMEM('INilPt','FREE','INTE',LINILPT,1)

      If (Do_SK) Call mma_deallocate(k_Vector)

      IF (CLOSEONE) THEN
         IRC=-1
         CALL CLSONE(IRC,0)
         IF (IRC.NE.0) Then
            WRITE (6,*) 'RASSI: Error opening file'
            CALL ABEND()
         END IF
      END IF

#ifdef _DMRG_
!     !> finalize MPS-SI interface
      if (doDMRG)then
        call finalize_dmrg()
        call qcmaquis_info_deinit
      end if
#endif
!     > free memory (if allocated at all - currently only for QD-NEVPT2 as ref wfn)
      call deinit_mspt2_eigenvectors()
*                                                                      *
************************************************************************
*                                                                      *
*     Close dafiles.
*
      Call DaClos(LuScr)
      IF ((SONATNSTATE.GT.0).OR.NATO.OR.Do_TMOM) Then
         Call DaClos(LuTDM)
         If (Allocated(JOB_INDEX)) Call mma_deallocate(JOB_INDEX)
         If (Allocated(CMO1)) Call mma_deallocate(CMO1)
         If (Allocated(CMO2)) Call mma_deallocate(CMO2)
         If (Allocated(DMAB)) Call mma_deallocate(DMAB)
      End If
      Call DaClos(LuExc)
*                                                                      *
************************************************************************
*                                                                      *
*     PRINT I/O STATISTICS:
      i=iPrintLevel(3)
      CALL FASTIO('STATUS')

      Call StatusLine('RASSI:','Finished.')
      IRETURN=0
      CALL QEXIT(ROUTINE)
      RETURN
      END
