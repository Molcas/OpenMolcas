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
      use rassi_global_arrays, only: HAM, SFDYS, SODYSAMPS, EIGVEC,
     &                               SODYSAMPSR, SODYSAMPSI,
     &                               PROP, ESHFT, HDIAG, JBNUM, LROOT
      use rassi_aux
      use kVectors
#ifdef _HDF5_
      use Dens2HDF5
#endif
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
      INTEGER IRC
      Real*8, Allocatable:: USOR(:,:),
     &                      USOI(:,:), OVLP(:,:), DYSAMPS(:,:),
     &                      ENERGY(:), DMAT(:), TDMZZ(:),
     &                      VNAT(:),OCC(:), SOENE(:)
      Integer, Allocatable:: IDDET1(:)
*                                                                      *
************************************************************************
*                                                                      *
*     Prolouge
*
      IRETURN=20

      Call StatusLine('RASSI:','Starting calculation')

      CALL GETPRINTLEVEL


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
      Call mma_allocate(OVLP,NSTATE,NSTATE,Label='OVLP')
      Call mma_allocate(DYSAMPS,NSTATE,NSTATE,Label='DYSAMPS')
      Call mma_allocate(EigVec,nState,nState,Label='EigVec')
      Call mma_allocate(ENERGY,nState,Label='Energy')
      Call mma_allocate(TocM,NSTATE*(NSTATE+1)/2,Label='TocM')

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
      IF (DYSO) Call mma_allocate(SFDYS,nZ,nState,nState,Label='SFDYS')

*
C Loop over jobiphs JOB1:
      Call mma_allocate(IDDET1,nState,Label='IDDET1')
      IDISK=0  ! Initialize disk address for TDMs.
      DO JOB1=1,NJOB
        DO JOB2=1,JOB1

        Fake_CMO2 = JOB1.eq.JOB2  ! MOs1 = MOs2  ==> Fake_CMO2=.true.

C Compute generalized transition density matrices, as needed:
          CALL GTDMCTL(PROP,JOB1,JOB2,OVLP,DYSAMPS,NZ,IDDET1,IDISK)
        END DO
      END DO
      Call mma_deallocate(IDDET1)

#ifdef _HDF5_
      CALL mh5_put_dset_array_real(wfn_overlap,
     &     OVLP,[NSTATE,NSTATE],[0,0])
#endif
      Call Put_dArray('State Overlaps',OVLP,NSTATE*NSTATE)

      IF(TRACK) CALL TRACK_STATE(OVLP)
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
      CALL MECTL(PROP,OVLP,HAM,ESHFT)
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
C and perhaps GTDMs.

C Hamiltonian matrix elements, eigenvectors:
      IF(IFHAM) THEN
        Call StatusLine('RASSI:','Computing Hamiltonian.')
        CALL EIGCTL(PROP,OVLP,DYSAMPS,HAM,EIGVEC,ENERGY)
      END IF


! +++ J. Creutzberg, J. Norell - 2018
! Write the spin-free Dyson orbitals to .DysOrb and .molden
! files if requested
*----------------------------------------------------------------

      IF (DYSEXPORT) THEN

       CALL WRITEDYS(DYSAMPS,SFDYS,NZ,ENERGY)

      END IF
! +++


*---------------------------------------------------------------------*
C Natural orbitals, if requested:
      IF(NATO) THEN
C CALCULATE AND WRITE OUT NATURAL ORBITALS.
        Call mma_allocate(DMAT,nBSQ,Label='DMAT')
        Call mma_allocate(TDMZZ,nTDMZZ,Label='TDMZZ')
        Call mma_allocate(VNAT,nBSQ,Label='VNAT')
        Call mma_allocate(OCC,nBST,Label='OCC')
*
        CALL NATORB_RASSI(DMAT,TDMZZ,VNAT,OCC,EIGVEC)
        CALL NATSPIN_RASSI(DMAT,TDMZZ,VNAT,OCC,EIGVEC)
*
        Call mma_deallocate(DMAT)
        Call mma_deallocate(TDMZZ)
        Call mma_deallocate(VNAT)
        Call mma_deallocate(OCC)
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
      DO ISTATE=1,NSTATE
         JOB=JBNUM(ISTATE)
         MPLET=MLTPLT(JOB)
         NSS=NSS+MPLET
      END DO

      Call mma_allocate(USOR,NSS,NSS,Label='USOR')
      USOR(:,:)=0.0D0
      For All (i=1:NSS) USOR(i,i)=1.0D0
      Call mma_allocate(USOI,NSS,NSS,Label='USOI')
      USOI(:,:)=0.0D0
      Call mma_allocate(SOENE,nSS,Label='SOENE')
      SOENE(:)=0.0D0

      IF(IFSO) THEN
        Call StatusLine('RASSI:','Computing SO Hamiltonian.')
        CALL SOEIG(PROP,USOR,USOI,SOENE,NSS,ENERGY)
      END IF

C Store TDMs in HDF5
#ifdef _HDF5_
      Call StoreDens(EigVec)
#endif

! +++ J. Norell - 2018
C Make the SO Dyson orbitals and amplitudes from the SF ones

      IF (DYSO.AND.IFSO) THEN
         Call mma_allocate(SODYSAMPS,NSS,NSS,Label='SODYSAMPS')
         Call mma_allocate(SODYSAMPSR,NSS,NSS,Label='SODYSAMPSR')
         Call mma_allocate(SODYSAMPSI,NSS,NSS,Label='SODYSAMPSI')

         CALL SODYSORB(NSS,USOR,USOI,DYSAMPS,NZ,SOENE)
      END IF

      IF (Allocated(SFDYS)) Call mma_deallocate(SFDYS)
! +++

      CALL PRPROP(PROP,USOR,USOI,SOENE,NSS,OVLP,
     &            ENERGY,JBNUM,EigVec)

C Plot SO-Natural Orbitals if requested
C Will also handle mixing of states (sodiag.f)
      IF(SONATNSTATE.GT.0) THEN
        CALL DO_SONATORB(NSS,USOR,USOI)
      END IF

      Call mma_deallocate(USOR)
      Call mma_deallocate(USOI)
      Call mma_deallocate(SOENE)
      IF (DYSO.AND.IFSO) THEN
         Call mma_deallocate(SODYSAMPS)
         Call mma_deallocate(SODYSAMPSR)
         Call mma_deallocate(SODYSAMPSI)
      END IF
*                                                                      *
************************************************************************
*                                                                      *
*   Trajectory Surface Hopping                                         *
*                                                                      *
*   Turns on the procedure if the Keyword HOP was specified.           *
*                                                                      *
      IF (HOP) then
        Call StatusLine('RASSI:','Trajectory Surface Hopping')
        CALL TSHinit(ENERGY)
      END IF
*                                                                      *
************************************************************************
*                                                                      *
* CEH April 2015 DQV diabatization scheme
* This passes the PROP matrix into the DQV diabatization subroutine
* The subroutine will compute a transformation matrix, which is used
* to compute diabats.
*
* The user has to compute x, y, z, xx, yy, zz, 1/r in this order.

      IF (DQVD) then
        Call StatusLine('RASSI:', 'DQV Diabatization')
        CALL DQVDiabat(PROP,HAM)
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
      Call mma_deallocate(HAM)
      Call mma_deallocate(EigVec)
      Call mma_deallocate(Energy)
      Call mma_deallocate(ESHFT)
      Call mma_deallocate(HDIAG)
      Call mma_deallocate(jDisk_TDM)
      Call mma_deallocate(JBNUM)
      Call mma_deallocate(LROOT)
      Call mma_deallocate(TocM)
      Call mma_deallocate(Prop)

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
      IF (SaveDens) Then
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
      RETURN
      END
