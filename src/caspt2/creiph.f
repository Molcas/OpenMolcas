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
* Copyright (C) 1997, Per Ake Malmqvist                                *
************************************************************************
      SUBROUTINE CREIPH_CASPT2(Heff,Ueff,U0)
      USE REFWFN
      IMPLICIT REAL*8 (A-H,O-Z)
C Normal operation: A new file, 'JOBMIX', will be created, with the
C CMO's and CI arrays of the JOBIPH, except that the CI arrays have
C been modified. They are now linear combinations of the original ones,
C using coefficients taken from the eigenvectors of the effective
C Hamiltonian.
C Also, replace the original CASSCF energies with CASPT2 or MS-CASPT2
C energies.
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "stdalloc.fh"
      integer JOBIPH, JOBMIX
      real(8) Weight(MxRoot)
      real(8) Heff(Nstate,Nstate),Ueff(Nstate,Nstate),U0(Nstate,Nstate)
      real(8),allocatable :: U0transpose(:,:),Utmp(:,:)

      CALL QENTER('CREIPH')

* First we need to back-transform the effective Hamiltonian in the
* basis of original CASSCF states by U0 * Heff * U0^T
* Note that in the case of a normal MS-CASPT2 this and the next step
* do not have any effect on Heff and Ueff
      call mma_allocate(U0transpose,Nstate,Nstate)
      call trnsps(Nstate,Nstate,U0,U0transpose)
      call transmat(Heff,U0transpose,Nstate)
      call mma_deallocate(U0transpose)

* Compute transformation matrix that diagonalizes the effective
* Hamiltonian expressed in the basis of original CASSCF states,
* i.e. simply combine the two transf matrices: Ueff = U0 * Ueff
      call mma_allocate(Utmp,Nstate,Nstate)
      call dgemm_('N','N',Nstate,Nstate,Nstate,
     &             1.0d0,U0,Nstate,Ueff,Nstate,
     &             0.0d0,Utmp,Nstate)
      Ueff=Utmp
      call mma_deallocate(Utmp)

C Not called, if .NOT.IFMIX, then only the new CI coefficients are
C printed, no JOBMIX file is created.
      IF (IFMSCOUP.AND.(ISCF.EQ.0)) THEN
        IF (.NOT.IFMIX) THEN
          IF (IPRGLB.GE.USUAL) CALL PRINT_CI_MIX(Ueff)
          RETURN
        END IF
      END IF
      IF (DOCUMULANT .OR. (.NOT.IFMIX)) RETURN

      IF(IFMSCOUP) THEN
        IF(IPRGLB.GE.USUAL) THEN
          WRITE(6,*)' THE ORIGINAL CI ARRAYS ARE NOW MIXED AS LINEAR'
          WRITE(6,*)' COMBINATIONS, GIVEN BY THE EIGENVECTORS.'
        END IF
      ELSE
        Ueff(1,1)=1.0D0
      END IF

      IF(IPRGLB.GE.USUAL) THEN
        WRITE(6,*)' A NEW JOBIPH FILE NAMED ''JOBMIX'' IS PREPARED.'
        WRITE(6,'(20A4)')('****',I=1,20)
      END IF

      CALL GETMEM('LCI1','ALLO','REAL',LCI1,MXCI)
      CALL GETMEM('LCI2','ALLO','REAL',LCI2,MXCI)
      JOBIPH=15
      CALL DANAME(JOBIPH,refwfn_filename)
      JOBMIX=11
      CALL DANAME(JOBMIX,'JOBMIX')
C IADR15 is already known (it is the table of contents of the
C JOBIPH file). When copying/modifying selected data from JOBIPH
C to JOBMIX, we use the same TOC array, IADR15.
      IAD15=0
      CALL IDAFILE(JOBIPH,2,IADR15,30,IAD15)
      IAD15=0
      CALL IDAFILE(JOBMIX,1,IADR15,30,IAD15)
      IAD15=IADR15(1)
* Initialize WEIGHT() (which is unused) just so detection
* of uninitialized memory does not get its knickers twisted
      CALL DCOPY_(MXROOT,[0.0D0],0,WEIGHT,1)
      CALL WR_RASSCF_INFO(JOBMIX,1,iAd15,
     &                    NACTEL,ISPIN,NSYM,LSYM,
     &                    NFRO,NISH,NASH,NDEL,NBAS,8,
     &                    NAME,LENIN8*MXORB,NCONF,HEADER,144,
     &                    TITLE,4*18*MXTIT,POTNUC,
     &                    LROOTS,NROOTS,IROOT,MXROOT,NRAS1,
     &                    NRAS2,NRAS3,NHOLE1,NELE3,IFQCAN,
     &                    Weight)
* Copy MO coefficients from JOBIPH to JOBMIX
      NCMO=NBSQT
      CALL GETMEM('LCMO','ALLO','REAL',LCMO,NCMO)
      IAD15=IADR15(9)
      CALL DDAFILE(JOBIPH,2,WORK(LCMO),NCMO,IAD15)
      IAD15=IADR15(9)
      CALL DDAFILE(JOBMIX,1,WORK(LCMO),NCMO,IAD15)
* If IFQCAN.EQ.0, there is also an additional CMO set:
      IF(IFQCAN.EQ.0) THEN
        IAD15=IADR15(2)
        CALL DDAFILE(JOBIPH,2,WORK(LCMO),NCMO,IAD15)
        IAD15=IADR15(2)
        CALL DDAFILE(JOBMIX,1,WORK(LCMO),NCMO,IAD15)
      END IF
      CALL GETMEM('LCMO','FREE','REAL',LCMO,NCMO)
* Copy all CI coefficients
      IDR=IADR15(4)
      IDW=IADR15(4)
      DO I=1,LROOTS
      CALL DDAFILE(JOBIPH,2,WORK(LCI1),NCONF,IDR)
      CALL DDAFILE(JOBMIX,1,WORK(LCI1),NCONF,IDW)
      END DO
* Replace old energy array with (MS-)CASPT2 energy values:
      NOLDE=MXROOT*MXITER
      CALL GETMEM('OLDE','ALLO','REAL',LOLDE,NOLDE)
      CALL DCOPY_(NOLDE,[0.0D0],0,WORK(LOLDE),1)
*      CALL DCOPY_(NSTATE,ENERGY,1,WORK(LOLDE),1)
      DO ISTATE=1,NSTATE
        ISNUM=MSTATE(ISTATE)
        WORK(LOLDE-1+ISNUM)=ENERGY(ISTATE)
      END DO
      IAD15=IADR15(6)
      CALL DDAFILE(JOBMIX,1,WORK(LOLDE),NOLDE,IAD15)
      CALL GETMEM('OLDE','Free','REAL',LOLDE,NOLDE)
      IAD15=IADR15(18)
CSVC: translates levels to orbital index
      CALL IDAFILE(JOBMIX,1,L2ACT,mxAct,IAD15)
CSVC: translates orbital index to levels
      CALL IDAFILE(JOBMIX,1,LEVEL,mxAct,IAD15)

* PAM07: Eliminate unsafe IPOSFILE calls, use instead dummy i/o operations
* to find disk addresses to CI arrays:
      NIDIST=0
      DO JSTATE=1,NSTATE
        JSNUM=MSTATE(JSTATE)
        NIDIST=MAX(NIDIST,JSNUM)
      END DO
      CALL GETMEM('DIST','ALLO','INTE',LIDIST,NIDIST)
      ID=IADR15(4)
      DO JSNUM=1,NIDIST
        IWORK(LIDIST-1+JSNUM)=ID
* This dummy operation does nothing, merely updates file pointer ID
        CALL DDAFILE(JOBIPH,0,WORK(LCI1),NCONF,ID)
      END DO
* PAM07: Now IDIST() is used, instead of IPOSFILE, below!

C PAM05: Now CREIPH is called also in the NOMULT=1 case, to allow making
C a JOBMIX file also when NOMULT was ordered. Then the energies
C will be state-specific, of course. But no mixing of CI vectors.
      IF (IFMSCOUP) THEN
C Also write effective Hamiltonian on Jobiph file:
        CALL GETMEM('EFFCP','ALLO','REAL',LEFFCP,LROOTS**2)
C Read the effective Hamiltonian on JobIph file:
        IAD15=IADR15(17)
        CALL DDAFILE(JOBIPH,2,WORK(LEFFCP),LROOTS**2,IAD15)
C Replace the relevant elements:
        DO I=1,NSTATE
          ISNUM=MSTATE(I)
          DO J=1,NSTATE
            JSNUM=MSTATE(J)
            WORK(LEFFCP-1+ISNUM+LROOTS*(JSNUM-1))=HEFF(I,J)
          END DO
        END DO
C Write the present effective Hamiltonian:
        IAD15=IADR15(17)
        CALL DDAFILE(JOBIPH,1,WORK(LEFFCP),LROOTS**2,IAD15)
C Write a diagonal Hamiltonian in the JOBMIX:
        IAD15=IADR15(17)
        CALL DCOPY_(LROOTS**2,[0.0D0],0,WORK(LEFFCP),1)
        DO ISTATE=1,NSTATE
          ISNUM=MSTATE(ISTATE)
          WORK(LEFFCP+(ISNUM-1)*LROOTS+ISNUM-1)=ENERGY(ISTATE)
        END DO
        CALL DDAFILE(JOBMIX,1,WORK(LEFFCP),LROOTS**2,IAD15)
        CALL GETMEM('EFFCP','FREE','REAL',LEFFCP,LROOTS**2)
* Now 'mix' those states that were treated in the multi-state CASPT2
        IF (IPRGLB.GE.USUAL) THEN
          WRITE(6,*)
          CALL CollapseOutput(1,'Mixed CI coefficients:')
        END IF
        DO ISTATE=1,NSTATE
          ISNUM=MSTATE(ISTATE)
          CALL DCOPY_(MXCI,[0.0D0],0,WORK(LCI2),1)
          DO JSTATE=1,NSTATE
            JSNUM=MSTATE(JSTATE)
*PAM07           IDISK=IADR15(4)+iPosFile(NCONF)*(JSNUM-1)
            IDISK=IWORK(LIDIST-1+JSNUM)
            CALL DDAFILE(JOBIPH,2,WORK(LCI1),NCONF,IDISK)
            X=Ueff(JSTATE,ISTATE)
            CALL DAXPY_(NCONF,X,WORK(LCI1),1,WORK(LCI2),1)
          END DO
          IF(ISCF.EQ.0) THEN
            IF(IPRGLB.GE.USUAL) THEN
              WRITE(6,'(1x,a,i3)')
     &        ' The CI coefficients for the MIXED state nr. ',ISTATE
              CALL PRWF_CP2(LSYM,NCONF,WORK(LCI2),CITHR)
            END IF
          END IF
*PAM07         IDISK=IADR15(4)+iPosFile(NCONF)*(ISNUM-1)
          IDISK=IWORK(LIDIST-1+ISNUM)
          CALL DDAFILE(JOBMIX,1,WORK(LCI2),NCONF,IDISK)
        END DO
        IF(IPRGLB.GE.USUAL) THEN
          CALL CollapseOutput(0,'Mixed CI coefficients:')
          WRITE(6,*)
        END IF
      END IF
      CALL GETMEM('DIST','FREE','INTE',LIDIST,NIDIST)
      CALL GETMEM('LCI1','FREE','REAL',LCI1,MXCI)
      CALL GETMEM('LCI2','FREE','REAL',LCI2,MXCI)

      CALL DACLOS(JOBIPH)
      CALL DACLOS(JOBMIX)

      CALL QEXIT('CREIPH')

      RETURN
      END
