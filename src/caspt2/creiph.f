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
      use output, only:silent,terse,usual,verbose,debug,insane,iPrGlb
      USE REFWFN, ONLY: REFWFN_FILENAME, IADR15
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
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
      integer JOBIPH, JOBMIX
      real(8) Weight(MxRoot)
      real(8) Heff(Nstate,Nstate),Ueff(Nstate,Nstate),U0(Nstate,Nstate)


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
      END IF

      IF(IPRGLB.GE.USUAL) THEN
        WRITE(6,*)' A NEW JOBIPH FILE NAMED ''JOBMIX'' IS PREPARED.'
        WRITE(6,'(20A4)')('****',I=1,20)
      END IF

* Note that JOBIPH file will contain all the RASSCF CI vectors
* plus a CASPT2 effective Hamiltonian for the selected states.
* The effective Hamiltonian for the states not included in the
* CASPT2 treatment will be diagonal with RASSCF energies!

* The JOBMIX will contain the (possibly mixed) CI vectors,
* with CASPT2 energies for the selected states and zero energy
* for states not included in the CASPT2 treatment
* If NoMulti was specified, the original state indexing is
* maintained, otherwise the new states are just 1, 2, 3...

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
* Modify root index in case of MS
      CALL GETMEM('JROOT','ALLO','INTE',LJROOT,MXROOT)
      IF (IFMSCOUP) THEN
        CALL ICOPY(MXROOT,[0],0,IWORK(LJROOT),1)
        DO ISTATE=1,NSTATE
          IWORK(LJROOT+ISTATE-1)=ISTATE
        END DO
        MROOTS=NSTATE
      ELSE
        CALL ICOPY(MXROOT,IROOT,1,IWORK(LJROOT),1)
        MROOTS=NROOTS
      END IF
* Initialize WEIGHT() (which is unused) just so detection
* of uninitialized memory does not get its knickers twisted
      CALL DCOPY_(MXROOT,[0.0D0],0,WEIGHT,1)
      CALL WR_RASSCF_INFO(JOBMIX,1,iAd15,
     &                    NACTEL,ISPIN,NSYM,STSYM,
     &                    NFRO,NISH,NASH,NDEL,NBAS,8,
     &                    NAME,LENIN8*MXORB,NCONF,HEADER,144,
     &                    TITLE,4*18*MXTIT,POTNUC,
     &                    LROOTS,MROOTS,IWORK(LJROOT),MXROOT,NRAS1,
     &                    NRAS2,NRAS3,NHOLE1,NELE3,IFQCAN,
     &                    Weight)
      CALL GETMEM('JROOT','FREE','INTE',LJROOT,MXROOT)
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
      IF (IFMSCOUP) THEN
        CALL DCOPY_(NSTATE,ENERGY,1,WORK(LOLDE),1)
      ELSE
        DO ISTATE=1,NSTATE
          ISNUM=MSTATE(ISTATE)
          WORK(LOLDE-1+ISNUM)=ENERGY(ISTATE)
        END DO
      END IF
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
      DO ISTATE=1,NSTATE
        JSNUM=MSTATE(ISTATE)
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
          WORK(LEFFCP+(ISTATE-1)*LROOTS+ISTATE-1)=ENERGY(ISTATE)
        END DO
        CALL DDAFILE(JOBMIX,1,WORK(LEFFCP),LROOTS**2,IAD15)
        CALL GETMEM('EFFCP','FREE','REAL',LEFFCP,LROOTS**2)
* Now 'mix' those states that were treated in the multi-state CASPT2
        IF (IPRGLB.GE.USUAL) THEN
          WRITE(6,*)
          CALL CollapseOutput(1,'Mixed CI coefficients:')
        END IF
        DO ISTATE=1,NSTATE
          CALL DCOPY_(MXCI,[0.0D0],0,WORK(LCI2),1)
          DO IISTATE=1,NSTATE
            JSNUM=MSTATE(IISTATE)
            IDISK=IWORK(LIDIST-1+JSNUM)
            CALL DDAFILE(JOBIPH,2,WORK(LCI1),NCONF,IDISK)
            X=Ueff(IISTATE,ISTATE)
            CALL DAXPY_(NCONF,X,WORK(LCI1),1,WORK(LCI2),1)
          END DO
          IF(ISCF.EQ.0) THEN
            IF(IPRGLB.GE.USUAL) THEN
              WRITE(6,'(1x,a,i3)')
     &        ' The CI coefficients for the MIXED state nr. ',ISTATE
              CALL PRWF_CP2(STSYM,NCONF,WORK(LCI2),CITHR)
            END IF
          END IF
          IDISK=IWORK(LIDIST-1+ISTATE)
          CALL DDAFILE(JOBMIX,1,WORK(LCI2),NCONF,IDISK)
        END DO
        IF(IPRGLB.GE.USUAL) THEN
          CALL CollapseOutput(0,'Mixed CI coefficients:')
          WRITE(6,*)
        END IF
C In case of XMS/XDW/RMS and NOMUL, the CI vectors are replaced by the
C rotated zeroth-order states (they should have been printed earlier,
C in grpini)
      ELSE IF (IFXMS.or.IFRMS) THEN
        DO ISTATE=1,NSTATE
          ISNUM=MSTATE(ISTATE)
          CALL DCOPY_(MXCI,[0.0D0],0,WORK(LCI2),1)
          DO IISTATE=1,NSTATE
            JSNUM=MSTATE(IISTATE)
            IDISK=IWORK(LIDIST-1+JSNUM)
            CALL DDAFILE(JOBIPH,2,WORK(LCI1),NCONF,IDISK)
            X=U0(IISTATE,ISTATE)
            CALL DAXPY_(NCONF,X,WORK(LCI1),1,WORK(LCI2),1)
          END DO
          IDISK=IWORK(LIDIST-1+ISNUM)
          CALL DDAFILE(JOBMIX,1,WORK(LCI2),NCONF,IDISK)
        END DO
      END IF
      CALL GETMEM('DIST','FREE','INTE',LIDIST,NIDIST)
      CALL GETMEM('LCI1','FREE','REAL',LCI1,MXCI)
      CALL GETMEM('LCI2','FREE','REAL',LCI2,MXCI)

      CALL DACLOS(JOBIPH)
      CALL DACLOS(JOBMIX)


      RETURN
      END
