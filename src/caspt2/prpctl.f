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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE PRPCTL
      USE PT2WFN
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "SysDef.fh"
      Logical FullMlk,lSave,Do_ESPF
#ifdef _MOLCAS_MPP_
      LOGICAL Is_Real_Par
#endif

      Character(Len=8) Label
      Character(Len=128) FILENAME,MDNAME
      Character(Len=80) Note
      Integer IndType(56)
      Real*8 Dummy(2),DUM(1)

      CALL QENTER('PRPCTL')

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        WRITE(6,'(1X,A)') ' ====================================='
        WRITE(6,'(1X,A)') ' CASPT2 properties were requested, but'
        WRITE(6,'(1X,A)') ' these are not efficiently implemented'
        WRITE(6,'(1X,A)') ' in parallel. If you do not need them,'
        WRITE(6,'(1X,A)') ' not using the PROP keyword could lead'
        WRITE(6,'(1X,A)') ' to a significant speed up.'
        WRITE(6,'(1X,A)') ' ====================================='
      END IF
#endif

* PAM2008 When this subroutine is called, the calculation has been done
* for the (individual) state nr JSTATE in 1,2,..,NSTATE.
* The corresponding CI-root from rasscf, the root state for this PT2,
* is number MSTATE(JSTATE) on the input JOBIPH file.
* JSTATE,NSTATE and MSTATE() are in common in the caspt2.fh file.
      IERR=0
      IF(NSTATE.GT.1) THEN
       N=MSTATE(JSTATE)
       IF(N.LE.0 .or. N.GT.999) THEN
        WRITE(6,*)' Subroutine PRPCTL fails -- It seems to get data'
        WRITE(6,*)' computed for a root nr ',N
        WRITE(6,*)' which is surely wrong.'
        WRITE(6,*)' PRPCTL gives up, there will be no calculations'
        WRITE(6,*)' done of orbitals, properties, etc for this state.'
        WRITE(6,*)' This was state nr JSTATE=',JSTATE
        WRITE(6,*)' in the MS-CASPT2 calculation.'
        IERR=1
       END IF
      END IF
      IF(IERR.GT.0) GOTO 999

      FullMlk=.True.
      IF (.NOT.PRORB ) FullMlk=.False.

      IF ( IPRGLB.GE.USUAL ) THEN
       WRITE(6,'(20A4)')('----',I=1,20)
      END IF

C Compute density matrix, output orbitals, and properties.

C Compute density matrix of CASPT2 wave function, in MO basis,
C to produce output orbitals.
C This density matrix may be approximated in several ways, see DENS.
      NDMAT=0
      NOCC=0
      DO ISYM=1,NSYM
        NO=NORB(ISYM)
        NDMAT=NDMAT+(NO**2+NO)/2
        NOCC=NOCC+NBAS(ISYM)
      END DO
      CALL GETMEM('DMAT','ALLO','REAL',LDMAT,NDMAT)
      CALL DCOPY_(NDMAT,[0.0D0],0,WORK(LDMAT),1)
      CALL GETMEM('LISTS','ALLO','INTE',LLISTS,NLSTOT)
      CALL MKLIST(iWORK(LLISTS))
      CALL DENS(IVECX,WORK(LDMAT))
      CALL GETMEM('LISTS','FREE','INTE',LLISTS,NLSTOT)

C Compute natural orbitals of CASPT2 wave function.
      CALL GETMEM('CMO','ALLO','REAL',LCMO,NCMO)
      CALL DCOPY_(NCMO,WORK(LCMOPT2),1,WORK(LCMO),1)
      CALL GETMEM('CNAT','ALLO','REAL',LCNAT,NCMO)
      CALL GETMEM('OCC','ALLO','REAL',LOCC,NOCC)
      CALL NATORB_CASPT2(WORK(LDMAT),WORK(LCMO),WORK(LOCC),WORK(LCNAT))
      CALL GETMEM('LCMO','FREE','REAL',LCMO,NCMO)
C Backtransform density matrix to original MO basis before storing
      CALL TRANSFOCK(WORK(LTORB),WORK(LDMAT),-1)
      CALL PT2WFN_DENSSTORE(WORK(LDMAT),NDMAT)
      CALL GETMEM('DMAT','FREE','REAL',LDMAT,NDMAT)
C Write natural orbitals as standard orbital file on LUMORB.
* PAM2008: Separate PT2ORB files for each state:
      FILENAME='PT2ORB'
      MDNAME='MD_PT2'
      IF(NSTATE.GT.1) THEN
       FILENAME='PT2ORB.x'
       MDNAME='MD_PT2.x'
       N=MSTATE(JSTATE)
       IF(N.LE.9) THEN
        WRITE(FILENAME(8:8),'(I1)') N
        WRITE(MDNAME(8:8),'(I1)') N
       ELSE IF(N.LE.99) THEN
        WRITE(FILENAME(8:9),'(I2)') N
        WRITE(MDNAME(8:9),'(I2)') N
       ELSE IF(N.LE.999) THEN
        WRITE(FILENAME(8:10),'(I3)') N
        WRITE(MDNAME(8:10),'(I3)') N
       END IF
      END IF
* PAM2008: For MS-CASPT2 with more than one state, orbital files will
* now be numbered PT2ORB.1 ... PT2ORB.999
* depending on which CI-root of the rasscf calculation that was
* the root state of the PT2.
      LUTMP=19
      LUTMP=ISFREEUNIT(LUTMP)
* PAM 2008: Add typeindex information:
*----------------------------------------------------------------------*
*     Make typeindex information                                       *
*----------------------------------------------------------------------*
      iShift=0
      DO ISYM=1,NSYM
        IndT=0
        IndType(1+iShift)= NFRO(ISYM)
        IndT=IndT+NFRO(ISYM)
        IndType(2+iShift)= NISH(ISYM)
        IndT=IndT+NISH(ISYM)
        IndType(3+iShift)= NRAS1(ISYM)
        IndT=IndT+NRAS1(ISYM)
        IndType(4+iShift)= NRAS2(ISYM)
        IndT=IndT+NRAS2(ISYM)
        IndType(5+iShift)= NRAS3(ISYM)
        IndT=IndT+NRAS3(ISYM)
        IndType(7+iShift)= NDEL(ISYM)
        IndT=IndT+NDEL(ISYM)
        IndType(6+iShift)= NBAS(ISYM)-IndT
        iShift=iShift+7
      EndDo
      If (NSTATE.GT.1) THEN
        Write(Note,'(A41,I3,A3,f22.12)')
     &   '* CASPT2 natural orbitals for root number',N,
     &   ' E=',Energy(JSTATE)
      Else
        Note='* CASPT2 natural orbitals'
      End If

      CALL WRVEC(FILENAME,LUTMP,'COI',NSYM,NBAS,NBAS,
     &  WORK(LCNAT), WORK(LOCC),Dummy  ,IndType,Note)
      iUHF=0
      Call Molden_Interface(iUHF,FILENAME,MDNAME)

C Write natural orbitals to standard output.
      IF ( IPRGLB.GE.VERBOSE) THEN
       WRITE(6,*)
       WRITE(6,'(A)')'  The CASPT2 orbitals are computed as natural '//
     &          'orbitals of a density matrix'
       WRITE(6,'(A)')'  defined as:'
       WRITE(6,'(A)')'   D = (D0 + D1 + D2)/<PSI|PSI>'
       WRITE(6,'(A)')' where D0..D2 are zeroth..2nd order contributions'
       WRITE(6,'(A)')' and |PSI> is the total wave function.'
       WRITE(6,'(A)')' A new RasOrb file named PT2ORB is prepared.'
       IF (PRORB) THEN
         IF ( OUTFMT.EQ.'LONG    ' ) THEN
           THRENE=2.0d0**31
           THROCC=-2.0d0**31
         ELSE IF ( OUTFMT.EQ.'DEFAULT ' ) THEN
           THRENE=5.0d+00
           THROCC=5.0d-04
         END IF
         CALL PRIMO('Output orbitals from CASPT2',
     &           .TRUE.,.FALSE.,THROCC,THRENE,NSYM,NBAS,
     &            NBAS,NAME,DUM,WORK(LOCC),WORK(LCNAT),-1)
       END IF
      END IF

* compute Mulliken's orbital populations

      IF ( IPRGLB.GE.USUAL ) THEN
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,'(6X,A)') 'Mulliken population Analysis:'
        WRITE(6,'(6X,A)') '-----------------------------'

        Call GetMem('Scr1','Allo','Real',LXXX,NBAST**2)
        iRc=-1
        iOpt=6
        iComp=1
        iSyLbl=1
        Label='Mltpl  0'
        Call RdOne(iRc,iOpt,Label,iComp,Work(LXXX),iSyLbl)
        If ( iRc.eq.0 ) then
           lSave = MSTATE(JSTATE) .eq. irlxroot
           Call Charge(nSym,nBas,Name,
     &               WORK(LCNAT),WORK(LOCC),WORK(LXXX),2,FullMlk,lSave)
        End If
        Call GetMem('Scr1','Free','Real',lXXX,NBAST**2)
      END IF

* compute one-electron properties

      IF ( IPRGLB.GE.USUAL ) THEN
        WRITE(6,*)
        WRITE(6,'(6X,A)') 'Expectation values of various properties:'
        WRITE(6,'(6X,A)') '-----------------------------------------'
      END IF

* The PRPT source code gives the following formula for the
* scratch space needed:
      NCOMP=6
      NTCOMP=15
      NSCR=(NBSQT+NBAST)/2+6+4*NCOMP+(NBAST*(NBAST+1))/2
     &      +4+2*NTCOMP*(NTCOMP+1)

      nDens=0
      Do i = 1, nSym
         nDens=nDens+nBas(i)*(nBas(i)+1)/2
      End Do
      Call GetMem('Scr2','Allo','Real',LXXX,NDENS)
*
      Call DOne_Caspt2(WORK(LCNAT),WORK(LOCC),Work(LXXX))
      Call Put_D1AO(WORK(LXXX),nDens)
*
      Note='Temporary orbital file used by prpt.'
      LuTmp=50
      LuTmp=IsFreeUnit(LuTmp)
      Call WrVec('TMPORB',LuTmp,'CO',nSym,nBas,nBas,
     &            WORK(LCNAT),WORK(LOCC),Dummy,IndType,Note)
      Call Prpt()
cnf
*
*------- ESPF analysis
      Call DecideOnESPF(Do_ESPF)
      lSave = MSTATE(JSTATE) .eq. irlxroot
      If (Do_ESPF) Call espf_analysis(lSave)
cnf
*
*---- On return from PrPt the 1-particle matrix is stored
*     in the beginning of the scratch array.
*
      Call GetMem('Scr2','Free','Real',lXXX,NDENS)
*
      CALL GETMEM('OCC','FREE','REAL',LOCC,NOCC)
      CALL GETMEM('CNAT','FREE','REAL',LCNAT,NCMO)

 999  CONTINUE

      CALL QEXIT('PRPCTL')
      RETURN
      END
