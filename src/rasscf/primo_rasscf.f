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
      Subroutine Primo_RasScf(VecTit,Ene,Occ,CMO)
************************************************************************
* purpose:
* Print MO-coefficients
*
************************************************************************

      Implicit Real*8 (A-H,O-Z)

      Character*(*) VecTit
      Dimension     CMO(*),Occ(*),Ene(*)

#include "rasdim.fh"
#include "general.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='PRIMO   ')
#include "rasscf.fh"
#include "WrkSpc.fh"

      DIMENSION NSLCT(8)
      Logical   PrOcc,PrEne
      Character*3 lIrrep(8)

      Character*8 Fmt1,Fmt2
      Character*132 Line,Blank
      Character*(LENIN8) Clean_BName
      External Clean_BName
* PAM Nov 05: Non-valence orbitals
      Dimension NVSH(8)

      Call Get_cArray('Irreps',lIrrep,24)

* set print format
      nCols=10
      lLine=120
      lPaper=132
      left=(lPaper-lLine)/2
      Write(Fmt1,'(A,I3.3,A)') '(',left,'X,A)'
      Write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
      Blank = ' '
* PAM Krapperup Nov 05: For the moment, selection of orbitals to be
* printed is ultimately determined here on the basis of PRETHR (and
* PROTHR) thresholds. These have either been set by the user, or
* determined in the CHKINP subroutine, possibly then on basis of
* user specification of OutFmt1 flags.
* Exception: OutFmt2='NOCORE  ' will inhibit printing of non-valece orbs.


* PAM Nov 09: Output marked as collapsible.
* print header
      Write(LF,*)
      Call CollapseOutput(1,'   Molecular orbitals:')
      Write(LF,'(6X,A)') '-------------------'
      Write(LF,*)
      Write(LF,Fmt2//'A)') VecTit(:mylen(VecTit))
      Write(LF,*)

* Flag ipt2 in common in src/Include/rasscf.fh
*   ipt2=0 means usual MO's, quasicanonical for
* inactives and virtuals, natural for active.
*   ipt2=1 means quasicanonical for actives also.
* PAM Apr 05: The rules for selecting orbitals for print
* appear a bit confused. I just follow the rules for now:
      If ( iPT2.eq.0 ) then
        PrOcc  = .true.
        PrEne  = .true.
      Else
        PrOcc  = .false.
        PrEne  = .true.
      End if

* Select orbitals to be printed.
* By default at least all occupied are printed.
* If iPT2==0 all orbitals with energy less than PrEThr
* and occupation bigger then PrOThr are also printed.
* If iPT2==1 all orbitals with occupations greater than
* PrOThr are also printed.

* Nr of orbitals=nr of basis functions
      NBTOT=0
      DO ISYM=1,NSYM
       NBTOT=NBTOT+NBAS(ISYM)
      END DO
*
* Put rasscf orbital energies on the runfile
      Call Put_darray('RASSCF OrbE',ENE,NBTOT)
*
* Initialize MARKIT.
      CALL GETMEM('MARKIT','ALLO','INTE',LMRKIT,NBTOT)
      IORB=0
      DO ISYM=1,NSYM
       NB=NBAS(ISYM)
       DO I=1,NB
        IORB=IORB+1
        IWORK(LMRKIT-1+IORB)=0
       END DO
      END DO
      IF (OutFmt2.ne.'NOCORE  ') THEN
* Mark MARKIT as selected for occupied orbitals.
       IORB=0
       DO ISYM=1,NSYM
        NFIA=NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
        DO I=1,NFIA
         IORB=IORB+1
         IWORK(LMRKIT-1+IORB)=1
        END DO
        IORB=IORB+NBAS(ISYM)-NFIA
       END DO
      ELSE
* Mark MARKIT as selected for valence or active orbitals.
       IORB=0
       Call Get_iArray('Non valence orbitals',NVSH,nSym)
       DO ISYM=1,NSYM
        NFIA=NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
        NSKIP=MIN(NFRO(ISYM)+NISH(ISYM),NVSH(ISYM))
        IORB=IORB+NSKIP
        DO I=NSKIP+1,NFIA
         IORB=IORB+1
         IWORK(LMRKIT-1+IORB)=1
        END DO
        IORB=IORB+NBAS(ISYM)-NFIA
       END DO
      END IF
* If PROCC, then only those orbitals that have occ no larger
* than or equal to PrOThr:
      IF (PROCC.and. PrOThr.ge.0.0D0) THEN
       IORB=0
       DO ISYM=1,NSYM
        NFIA=NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
        DO I=1,NFIA
         IORB=IORB+1
         IF(OCC(IORB).lt.PROTHR) IWORK(LMRKIT-1+IORB)=0
        END DO
        IORB=IORB+NBAS(ISYM)-NFIA
       END DO
      END IF
* But if PRENE, then also those orbitals that have energy less
* than or equal to PrEThr, skipping deleted orbitals of course.
      IF (PRENE) THEN
       IORB=0
       DO ISYM=1,NSYM
        NFIA=NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
        IORB=IORB+NFIA
        ND=NDEL(ISYM)
        DO I=NFIA+1,NBAS(ISYM)-ND
         IORB=IORB+1
         IF(ENE(IORB).le.PRETHR) IWORK(LMRKIT-1+IORB)=1
        END DO
        IORB=IORB+ND
       END DO
      END IF
* Let ISELECT enumerate the orbitals to be printed, rather than
* just marking them:
      NSLCTT=0
      NO=0
      DO ISYM=1,NSYM
       NS=0
       NB=NBAS(ISYM)
       DO IO=NO+1,NO+NB
         IF(IWORK(LMRKIT-1+IO).EQ.1) NS=NS+1
       END DO
       NO=NO+NB
       NSLCT(ISYM)=NS
       NSLCTT=NSLCTT+NS
      END DO
      CALL GETMEM('SELECT','ALLO','INTE',LSLCT,NSLCTT)
      IS=0
      NO=0
      DO ISYM=1,NSYM
       NB=NBAS(ISYM)
       DO IO=NO+1,NO+NB
        IF(IWORK(LMRKIT-1+IO).EQ.1) THEN
          IS=IS+1
          IWORK(LSLCT-1+IS)=IO
        END IF
       END DO
       NO=NO+NB
      END DO
* Get rid of MARKIT.
      CALL GETMEM('MARKIT','FREE','INTE',LMRKIT,NBTOT)

* finally, print, the MOs
      If ( OutFmt2.eq.'FULL    ' ) then

* print orbitals, using default format.
        ISOFF=0
        IBOFF=0
        ICOFF=0
        DO ISYM=1,NSYM
         NB=NBAS(ISYM)
         NS=NSLCT(ISYM)
         IF(NS.GT.0) THEN
          Write(LF,*)
          Write(LF,*)
          Write(LF,*)
          Write(LF,Fmt2//'A,I2,A,A)')
     &           'Molecular orbitals for symmetry species',iSym,
     &           ': ',lIrrep(iSym)
          DO ISSTART=1,NS,NCOLS
           ISEND=MIN(ISSTART+NCOLS-1,NS)
           Write(LF,*)
           Write(LF,*)
           Write(LF,Fmt2//'A,7X,10I10)')'Orbital ',
     &               (IWORK(LSLCT-1+ISOFF+I)-IBOFF,I=ISSTART,ISEND)
           IF (PRENE) THEN
             Write(LF,Fmt2//'A,7X,10F10.4)')'Energy  ',
     &               (ENE(IWORK(LSLCT-1+ISOFF+I)),I=ISSTART,ISEND)
           END IF
           IF (PROCC) THEN
             Write(LF,Fmt2//'A,7X,10F10.4)')'Occ. No.',
     &               (OCC(IWORK(LSLCT-1+ISOFF+I)),I=ISSTART,ISEND)
           END IF
           Write(LF,*)
           DO IB=1,NB
            Write(LF,'(2X,I4,1X,A,10F10.4)') IB,
     &        Clean_BName(NAME(IBOFF+IB),LENIN),
     &        (CMO(ICOFF+(IWORK(LSLCT-1+ISOFF+I)-1-IBOFF)*NB+IB),
     &        I=ISSTART,ISEND)
           END DO
          END DO
         END IF
         ISOFF=ISOFF+NS
         IBOFF=IBOFF+NB
         ICOFF=ICOFF+NB**2
        END DO

      Else If ( OutFmt2.eq.'COMPACT ' ) then

* print orbitals, using compact format.

        ISOFF=0
        IBOFF=0
        ICOFF=0
        DO ISYM=1,NSYM
          NB=NBAS(ISYM)
          IF ( NSLCT(ISYM).NE.0 ) THEN
            Write(LF,*)
            Write(LF,FMT2//'A,I2,A,A)')
     &           'MOLECULAR ORBITALS FOR SYMMETRY SPECIES',ISYM,
     &           ': ',LIRREP(ISYM)
            Write(LF,*)
            IF ( PROCC.AND.PRENE ) THEN
              Write(LF,FMT2//'A)')
     &          'INDEX  ENERGY  OCCUPATION COEFFICIENTS ...'
            ELSE IF ( PROCC ) THEN
              Write(LF,FMT2//'A)')
     &          'INDEX  ENERGY  COEFFICIENTS ...'
            ELSE IF ( PRENE ) THEN
              Write(LF,FMT2//'A)')
     &          'INDEX  OCCUPATION  COEFFICIENTS ...'
            ELSE
              Write(LF,FMT2//'A)')
     &          'INDEX  COEFFICIENTS ...'
            END IF
            DO IS = 1,NSLCT(ISYM)
              IORB=IWORK(LSLCT-1+ISOFF+IS)
              ICOL=IORB-IBOFF
              LINE = BLANK
              IST = 1
              Write(LINE(IST:132),'(I5)') ICOL
              IST = IST+5
              IF ( PRENE ) THEN
                Write(LINE(IST:132),'(F10.4)') ENE(IORB)
                 IST = IST+10
              END IF
              IF ( PROCC ) THEN
                Write(LINE(IST:132),'(F10.4)') OCC(IORB)
                IST = IST+10
              END IF
              Write(LF,FMT2//'A)') LINE
              LINE = BLANK
              IST = 9
              DO IBAS = 1,NB
                CC = CMO(ICOFF+(ICOL-1)*NB+IBAS)
                IF ( ABS(CC).GE.0.1D0 ) THEN
                  Write(LINE(IST:132),'(I4,1X,A,A,F7.4,A)')
     &              IBAS,Clean_BName(NAME(IBOFF+IBAS),LENIN),
     &              '(',CC,')'
                  IST = IST+28
                  IF ( IST.GT.(132-LEFT-28) ) THEN
                    Write(LF,FMT2//'A)') TRIM(LINE)
                    LINE = BLANK
                    IST = 9
                  END IF
                END IF
              END DO
              Write(LF,FMT2//'A)') TRIM(LINE)
              LINE = BLANK
              IST = 9
            END DO
          END IF
          ISOFF=ISOFF+NSLCT(ISYM)
          IBOFF=IBOFF+NB
          ICOFF=ICOFF+NB**2
        END DO

      End If

      CALL GETMEM('SELECT','FREE','INTE',LSLCT,NSLCTT)
* PAM 09: Reset to non-collapsible output before return.
      Call CollapseOutput(0,'   Molecular orbitals:')
      Write(6,*)
*
      Return
      End
