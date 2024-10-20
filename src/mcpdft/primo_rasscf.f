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
      Subroutine Primo_RasScf_m(VecTit,Ene,Occ,CMO)
************************************************************************
* purpose:
* Print MO-coefficients
*
************************************************************************

      use mcpdft_output, only: lf
      use stdalloc, only: mma_allocate, mma_deallocate

      Implicit Real*8 (A-H,O-Z)

      Character(LEN=*) VecTit
      Real*8 CMO(*),Occ(*),Ene(*)

#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"

      Integer NSLCT(8)
      Character(LEN=3) lIrrep(8)

      Character(LEN=8) Fmt1,Fmt2
      Character(LEN=132) Line,Blank
      Character(LEN=LENIN8), External:: Clean_BName

      Integer, Allocatable:: MrkIt(:), SLCT(:)

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
* determined in the CHKINP subroutine


* PAM Nov 09: Output marked as collapsible.
* print header
      Write(LF,*)
      Call CollapseOutput(1,'   Molecular orbitals:')
      Write(LF,'(6X,A)') '-------------------'
      Write(LF,*)
      Write(LF,Fmt2//'A)') trim(VecTit)
      Write(LF,*)

* Flag ipt2 in common in src/Include/rasscf.fh
*   ipt2=0 means usual MO's, quasicanonical for
! In MC-PDFT, ipt2 is always 0
* inactives and virtuals, natural for active.
* PAM Apr 05: The rules for selecting orbitals for print
* appear a bit confused. I just follow the rules for now:



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
      CALL mma_allocate(MRKIT,NBTOT,Label='MrkIt')
      IORB=0
      DO ISYM=1,NSYM
       NB=NBAS(ISYM)
       DO I=1,NB
        IORB=IORB+1
        MRKIT(IORB)=0
       END DO
      END DO
* Mark MARKIT as selected for occupied orbitals.
      IORB=0
      DO ISYM=1,NSYM
       NFIA=NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
       DO I=1,NFIA
        IORB=IORB+1
        MRKIT(IORB)=1
       END DO
       IORB=IORB+NBAS(ISYM)-NFIA
      END DO

* only those orbitals that have occ no larger
* than or equal to PrOThr:
       IORB=0
       DO ISYM=1,NSYM
        NFIA=NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
        DO I=1,NFIA
         IORB=IORB+1
         ! TODO(matthew hennefarth) : if this is occupancy, then the
         ! following line will never occur
         IF(OCC(IORB).lt. 0.0D0) MRKIT(IORB)=0
        END DO
        IORB=IORB+NBAS(ISYM)-NFIA
       END DO
* also those orbitals that have energy less
* than or equal to PrEThr, skipping deleted orbitals of course.
       IORB=0
       DO ISYM=1,NSYM
        NFIA=NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
        IORB=IORB+NFIA
        ND=NDEL(ISYM)
        DO I=NFIA+1,NBAS(ISYM)-ND
         IORB=IORB+1
         IF(ENE(IORB).le.PRETHR) MRKIT(IORB)=1
        END DO
        IORB=IORB+ND
       END DO

* Let ISELECT enumerate the orbitals to be printed, rather than
* just marking them:
      NSLCTT=0
      NO=0
      DO ISYM=1,NSYM
       NS=0
       NB=NBAS(ISYM)
       DO IO=NO+1,NO+NB
         IF(MRKIT(IO).EQ.1) NS=NS+1
       END DO
       NO=NO+NB
       NSLCT(ISYM)=NS
       NSLCTT=NSLCTT+NS
      END DO
      CALL mma_allocate(SLCT,NSLCTT,Label='SLCT')
      IS=0
      NO=0
      DO ISYM=1,NSYM
       NB=NBAS(ISYM)
       DO IO=NO+1,NO+NB
        IF(MRKIT(IO).EQ.1) THEN
          IS=IS+1
          SLCT(IS)=IO
        END IF
       END DO
       NO=NO+NB
      END DO
* Get rid of MARKIT.
      CALL mma_deallocate(MRKIT)

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
           Write(LF,Fmt2//'A,6X,10I10)')'Orbital ',
     &               (SLCT(ISOFF+I)-IBOFF,I=ISSTART,ISEND)
           Write(LF,Fmt2//'A,6X,10F10.4)')'Energy  ',
     &               (ENE(SLCT(ISOFF+I)),I=ISSTART,ISEND)
           Write(LF,Fmt2//'A,6X,10F10.4)')'Occ. No.',
     &               (OCC(SLCT(ISOFF+I)),I=ISSTART,ISEND)
           Write(LF,*)
           DO IB=1,NB
            Write(LF,'(2X,I3,1X,A,10F10.4)') IB,
     &        Clean_BName(NAME(IBOFF+IB),LENIN),
     &        (CMO(ICOFF+(SLCT(ISOFF+I)-1-IBOFF)*NB+IB),
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
            Write(LF,FMT2//'A)')
     &          'INDEX  ENERGY  OCCUPATION COEFFICIENTS ...'

            DO IS = 1,NSLCT(ISYM)
              IORB=SLCT(ISOFF+IS)
              ICOL=IORB-IBOFF
              LINE = BLANK
              IST = 1
              Write(LINE(IST:132),'(I5)') ICOL
              IST = IST+5
              Write(LINE(IST:132),'(F10.4)') ENE(IORB)
               IST = IST+10
              Write(LINE(IST:132),'(F10.4)') OCC(IORB)
              IST = IST+10
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

      CALL mma_deallocate(SLCT)
* PAM 09: Reset to non-collapsible output before return.
      Call CollapseOutput(0,'   Molecular orbitals:')
      Write(6,*)
*
      Return
      End
