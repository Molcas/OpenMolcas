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
      SUBROUTINE MECTL(PROP,OVLP,HAM,ESHFT)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='MECTL')
#include "symmul.fh"
#include "rassi.fh"
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "Files.fh"
#include "SysDef.fh"
      REAL*8 PROP(NSTATE,NSTATE,NPROP),OVLP(NSTATE,NSTATE),
     &       HAM(NSTATE,NSTATE),ESHFT(NSTATE)
*


      CALL QENTER(ROUTINE)

* PAM Sep 2014: this section moved to within GTDMCTL loops
*C Transition density matrices, TDMZZ, in AO basis.
*C WDMZZ similar, but WE-reduced 'triplet' densities.
*      CALL GETMEM('TDMZZ','Allo','Real',LTDMZZ,NTDMZZ)
*      CALL GETMEM('TSDMZZ','Allo','Real',LTSDMZZ,NTDMZZ)
*      CALL GETMEM('WDMZZ','Allo','Real',LWDMZZ,NTDMZZ)
*C Double loop over the states
*      DO ISTATE=1,NSTATE
*        DO JSTATE=1,ISTATE
*C IDTDM: TOC array for transition 1-matrices
*          IDISK=IDTDM(ISTATE,JSTATE)
*          CALL DDAFILE(LUTDM,2,WORK(LTDMZZ),NTDMZZ,IDISK)
*          CALL DDAFILE(LUTDM,2,WORK(LTSDMZZ),NTDMZZ,IDISK)
*          CALL DDAFILE(LUTDM,2,WORK(LWDMZZ),NTDMZZ,IDISK)
*          CALL PROPER (PROP,ISTATE,JSTATE,WORK(LTDMZZ),WORK(LWDMZZ))
*        END DO
*      END DO
*      CALL GETMEM('TDMZZ','Free','Real',LTDMZZ,NTDMZZ)
*      CALL GETMEM('TSDMZZ','Free','Real',LTSDMZZ,NTDMZZ)
*      CALL GETMEM('WDMZZ','Free','Real',LWDMZZ,NTDMZZ)


C Print results:
      NCOL=4
      IF(PRXVR.and.IPGLOB.ge.SILENT) THEN
      WRITE(6,*)
      Call CollapseOutput(1,'Expectation values for input states')
      WRITE(6,'(3X,A)')     '-----------------------------------'
      WRITE(6,*)
      WRITE(6,*)' EXPECTATION VALUES OF 1-ELECTRON OPERATORS'
      WRITE(6,*)' FOR THE RASSCF INPUT WAVE FUNCTIONS:'
      WRITE(6,*)
      DO IPROP=1,NPROP
       IF(IPUSED(IPROP).NE.0) THEN

* Skip printing if all the diagonal values are very small
*  (presumed zero for reasons of selection rules)
        PLIMIT=1.0D-10
        PMAX=0.0D0
        DO I=1,NSTATE
         PMAX=MAX(PMAX,ABS(PROP(I,I,IPROP)+PNUC(IPROP)*OVLP(I,I)))
        END DO
        IF(PMAX.GE.PLIMIT) THEN

        DO ISTA=1,NSTATE,NCOL
          IEND=MIN(NSTATE,ISTA+NCOL-1)
          WRITE(6,*)
        WRITE(6,'(1X,A,A8,A,I4)')
     &  'PROPERTY: ',PNAME(IPROP),'   COMPONENT:**',ICOMP(IPROP)
          WRITE(6,'(1X,A,3D17.8)')
     &'ORIGIN    :',(PORIG(I,IPROP),I=1,3)
          WRITE(6,'(1X,A,I8,3I17)')
     &'STATE     :',(I,I=ISTA,IEND)
          WRITE(6,*)
          WRITE(6,'(1X,A,4(1x,G16.9))')
     &'ELECTRONIC:',(PROP(I,I,IPROP),I=ISTA,IEND)
          WRITE(6,'(1X,A,4(1x,G16.9))')
     &'NUCLEAR   :',(PNUC(IPROP)*OVLP(I,I),I=ISTA,IEND)
          WRITE(6,'(1X,A,4(1x,G16.9))')
     &'TOTAL     :',(PROP(I,I,IPROP)+PNUC(IPROP)*OVLP(I,I),I=ISTA,IEND)
          WRITE(6,*)
        END DO

        END IF
       END IF
      END DO
      Call CollapseOutput(0,'Expectation values for input states')
      END IF

      IFON=1
      X=0.0D0
      DO I=2,NSTATE
       DO J=1,I-1
        X=MAX(X,ABS(OVLP(I,J)))
       END DO
      END DO
      IF (X.GE.1.0D-6) IFON=0
      IFHD=1
      X=0.0D0
      DO I=2,NSTATE
       DO J=1,I-1
        X=MAX(X,ABS(HAM(I,J)))
       END DO
      END DO
      IF (X.GE.1.0D-6) IFHD=0
      IF(IFON.eq.0) IFHD=0

      IF(IPGLOB.GE.USUAL) THEN
       IF(IFHAM) THEN
         WRITE(6,*)
         WRITE(6,*)' HAMILTONIAN MATRIX FOR THE ORIGINAL STATES:'
         WRITE(6,*)
         IF(IFHD.eq.1) THEN
          WRITE(6,*)' Diagonal, with energies'
          WRITE(6,'(5(1X,F15.8))')(HAM(J,J),J=1,NSTATE)
         !do J=1,NSTATE
         !WRITE(6,*) HAM(J,J)
         !enddo
         ELSE
          DO ISTA=1,NSTATE,5
            IEND=MIN(ISTA+4,NSTATE)
            WRITE(6,'(10X,5(8X,A3,I4,A3))')
     &         (' | ', I, ' > ',I=ISTA,IEND)
            DO J=1,NSTATE
              WRITE(6,'(A3,I4,A3,5(2X,F16.8))')
     &           ' < ',J,' | ', (HAM(I,J),I=ISTA,IEND)
            END DO
          ENDDO
         END IF
       END IF
      END IF

      IF(IPGLOB.GE.USUAL) THEN
        WRITE(6,*)
        WRITE(6,*)'     OVERLAP MATRIX FOR THE ORIGINAL STATES:'
        WRITE(6,*)
        IF(IFON.eq.1) THEN
         WRITE(6,*)' Diagonal, with elements'
         WRITE(6,'(5(1X,F15.8))')(OVLP(J,J),J=1,NSTATE)
        ELSE
         DO ISTATE=1,NSTATE
          WRITE(6,'(5(1X,F15.8))')(OVLP(ISTATE,J),J=1,ISTATE)
         END DO
        END IF
      END IF

* Addition by A.Ohrn. If ToFile keyword has been specified, we put
* numbers on the auxiliary rassi-to-qmstat file.
      If(ToFile) then
        If(.not.IfHam) then
          Write(6,*)
          Write(6,*)'You ask me to print hamiltonian, but there is '
     &//'none to print!'
          Call Abend()
        Endif
        Call DaName(LuEig,FnEig)
        iDisk=0
        Do iState=1,nState
          Do jState=1,iState
            Call dDaFile(LuEig,1,Ham(iState,jState),1,iDisk)
          Enddo
        Enddo
        Do iState=1,nState
          Do jState=1,iState
            Call dDaFile(LuEig,1,OvLp(iState,jState),1,iDisk)
          Enddo
        Enddo
*-- File is closed in eigctl.
        Call DaClos(LuEig)
      Endif
* End of addition by A.Ohrn.

      IF(IFHAM .AND. (IFHDIA.OR.IFSHFT)) THEN
        DO ISTATE=1,NSTATE
          IF(.NOT.IFSHFT) ESHFT(ISTATE)=0.0D0
          IF(IFHDIA) ESHFT(ISTATE)=ESHFT(ISTATE)+
     &              (Work(LHDIAG+ISTATE-1)-HAM(ISTATE,ISTATE))
        END DO
        DO ISTATE=1,NSTATE
         DO JSTATE=1,NSTATE
          HAM(ISTATE,JSTATE)=HAM(ISTATE,JSTATE)+
     &      0.5D0*(ESHFT(ISTATE)+ESHFT(JSTATE))*OVLP(ISTATE,JSTATE)
         END DO
        END DO
        IFHD=1
        DO I=2,NSTATE
         DO J=1,J-1
          IF(ABS(HAM(I,J)).GE.1.0D-10) IFHD=0
         END DO
        END DO
        IF(IFON.eq.0) IFHD=0
        IF(IPGLOB.GE.USUAL) THEN
         WRITE(6,*)
         WRITE(6,*)' USER-MODIFIED HAMILTONIAN FOR THE ORIGINAL STATES:'
         WRITE(6,*)'(With user shifts, and/or replaced diagonal'
         WRITE(6,*)' elements, including overlap corrections.)'
         IF(IFHD.eq.1) THEN
          WRITE(6,*)' Diagonal, with energies'
          WRITE(6,'(5(1X,F15.8))')(HAM(J,J),J=1,NSTATE)
         !do J=1,NSTATE
         !WRITE(6,*) HAM(J,J)
         !enddo
         ELSE
          DO ISTATE=1,NSTATE
            WRITE(6,'(5(1X,F15.8))')(HAM(ISTATE,J),J=1,ISTATE)
          END DO
         END IF
        END IF
      END IF
CPAM00 End of updated HDIA/SHIFT section.

      IF(IPGLOB.GT.SILENT .and. PRMER) THEN
      WRITE(6,*)
      Call CollapseOutput(1,'Matrix elements for input states')
      WRITE(6,'(3X,A)')     '--------------------------------'
      WRITE(6,*)
      WRITE(6,*)' MATRIX ELEMENTS OF 1-ELECTRON OPERATORS'
      WRITE(6,*)' FOR THE RASSCF INPUT WAVE FUNCTIONS:'
      WRITE(6,*)
      DO IPROP=1,NPROP
       IF(IPUSED(IPROP).NE.0) THEN
        DO ISTA=1,NSTATE,NCOL
          IEND=MIN(NSTATE,ISTA+NCOL-1)
          WRITE(6,*)
          WRITE(6,'(1X,A,A8,A,I4)')
     &  'PROPERTY: ',PNAME(IPROP),'   COMPONENT:*',ICOMP(IPROP)
          WRITE(6,'(1X,A,3D17.8)')
     &'ORIGIN    :',(PORIG(I,IPROP),I=1,3)
          WRITE(6,'(1X,A,I8,3I17)')
     &'STATE     :',(I,I=ISTA,IEND)
          WRITE(6,*)
          DO J=1,NSTATE
            WRITE(6,'(1X,I5,6X,4(1x,G16.9))')
     & J,(PROP(J,I,IPROP)+PNUC(IPROP)*OVLP(J,I),I=ISTA,IEND)
          END DO
          WRITE(6,*)
        END DO
       END IF
      END DO
      Call CollapseOutput(0,'Matrix elements for input states')
      END IF
cnf
      If (IfDCpl) Then
         Call Get_iScalar('Unique atoms',natom)
         Call Allocate_Work(iNucChg,natom)
         Call Get_dArray('Effective nuclear Charge',Work(iNucChg),nAtom)
         nST = nState*(nState+1)/2
         Call Allocate_Work(iDerCpl,3*natom*nST)
         Call AppDerCpl(natom,nST,Work(iNucChg),Prop,
     &                  Work(iDerCpl),HAM)
         Call Free_Work(iDerCpl)
         Call Free_Work(iNucChg)
      End If
cnf

      WRITE(6,*)
      CALL QEXIT(ROUTINE)
      RETURN
      END
