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
      SUBROUTINE PROPER (PROP,ISTATE,JSTATE,TDMZZ,WDMZZ)
      use rassi_global_arrays, only : JBNUM
      use RASSI_AUX
      use stdalloc, only: mma_allocate, mma_deallocate
      use Cntrl, only: NSTATE, NPROP, lSym1, lSym2, ToFile,
     &                 IRREP, PNAME, PTYPE
      use cntrl, only: FnTOM, LuTOM
      use Symmetry_Info, only: Mul, nSym=>nIrrep
      use rassi_data, only: NTDMZZ,NBST,NBASF

      IMPLICIT None
      REAL*8 PROP(NSTATE,NSTATE,NPROP)
      Real*8 TDMZZ(NTDMZZ),WDMZZ(NTDMZZ)
      Integer ISTATE, JSTATE
      Integer IOFF(8)
      CHARACTER(LEN=8) LABEL
      Save iDiskSav !(For ToFile)
      Integer, SAVE:: ICALL=0
      Real*8, Allocatable:: SCR(:,:)
      Real*8, Allocatable:: IP(:)
      Integer JOB1, JOB2, iSy12, Mask, NIP, nScr, iDisk, iDIskSav, I, J,
     &        IndCall, iProp, iType

C COMBINED SYMMETRY OF STATES:
      JOB1=JBNUM(ISTATE)
      JOB2=JBNUM(JSTATE)
      LSYM1=IRREP(JOB1)
      LSYM2=IRREP(JOB2)
      ISY12=MUL(LSYM1,LSYM2)
C THE SYMMETRY CHECK MASK:
      MASK=2**(ISY12-1)
C ALLOCATE A BUFFER FOR READING ONE-ELECTRON INTEGRALS
      NIP=4+(NBST*(NBST+1))/2
      CALL mma_allocate(IP,NIP,Label='IP')
C FIRST SET UP AN OFFSET TABLE FOR SYMMETRY BLOCKS OF TDMSCR
      Call mk_IOFF(IOFF,nSYM,NBASF,ISY12)
C CALCULATE THE SYMMETRIC AND ANTISYMMETRIC FOLDED TRANS D MATRICES
C AND SIMILAR WE-REDUCED SPIN DENSITY MATRICES
      NSCR=(NBST*(NBST+1))/2
      Call mma_allocate(SCR,nSCR,4,LABEL='SCR')
      SCR(:,:)=0.0D0
      Call MK_TWDM(nSym,TDMZZ,WDMZZ,nTDMZZ,SCR,nSCR,iOFF,NBASF,ISY12)
*
C AT THIS POINT, THE SYMMETRICALLY AND ANTISYMMETRICALLY FOLDED
C DENSITY MATRICES, AND WE-REDUCED SPIN DENSITY MATRICES, HAVE BEEN
C CALCULATED BEGINNING IN SCR.
C LOOP OVER ALL REQUIRED ONE-ELECTRON OPERATORS:
C
C-------------------------------------------
*If requested by user, put SCR in an unformatted file for later
*use by another program. (A.Ohrn)
      If(ToFile) then
        Call DaName(LuToM,FnToM)
        If(iCall.eq.0) then  !Make room for table-of-contents
          iDisk=0
          Call ICOPY(nState*(nState+1)/2,[-1],0,TocM,1)
          Call iDaFile(LuToM,1,TocM,nState*(nstate+1)/2,iDisk)
          TocM(1)=iDisk
          iDiskSav=iDisk
          iCall=1
        Endif
        If (iState.lt.jState) Then
*
*          For the rest of the code to work this cannot be violated.
*
           Write (6,*) 'Proper: iState.lt.jState'
           Call Abend()
        End If
        i=iState
        j=jState
        indCall=i*(i-1)/2+j  !Which call this is
        ToCM(indCall)=iDiskSav
        iDisk=iDiskSav
*       Write (*,*) 'IndCall,iDisk=',IndCall,iDisk
        Call dDaFile(LuToM,1,SCR,4*nSCR,iDisk) !The THING.
        iDiskSav=iDisk  !Save diskaddress.
        iDisk=0
        Call iDaFile(LuToM,1,TocM,nState*(nState+1)/2,iDisk)
                            !Put table of contents.
        Call DaClos(LuToM)
      Endif
*End of ToFile
*     Write (*,*) 'ISTATE,JSTATE=',ISTATE,JSTATE
      DO IPROP=1,NPROP
        PROP(ISTATE,JSTATE,IPROP)=0.0D00
        LABEL=PNAME(IPROP)

        CALL UPCASE(LABEL)

c If the user wants the ASD term, it is the same as
c the EF2 term without the nuclear contribution
c the new magnetic integrals are calculated from X2C
c the old spin-dependent part is thus denoted as ASDO
c O for old
        IF(LABEL(1:4).EQ.'ASDO') THEN
          LABEL(1:4) = 'EF2 '
          !write(6,*)"EF2---->ASD Here"
        END IF
        IF(LABEL(1:4).EQ.'TMOM') CYCLE

        IF(LABEL(1:4).EQ.'PSOP') THEN
          LABEL(1:4) = 'PSOI'
        !write(6,*)"PSOI---->PSOP Here"
        END IF

        IF(LABEL(1:6).EQ.'DMP   ') THEN
          LABEL(1:6) = 'DMS  1'

        END IF

        ITYPE=0
        IF(PTYPE(IPROP).EQ.'HERMSING') ITYPE=1
        IF(PTYPE(IPROP).EQ.'ANTISING') ITYPE=2
        IF(PTYPE(IPROP).EQ.'HERMTRIP') ITYPE=3
        IF(PTYPE(IPROP).EQ.'ANTITRIP') ITYPE=4
        IF(ITYPE.EQ.0) THEN
          WRITE(6,*)'RASSI/PROPER internal error.'
          WRITE(6,*)'Erroneous property type.'
          WRITE(6,*)'PTYPE(IPROP)=',PTYPE(IPROP)
          CALL ABEND()
        END IF
*
        Call MK_PROP(PROP,IPROP,ISTATE,JSTATE,LABEL,ITYPE,
     &               IP,NIP,SCR,NSCR,MASK,ISY12,IOFF)
*
      END DO
      Call mma_deallocate(SCR)
      Call mma_deallocate(IP)

      END SUBROUTINE PROPER
