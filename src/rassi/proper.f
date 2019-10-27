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
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='PROPER')
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Files.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      DIMENSION TDMZZ(NTDMZZ),WDMZZ(NTDMZZ)
      DIMENSION IOFF(8)
      CHARACTER*8 LABEL
      REAL*8 PROP(NSTATE,NSTATE,NPROP)
      Save iDiskSav !(For ToFile)
      SAVE ICALL
      DATA ICALL /0/
      Real*8, Allocatable:: SCR(:,:)

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
      CALL GETMEM('IP    ','ALLO','REAL',LIP,NIP)
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
CTL2004-start
C-------------------------------------------
      IF (NONA.AND.(ISTATE.EQ.MAX(NONA_ISTATE,NONA_ISTATE))
     &        .AND.(JSTATE.EQ.MIN(NONA_ISTATE,NONA_JSTATE))) THEN
C
C IF NONADIABATIC COUPLINGS ARE REQUIRED LET'S COMPUTE THEM RIGHT NOW!
C
         CALL COMP_NAC(ISTATE, JSTATE, SCR, 4*nSCR, ISY12, IOFF, LCI1)
      END IF
C WE CONTINUE WITH THE NORMAL CALCULATION OF PROPERTIES...
C-------------------------------------------
CTL2004-end
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
*          For the rest of the code to work this can not be violated.
*
           Write (6,*) 'Proper: iState.lt.jState'
           Call Abend()
        End If
        i=iState
        j=jState
        indCall=i*(i-1)/2+j  !Which call this is
        ToCM(indCall)=iDiskSav
        ind=indCall+1
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
        IF(LABEL(1:3).EQ.'ASD') THEN
          LABEL(1:3) = 'EF2'
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
     &               WORK(LIP),NIP,SCR,NSCR,MASK,ISY12,IOFF)
*
      END DO
      Call mma_deallocate(SCR)
      CALL GETMEM('      ','FREE','REAL',LIP,NIP)
      RETURN
      END
