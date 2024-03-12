!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!#define _DEBUGPRINT_
      SUBROUTINE MKGUGA(SGS,CIS)
!
!     PURPOSE: MAKE THE GUGA TABLES
!     NOTE:    TO RETAIN THE TABLES AVAILABLE FOR LATER PURPOSES
!              THE START ADRESSES OF OF THE ARRAYS ETC. ARE STORED IN
!              THREE USER DEFINED TYPES. Consult the struct.F90 and gugx.F90 files for the details.
!
#ifdef _DEBUGPRINT_
      use Definitions, only: LF => u6
#endif
      use stdalloc, only: mma_allocate
      use struct, only: SGStruct, CIStruct
      IMPLICIT None
      Type(SGStruct),Target:: SGS
      Type(CIStruct) CIS

!     CREATE THE SYMMETRY INDEX VECTOR
!
      CALL MKISM()
!
!     COMPUTE TOP ROW OF THE GUGA TABLE
!
      Call mknVert0()

!Note that we do not associate the arrays here since the are not allocated yet.
      Associate (nVert => SGS%nVert, nVert0=>SGS%nVert0, IFRAS=>SGS%IFRAS)
!
!     SET UP A FULL PALDUS DRT TABLE:
!     (INITIALLY NO RESTRICTIONS ARE PUT UP)
!
      NVERT=NVERT0
      IF(IFRAS.NE.0) THEN
         CALL mma_allocate(SGS%DRT0,NVERT0,5,Label='DRT0')
         CALL mma_allocate(SGS%DOWN0,[1,NVERT0],[0,3],Label='DOWN0')
         SGS%DRTP => SGS%DRT0
         SGS%DOWNP=> SGS%DOWN0
      ELSE
         CALL mma_allocate(SGS%DRT,NVERT,5,Label='SGS%DRT')
         CALL mma_allocate(SGS%DOWN,[1,NVERT],[0,3],Label='SGS%DOWN')
         SGS%DRTP => SGS%DRT
         SGS%DOWNP=> SGS%DOWN
      ENDIF

      CALL mkDRT0()
!
!     IF THIS IS A RAS CALCULATION PUT UP RESTRICTIONS BY DELETING
!     VERTICES WHICH VIOLATE THE FORMER.
!
      IF(IFRAS.NE.0) THEN
        CALL mkRAS()
!
!     REASSEMBLE THE DRT TABLE (REMOVE DISCONNECTED VERTICES)
!
        CALL mkDRT()

!     IF THIS IS A CAS CALCULATION PROCEED WITH THE UNRESTRICTED
!     DRT TABLE
!
      ENDIF

      SGS%DOWNP=> Null()
      SGS%DRTP => Null()
!
!     CALCULATE ARC WEIGHT.
!
      CALL MKDAW()
!
!     COMPUTE UPCHAIN TABLE AND REVERSE ARC WEIGHTS
!
      CALL MKRAW()
!
!     COMPUTE LTV TABLES.
!
      CALL MKLTV()
!
!     COMPUTE MIDLEVEL AND LIMITS ON MIDVERTICE.
!
      CALL MKMID()
!
      End Associate
!
!     EXIT
!
contains
      SUBROUTINE MKISM()
      use UnixInfo, only: ProgName
      Implicit None

      If (ProgName(1:6)=='rassi') Then

         Call mkISm_Rassi()

      Else If (ProgName(1:4)=='mclr') Then

         Call mkISm_mclr()

      Else If (ProgName(1:6)=='caspt2') Then

         Call mkISm_cp2()

      Else

         Call mkNSM()

      End If

      END SUBROUTINE MKISM

      SUBROUTINE MKISM_MCLR()
      use stdalloc, only: mma_allocate
      Implicit None

#include "Input.fh"
#include "detdim.fh"
#include "spinfo_mclr.fh"
      Integer :: iOrb, iSym, iBas

      SGS%NLEV=ntASh

      Call mma_allocate(SGS%ISM,SGS%nLev,Label='SGS%ISM')

      iOrb=0
      Do iSym=1,nSym
         Do iBas=1,nRs1(iSym)
            iOrb=iOrb+1
            SGS%ISM(iOrb)=iSym
         End Do
      End Do
      Do iSym=1,nSym
         Do iBas=1,nRs2(iSym)
            iOrb=iOrb+1
            SGS%ISM(iOrb)=iSym
         End Do
      End Do
      Do iSym=1,nSym
         Do iBas=1,nRs3(iSym)
            iOrb=iOrb+1
            SGS%ISM(iOrb)=iSym
         End Do
      End Do

      End SUBROUTINE MKISM_MCLR

      SUBROUTINE MKISM_RASSI()
      use gugx, only: LEVEL
      use stdalloc, only: mma_allocate
      Implicit None
#include "rassi.fh"
      Integer ITABS, ISYM, IT, ILEV, nSym


         nSym=SGS%nSym
         SGS%NLEV=NASHT ! Total number of active orbitals
! Allocate Level to Symmetry table ISm:
         Call mma_allocate(SGS%ISm,SGS%nLev,Label='SGS%ISm')
         ITABS=0
         DO ISYM=1,NSYM
           DO IT=1,NASH(ISYM)
             ITABS=ITABS+1
             ILEV=LEVEL(ITABS)
             SGS%ISM(ILEV)=ISYM
           END DO
         END DO

      END SUBROUTINE MKISM_RASSI

      SUBROUTINE mkism_cp2()

      use fciqmc_interface, only: DoFCIQMC
      use stdalloc, only: mma_allocate
      use gugx, only: L2ACT, LEVEL

      IMPLICIT NONE

#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"

#include "SysDef.fh"
      Integer nLev

      INTEGER IT,ITABS,ILEV,ISYM, iq

      NLEV=NASHT
      SGS%nLev = NLEV
      Call mma_allocate(SGS%ISM,NLEV,Label='ISM')
! ISM(LEV) IS SYMMETRY LABEL OF ACTIVE ORBITAL AT LEVEL LEV.
! PAM060612: With true RAS space, the orbitals must be ordered
! first by RAS type, then by symmetry.
      ITABS=0
      DO ISYM=1,NSYM
        DO IT=1,NASH(ISYM)
          ITABS=ITABS+1
! Quan: Bug in LEVEL(ITABS) and L2ACT
          if (DoCumulant .or. DoFCIQMC) then
             do iq=1,NLEV
               LEVEL(iq)=iq
               L2ACT(iq)=iq
             enddo
          endif
          ILEV=LEVEL(ITABS)
          SGS%ISM(ILEV)=ISYM
        END DO
      END DO

      END SUBROUTINE mkism_cp2

      SUBROUTINE MKNSM()
!     PUPROSE: CREATE THE SYMMETRY INDEX VECTOR
!
      use gugx, only: SGS
      use stdalloc, only: mma_allocate
      IMPLICIT None
!
! to get some dimensions
#include "rasdim.fh"
! NSM form rasscf,fh
#include "rasscf.fh"
! NSYM from general.fh
#include "general.fh"
! NGAS and NGSSH from gas.fh
#include "gas.fh"
!
      Integer IGAS, ISYM, LEV, NLEV, NSTA

      NLEV=0
      DO IGAS=1,NGAS
        DO ISYM=1,NSYM
          NSTA=NLEV+1
          NLEV=NLEV+NGSSH(IGAS,ISYM)
          DO LEV=NSTA,NLEV
            NSM(LEV)=ISYM
          END DO
        END DO
      END DO

      If (SGS%nSym/=0) Then
         SGS%nLev=nLev
         Call mma_allocate(SGS%ISM,nLev,Label='SGS%ISM')
         SGS%ISM(1:nLev)=NSM(1:nLev)
      End If

      END SUBROUTINE MKNSM

      Subroutine mknVert0()
      use Definitions, only: LF => u6
      Implicit None

      Integer IAC

      Associate (iSpin=> SGS%iSpin, nActEl=>SGS%nActEl, nLev=>SGS%nLev,  &
     &           IA0=>SGS%IA0, IB0=>SGS%IB0, IC0=>SGS%IC0, nVert0=>SGS%nVert0)

      IB0=ISPIN-1
      IA0=(NACTEL-IB0)/2
      IC0=NLEV-IA0-IB0

      IF ( ((2*IA0+IB0).NE.NACTEL) .OR.                                 &
     &     (IA0.LT.0) .OR.                                              &
     &     (IB0.LT.0) .OR.                                              &
     &     (IC0.LT.0) ) then
        Write(LF,*)'mknVert0 Error: Impossible specifications.'
        Write(LF,'(1x,a,3I8)')'NACTEL,NLEV,ISPIN:',NACTEL,NLEV,ISPIN
        Write(LF,'(1x,a,3I8)')'IA0,IB0,IC0:      ',IA0,IB0,IC0
        Write(LF,*)' This is a severe internal error, or possibly'
        Write(LF,*)' indicates a strange input which should have been'
        Write(LF,*)' diagnosed earlier. Please submit a bug report.'
        Call Abend()
      End If
      IAC=MIN(IA0,IC0)
      NVERT0=((IA0+1)*(IC0+1)*(2*IB0+IAC+2))/2-(IAC*(IAC+1)*(IAC+2))/6

      End Associate

      End Subroutine mknVert0

      SUBROUTINE mkDRT0()
#ifdef _DEBUGPRINT_
      use Definitions, only: u6
#endif
      use stdalloc, only: mma_allocate, mma_deallocate
!
!     PURPOSE: CONSTRUCT THE UNRESTRICTED GUGA TABLE
!
      IMPLICIT None
!
      Integer, PARAMETER:: LTAB=1,NTAB=2,ATAB=3,BTAB=4,CTAB=5
      Integer, Parameter :: DA(0:3)=[0,0, 1,1],                                     &
                            DB(0:3)=[0,1,-1,0],                                    &
                            DC(0:3)=[1,0, 1,0]
      Integer ADDR, ADWN, AUP, BC, BDWN, BUP, CDWN, CUP, DWN, I,  &
              MXADDR, NACTEL, STEP, VDWN, LEV, VEND, VERT, VSTA, VUP,   &
              VUPS, nTmp
      Integer, Allocatable:: TMP(:)

      Associate (nVert=>SGS%nVert, DRT=>SGS%DRTP, DOWN=>SGS%DownP,                &
                 A0=>SGS%IA0, B0=>SGS%IB0, C0=>SGS%IC0, nLev=>SGS%nLev)

      NTMP=((NLEV+1)*(NLEV+2))/2
      CALL mma_allocate(TMP,NTMP,Label='TMP')
!
!     SET UP TOP ROW
!
      NACTEL=2*A0+B0
      NLEV=A0+B0+C0
      DRT(1,LTAB)=NLEV
      DRT(1,NTAB)=NACTEL
      DRT(1,ATAB)=A0
      DRT(1,BTAB)=B0
      DRT(1,CTAB)=C0
      VSTA=1
      VEND=1!
#ifdef _DEBUGPRINT_
      Write(u6,*) 'A0,B0,C0,NVERT=',A0,B0,C0,NVERT
#endif
!
!     LOOP OVER ALL LEVELS
!
      DO LEV=NLEV,1,-1
        MXADDR=((LEV+1)*(LEV+2))/2
        DO I=1,MXADDR
          TMP(I)=0
        END DO
!
!     LOOP OVER VERTICES
!
        DO VUP=VSTA,VEND
          AUP=DRT(VUP,ATAB)
          BUP=DRT(VUP,BTAB)
          CUP=DRT(VUP,CTAB)
!
!     LOOP OVER CASES
!     AND STORE ONLY VALID CASE NUMBERS WITH ADRESSES
!
          DO STEP=0,3
            DOWN(VUP,STEP)=0
            ADWN=AUP-DA(STEP)
            IF(ADWN.LT.0) Cycle
            BDWN=BUP-DB(STEP)
            IF(BDWN.LT.0) Cycle
            CDWN=CUP-DC(STEP)
            IF(CDWN.LT.0) Cycle
            BC=BDWN+CDWN
            ADDR=1+(BC*(BC+1))/2 + CDWN
            TMP(ADDR)=4*VUP+STEP
            DOWN(VUP,STEP)=ADDR
          END DO
        END DO
        VDWN=VEND
!
!     NOW INSERT VALID CASES INTO DRT TABLE
!
        DO ADDR=1,MXADDR
          VUPS=TMP(ADDR)
          IF(VUPS.EQ.0) Cycle
          VUP=VUPS/4
          STEP=MOD(VUPS,4)
          VDWN=VDWN+1
          DRT(VDWN,ATAB)=DRT(VUP,ATAB)-DA(STEP)
          DRT(VDWN,BTAB)=DRT(VUP,BTAB)-DB(STEP)
          DRT(VDWN,CTAB)=DRT(VUP,CTAB)-DC(STEP)
          TMP(ADDR)=VDWN
        END DO
!
!     CREATE DOWN CHAIN TABLE
!
        DO VUP=VSTA,VEND
          DO STEP=0,3
            DWN=DOWN(VUP,STEP)
            IF(DWN.NE.0) DOWN(VUP,STEP)=TMP(DWN)
          END DO
        END DO
        VSTA=VEND+1
        VEND=VDWN
      END DO
! End of loop over levels.
!
!     ADDING THE ZERO LEVEL TO DRT AND DOWNCHAIN TABLE
!
      DO  I=1,5
        DRT(VEND,I)=0
      END DO
      DO STEP=0,3
        DOWN(VEND,STEP)=0
      END DO
!
!     COMPLETE DRT TABLE BY ADDING NO. OF ORBITALS AND ELECTRONS
!     INTO THE FIRST AND SECOND COLUMN
!
      DO VERT=1,VEND
        DRT(VERT,LTAB)=DRT(VERT,ATAB)+DRT(VERT,BTAB)+DRT(VERT,CTAB)
        DRT(VERT,NTAB)=2*DRT(VERT,ATAB)+DRT(VERT,BTAB)
      END DO
#ifdef _DEBUGPRINT_
      DO VERT=1,VEND
        Write (6,*) 'DRT0(:,LTAB)=',DRT(VERT,LTAB)
        Write (6,*) 'DRT0(:,NTAB)=',DRT(VERT,NTAB)
      END DO
#endif
!
      CALL mma_deallocate(TMP)

      End Associate

      END SUBROUTINE mkDRT0

      SUBROUTINE mkRAS()
      use UnixInfo, only: ProgName
      IMPLICIT None

      If (ProgName(1:5).eq.'rassi') Then
         Call rmvert(SGS)
      Else
         Call RESTR(SGS)
      End If

      END SUBROUTINE mkRAS

      SUBROUTINE mkDRT()
!
!     PURPOSE: USING THE UNRESTRICTED DRT TABLE GENERATED BY DRT0 AND
!              THE MASKING ARRAY PRODUCED BY RESTR COPY ALL VALID
!              VERTICES FROM THE OLD TO THE NEW DRT TABLE
!

      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
!
      Integer IV, IVNEW, ITAB, IC, ID, IDNEW
!
      CALL mma_allocate(SGS%DRT,SGS%nVert,5,Label='DRT')
      CALL mma_allocate(SGS%DOWN,[1,SGS%nVert],[0,3],Label='SGS%DOWN')

      Associate (nVert0=>SGS%nVert0, nVert=>SGS%nVert, &
                 iDRT0=>SGS%DRT0, iDOWN0=>SGS%Down0, &
                 iDRT =>SGS%DRT , iDOWN =>SGS%Down , &
                 iVer=>SGS%Ver)

      DO IV=1,NVERT0
        IVNEW=IVER(IV)
        IF(IVNEW.EQ.0) Cycle
        DO ITAB=1,5
          IDRT(IVNEW,ITAB)=IDRT0(IV,ITAB)
        END DO
        DO IC=0,3
          ID=IDOWN0(IV,IC)
          IDNEW=0
          IF(ID.NE.0) IDNEW=IVER(ID)
          IDOWN(IVNEW,IC)=IDNEW
        END DO
      END DO
#ifdef _DEBUGPRINT_
      DO IV=1,nVert
        Write (6,*) 'DRT(i,:)=',iDRT(IV,:)
      END DO
#endif
      End Associate

      CALL mma_deallocate(SGS%Ver)
      CALL mma_deallocate(SGS%DRT0)
      CALL mma_deallocate(SGS%DOWN0)

      END SUBROUTINE mkDRT

      SUBROUTINE MKDAW()
!     PURPOSE: CONSTRUCT DIRECT ARC WEIGHTS TABLE
!
#ifdef _DEBUGPRINT_
      use Definitions, only: LF => u6
#endif
      use stdalloc, only: mma_allocate
      IMPLICIT None
!
      Integer IC, IV, ISUM, IDWN

      CALL mma_allocate(SGS%DAW,[1,SGS%nVert],[0,4],Label='SGS%DAW')

      Associate (nVert => SGS%nVert, iDown => SGS%Down, iDaw => SGS%Daw)
!
!     BEGIN TO CONSTRUCT DOWN CHAIN TABLE
!
      DO IC=0,3
       IDAW(NVERT,IC)=0
      END DO
      IDAW(NVERT,4)=1
      DO IV=NVERT-1,1,-1
        ISUM=0
        DO IC=0,3
          IDAW(IV,IC)=0
          IDWN=IDOWN(IV,IC)
          IF(IDWN.EQ.0) Cycle
          IDAW(IV,IC)=ISUM
          ISUM=ISUM+IDAW(IDWN,4)
        END DO
        IDAW(IV,4)=ISUM
      END DO
!
#ifdef _DEBUGPRINT_
      Write(LF,*)
      Write(LF,*)' DIRECT ARC WEIGHTS:'
      DO IV=1,NVERT
        Write(LF,'(1X,I4,5X,5(1X,I6))') IV,(IDAW(IV,IC),IC=0,4)
      END DO
      Write(LF,*)
#endif

      End Associate
      END SUBROUTINE MKDAW

      SUBROUTINE MKRAW()
!
!     PURPOSE: CONSTRUCT UPCHAIN INDEX TABLE AND REVERSE ARC WEIGHTS
!
#ifdef _DEBUGPRINT_
      use Definitions, only: LF => u6
#endif
      use stdalloc, only: mma_allocate
      IMPLICIT None
!
      Integer IU, IC, IDWN, IV, ISUM

      CALL mma_allocate(SGS%UP,[1,SGS%nVert],[0,3],Label='SGS%UP')
      CALL mma_allocate(SGS%RAW,[1,SGS%nVert],[0,4],Label='SGS%RAW')

      Associate (nVert=>SGS%nVert, iDown=>SGS%Down, iUp=>SGS%UP, iRaw=>SGS%Raw)
!
!     BEGIN BY CONSTRUCTING THE UPCHAIN TABLE IUP:
!
      IUP(:,:)=0
      DO IU=1,NVERT-1
        DO IC=0,3
          IDWN=IDOWN(IU,IC)
          IF(IDWN.EQ.0) Cycle
          IUP(IDWN,IC)=IU
        END DO
      END DO
!
#ifdef _DEBUGPRINT_
      Write(LF,*)
      Write(LF,*)' THE UPCHAIN TABLE IN MKRAW:'
      DO IV=1,NVERT
        Write(LF,'(1X,I4,5X,4(1X,I6))') IV,(IUP(IV,IC),IC=0,3)
      END DO
      Write(LF,*)
#endif
!
!     USE UPCHAIN TABLE TO CALCULATE THE REVERSE ARC WEIGHT TABLE:
!
      IRAW(1,0:3)=0
      IRAW(1,4)=1
      DO IV=2,NVERT
        ISUM=0
        DO IC=0,3
          IRAW(IV,IC)=0
          IU=IUP(IV,IC)
          IF(IU.EQ.0) Cycle
          IRAW(IV,IC)=ISUM
          ISUM=ISUM+IRAW(IU,4)
        END DO
        IRAW(IV,4)=ISUM
      END DO
!
#ifdef _DEBUGPRINT_
      Write(LF,*)
      Write(LF,*)' THE REVERSE ARC WEIGHT TABLE IN MKRAW:'
      DO IV=1,NVERT
        Write(LF,'(1X,I4,5X,5(1X,I6))') IV,(IRAW(IV,IC),IC=0,4)
      END DO
      Write(LF,*)
#endif
      End Associate

      END SUBROUTINE MKRAW

      SUBROUTINE MKLTV()
!     PURPOSE: FIND THE MIDLEVEL
!
#ifdef _DEBUGPRINT_
      use Definitions, only: LF => u6
#endif
      use stdalloc, only: mma_allocate
      IMPLICIT None

      Integer, Parameter:: LTAB=1
      Integer IV, LEV

      CALL mma_allocate(SGS%LTV,[-1,SGS%NLEV],Label='LTV')

      Associate (nVert=>SGS%nVert, nLev=>SGS%nLev, iDRT=>SGS%DRT, LTV=>SGS%LTV)
!
!     SET UP A LEVEL-TO-VERTEX TABLE, LTV, AND IDENTIFY MIDVERTICES:
!
      LTV(:)=0
!
      DO IV=1,NVERT
        LEV=IDRT(IV,LTAB)
        LTV(LEV)=LTV(LEV)+1
      End Do
!
      DO LEV=NLEV,0,-1
        LTV(LEV-1)=LTV(LEV-1)+LTV(LEV)
      End Do
!
      DO LEV=-1,NLEV-1
        LTV(LEV)=1+LTV(LEV+1)
      End Do

      End Associate

      END SUBROUTINE MKLTV

      SUBROUTINE MKMID()
!     PURPOSE: FIND THE MIDLEVEL
!
#ifdef _DEBUGPRINT_
      use Definitions, only: LF => u6
#endif
      IMPLICIT None
!
      Integer IV, MINW, MV, NW, IL

      Associate (nVert=>SGS%nVert, nLev=>SGS%nLev, MidLev => SGS%MidLev, &
                 MvSta=>SGS%MVSta, MVEnd=>SGS%MVEnd,   &
                 MxUp=>SGS%MxUP, MxDWN=>SGS%MxDwn, iDAW=>SGS%Daw,        &
                 iRaw=>SGS%Raw, LTV=>SGS%LTV, nMidV=>CIS%nMidV)
!
!     USE IDAW,IRAW TABLES TO DETERMINE MIDLEV.
!     THE ASSUMPTION IS THAT A BALANCED NUMBER OF UPPER/LOWER WALKS
!     IS THE BEST CHOICE.
!
!hrl 980529 fix for nLev=0 (no orbitals in any active space)
!     Since LTV(-1:nLev)  and the statement after the loop
!     MVSta=LTV(MidLev) we have the condition MidLev>=nLev
!     Hence MidLev=1 is inappropriate for nLev=0
!     MIDLEV=1
!
      If (nLev==0) Then
         MIDLEV=0
      Else
         MIDLEV=1
      End If
      MINW=1000000
      DO IL=1,NLEV-1
        NW=0
        DO IV=LTV(IL),LTV(IL-1)-1
          NW=NW+IRAW(IV,4)-IDAW(IV,4)
        END DO
        NW=ABS(NW)
        IF(NW.GE.MINW) Cycle
        MIDLEV=IL
        MINW=NW
      END DO
      MVSta=LTV(MIDLEV)
      MVEnd=LTV(MIDLEV-1)-1
      nMidV=MVEnd-MVSta+1
!
!     NOW FIND THE MAX NUMBERS OF UPPER AND LOWER WALKS. RESPECTIVELY
!     (DISREGARDING SYMMETRY)
!
      MXUP=0
      MXDWN=0
      DO MV=MVSta,MVEnd
        if(MXUP<IRAW(MV,4)) MXUP=IRAW(MV,4)
        if(MXDWN<IDAW(MV,4)) MXDWN=IDAW(MV,4)
      END DO
!
#ifdef _DEBUGPRINT_
      Write(LF,*)
      Write(LF,'(A,I3)')' MIDLEVEL =             ',MIDLEV
      Write(LF,'(A,I3)')' NUMBER OF MIDVERTICES =',NMIDV
      Write(LF,'(A,I3)')' FIRST MIDVERTEX =      ',MVSta
      Write(LF,'(A,I3)')' LAST MIDVERTEX =       ',MVEnd
      Write(LF,'(A,I3)')' MAX. NO UPPER WALKS=   ',MXUP
      Write(LF,'(A,I3)')' MAX. NO LOWER WALKS=   ',MXDWN
      Write(LF,*)
#endif

      End Associate

      END SUBROUTINE MKMID
      END SUBROUTINE MKGUGA

      SUBROUTINE MKGUGA_FREE(SGS,CIS,EXS)
!
!     PURPOSE: FREE THE GUGA TABLES
!
      use struct, only: SGStruct, CIStruct, EXStruct
      IMPLICIT None
      Type(SGStruct),Target:: SGS
      Type(CIStruct) CIS
      Type(EXStruct) EXS
!
      Call sgclose()

      Call cxclose()

contains

Subroutine SGClose()
use stdalloc, only: mma_deallocate
Implicit None

If (Allocated(SGS%ISM)) Call mma_deallocate(SGS%ISM)
If (Allocated(SGS%DRT0)) Call mma_deallocate(SGS%DRT0)
If (Allocated(SGS%DOWN0)) Call mma_deallocate(SGS%DOWN0)
If (Allocated(SGS%DRT)) Call mma_deallocate(SGS%DRT)
If (Allocated(SGS%DOWN)) Call mma_deallocate(SGS%DOWN)
If (Allocated(SGS%UP)) Call mma_deallocate(SGS%UP)
If (Allocated(SGS%MAW)) Call mma_deallocate(SGS%MAW)
If (Allocated(SGS%LTV)) Call mma_deallocate(SGS%LTV)
If (Allocated(SGS%DAW)) Call mma_deallocate(SGS%DAW)
If (Allocated(SGS%RAW)) Call mma_deallocate(SGS%RAW)
If (Allocated(SGS%SCR)) Call mma_deallocate(SGS%SCR)
If (Allocated(SGS%Ver)) Call mma_deallocate(SGS%Ver)
SGS%DRTP => Null()
SGS%DOWNP => Null()

end Subroutine SGClose

Subroutine CXClose()
use stdalloc, only: mma_deallocate

IF (Allocated(CIS%NOW)) Call mma_deallocate(CIS%NOW)
IF (Allocated(CIS%IOW)) Call mma_deallocate(CIS%IOW)
IF (Allocated(CIS%NCSF)) Call mma_deallocate(CIS%NCSF)
IF (Allocated(CIS%NOCSF)) Call mma_deallocate(CIS%NOCSF)
IF (Allocated(CIS%IOCSF)) Call mma_deallocate(CIS%IOCSF)
IF (Allocated(CIS%ICase)) Call mma_deallocate(CIS%ICase)

IF (Allocated(EXS%NOCP)) Call mma_deallocate(EXS%NOCP)
IF (Allocated(EXS%IOCP)) Call mma_deallocate(EXS%IOCP)
IF (Allocated(EXS%ICoup)) Call mma_deallocate(EXS%ICoup)
IF (Allocated(EXS%VTab)) Call mma_deallocate(EXS%VTab)
IF (Allocated(EXS%SGTMP)) Call mma_deallocate(EXS%SGTMP)
IF (Allocated(EXS%MVL)) Call mma_deallocate(EXS%MVL)
IF (Allocated(EXS%MVR)) Call mma_deallocate(EXS%MVR)
If (Allocated(EXS%USGN)) Call mma_deallocate(EXS%USGN)
If (Allocated(EXS%LSGN)) Call mma_deallocate(EXS%LSGN)

end subroutine CXClose

      END SUBROUTINE MKGUGA_FREE
