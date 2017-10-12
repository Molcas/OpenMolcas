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
* Copyright (C) Valera Veryazov                                        *
************************************************************************
*  WRVEC
*
*> @brief
*>   A routine to write MO coefficients, occupation numbers, one-electron energies and type index information to ``INPORB`` file
*> @author V. Veryazov
*>
*> @details
*> New version of ::wrvec routine.
*> ::WRVEC is a wrapper to ::WRVEC_, which writes UHF
*> information to ``INPORB`` file.
*>
*> \p Label defines the type of information to write to ``INPORB`` file.
*> Valid targets are: ``C``---CMO, ``O``---OCC, ``E``---EORB, ``I``---INDT, ``A``---Append Index
*>
*> Example: Write CMO coeff. for RHF:
*>
*> \code
*> Call WrVec('INPORB',Lu,'C',NSYM,NBAS,NBAS,CMO,Dummy,Dummy,iDummy,Title)
*> \endcode
*>
*> @param[in] Name  File name
*> @param[in] LU_   Unit number
*> @param[in] LABEL Task
*> @param[in] NSYM  N symmetries
*> @param[in] NBAS  N basis functions
*> @param[in] NORB  N orbitals
*> @param[in] CMO   MO coefficients
*> @param[in] OCC   Occupations
*> @param[in] EORB  One electron energies
*> @param[in] INDT  Type Index information
*> @param[in] TITLE Title of orbitals
************************************************************************
      SUBROUTINE WRVEC(Name,LU_,LABEL,NSYM,NBAS,NORB,CMO,
     & OCC, EORB, INDT,TITLE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NBAS(NSYM),NORB(NSYM),CMO(*),OCC(*),EORB(*),INDT(7,8)
      CHARACTER*(*) TITLE, Name,LABEL
      Dimension vDum(2)
      CALL WrVec_(Name,LU_,LABEL,0,NSYM,NBAS,NORB,CMO,
     & vDum, OCC, vDum,
     & EORB,vDum, INDT,TITLE,0)
       RETURN
       END

      SUBROUTINE WRVEC_(Name,LU_,LABEL,IUHF,NSYM,NBAS,NORB,CMO,
     & CMO_ab, OCC, OCC_ab,
     & EORB,EORB_ab, INDT,TITLE,iWFtype)
*
* The routine to dump information to INPORB
*
* --------------------------------------------------------------------------------
* iWFtype =  0  -- Unknown origin of orbitals
*            1  -- Orbitals for Guessorb
*            2  -- Orbitals for closed shell HF
*            3  -- Orbitals for closed shell DFT
*            4  -- Orbitals for unrestricted HF
*            5  -- Orbitals for unrestricted DFT
*            6  -- Natural orbitals for unrestricted HF
*            7  -- Natural orbitals for unrestricted DFT
*            8  --
* --------------------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NBAS(NSYM),NORB(NSYM),CMO(*),OCC(*),EORB(*),INDT(*)
      DIMENSION CMO_ab(*),OCC_ab(*), EORB_ab(*)
      CHARACTER*(*) TITLE, Name,LABEL
      CHARACTER*8 Line
      CHARACTER FMT*40
      Logical Exist
      Character*10 Buff
      Dimension iBuff(0:7)
      Character*7 Crypt
C-SVC: variable to hold birth certificate
      Character cDNA*256
      Logical IsBorn
      Data Crypt/'fi123sd'/
#include "inporbfmt.fh"
*
*

* Analyze Label
      iCMO=0
      iOCC=0
      iEne=0
      iInd=0
      iTwoE=0
      iAppend=0
      iExtras=0
      if(index(Label,'C').ne.0) iCMO=1
      if(index(Label,'O').ne.0) iOcc=1
      if(index(Label,'E').ne.0) iEne=1
      if(index(Label,'I').ne.0) iInd=1
      if(index(Label,'T').ne.0) iTwoE=1
      if(index(Label,'A').ne.0) iAppend=1
      if(index(Label,'X').ne.0) iExtras=1
      isIndex=0
      Lu=Lu_
      Call OpnFl(Name,Lu,Exist)
      Rewind(Lu)
      if(iAppend.eq.1) then
       iInd=1
50     read(Lu,'(A)',end=102,err=102) Line
       if(Line.eq.'#INDEX') then
         isIndex=1
         goto 100
        endif
       goto 50
102    continue
       Call Append_file(LU)
c#ifdef NAGFOR
c        BACKSPACE(LU)
c#endif
       goto 100
      endif


* Use version 2.1!
      iVer=iVer21
*
*  Write INFO header
*
      WRITE(LU,'(A)') Magic(iVer)
      Write(Lu,'(A)') '#INFO'

      KCMO  = 0
      IF(TITLE(1:1).NE.'*') TITLE='*'//TITLE(:LEN(TITLE)-1)
      WRITE(LU,'(A)') TITLE(:mylen(TITLE))
      Write(LU,'(3i8)') IUHF, NSYM, iWFtype
      WRITE(LU,'(8i8)') (NBAS(I),I=1,NSYM)
      WRITE(LU,'(8i8)') (NORB(I),I=1,NSYM)
      Call qpg_cArray('BirthCertificate',IsBorn,nDNA)
      IF (.NOT.IsBorn) THEN
        WRITE(6,*) 'RunFile has no Birth Certificate'
      ELSE
        cDNA=' '
        Call Get_cArray('BirthCertificate',cDNA(:nDNA),nDNA)
        WRITE(LU,'(A)') '*BC:'//cDNA(:mylen(cDNA))
      ENDIF

* Extras section
      If (iTwoE.eq.1 ) iExtras=1  ! so far only case

      If(iExtras.eq.1) then
        Write(Lu,'(A)') '#EXTRAS'
        If(iTwoE.eq.1) then
          WRITE(LU,'(A)') '* ACTIVE TWO-EL ENERGY'
          WRITE(LU,'(E18.12)') EORB_ab(1)
        EndIf
      EndIf

* ORB section
      if(iCMO.eq.1) then
      NDIV=nDivOrb(iVer)
      FMT=FmtOrb(iVer)
      Write(Lu,'(A)') '#ORB'
      KCMO  = 0
      DO ISYM=1,NSYM
         DO IORB=1,NORB(ISYM)
           WRITE(LU,'(A,2I5)') '* ORBITAL',ISYM,IORB
            DO IBAS=1,NBAS(ISYM),NDIV
              IBASEND=MIN(IBAS+NDIV-1,NBAS(ISYM))
            WRITE(LU,FMT) (CMO(I+KCMO),I=IBAS,IBASEND)
            EndDo
           KCMO=KCMO+NBAS(ISYM)
         EndDo
      EndDo
      if(iUHF.eq.1) then
      Write(Lu,'(A)') '#UORB'
      KCMO  = 0
      DO ISYM=1,NSYM
         DO IORB=1,NORB(ISYM)
           WRITE(LU,'(A,2I5)') '* ORBITAL',ISYM,IORB
            DO IBAS=1,NBAS(ISYM),NDIV
              IBASEND=MIN(IBAS+NDIV-1,NBAS(ISYM))
          WRITE(LU,FMT) (CMO_ab(I+KCMO),I=IBAS,IBASEND)
            EndDo
           KCMO=KCMO+NBAS(ISYM)
         EndDo
      EndDo
      Endif  ! UHF
      Endif  ! iCMO
* OCC section

      if(iOcc.eq.1) then
      NDIV=nDivOccMR(iVer)
      FMT=FmtOccMR(iVer)
      Write(Lu,'(A)') '#OCC'
      WRITE(LU,'(A)') '* OCCUPATION NUMBERS'
      KOCC=0
      DO ISYM=1,NSYM
         DO IORB=1,NORB(ISYM),NDIV
           IORBEND=MIN(IORB+NDIV-1,NORB(ISYM))
            WRITE(LU,FMT) (OCC(I+KOCC),I=IORB,IORBEND)
         EndDo
         KOCC=KOCC+NORB(ISYM)
      EndDo

      if(iUHF.eq.1) then
      Write(Lu,'(A)') '#UOCC'
      WRITE(LU,'(A)') '* Beta OCCUPATION NUMBERS'
      KOCC=0
      DO ISYM=1,NSYM
         DO IORB=1,NORB(ISYM),NDIV
           IORBEND=MIN(IORB+NDIV-1,NORB(ISYM))
         WRITE(LU,FMT) (OCC_ab(I+KOCC),I=IORB,IORBEND)
         EndDo
         KOCC=KOCC+NORB(ISYM)
      EndDo
      Endif  ! UHF
      
      
      NDIV=nDivOcc(iVer)
      FMT=FmtOcc(iVer)
      Write(Lu,'(A)') '#OCHR'
      WRITE(LU,'(A)') '* OCCUPATION NUMBERS'
      KOCC=0
      DO ISYM=1,NSYM
         DO IORB=1,NORB(ISYM),NDIV
           IORBEND=MIN(IORB+NDIV-1,NORB(ISYM))
            WRITE(LU,FMT) (OCC(I+KOCC),I=IORB,IORBEND)
         EndDo
         KOCC=KOCC+NORB(ISYM)
      EndDo

      if(iUHF.eq.1) then
      Write(Lu,'(A)') '#UOCHR'
      WRITE(LU,'(A)') '* Beta OCCUPATION NUMBERS'
      KOCC=0
      DO ISYM=1,NSYM
         DO IORB=1,NORB(ISYM),NDIV
           IORBEND=MIN(IORB+NDIV-1,NORB(ISYM))
         WRITE(LU,FMT) (OCC_ab(I+KOCC),I=IORB,IORBEND)
         EndDo
         KOCC=KOCC+NORB(ISYM)
      EndDo
      Endif  ! UHF
      Endif  ! iOcc

* ONE section
      if(iEne.eq.1) then
      NDIV=nDivEne(iVer)
      FMT=FmtEne(iVer)
      Write(Lu,'(A)') '#ONE'
      WRITE(LU,'(A)') '* ONE ELECTRON ENERGIES'
      KOCC=0
      DO ISYM=1,NSYM
         DO IORB=1,NORB(ISYM),NDIV
           IORBEND=MIN(IORB+NDIV-1,NORB(ISYM))
            WRITE(LU,FMT) (EORB(I+KOCC),I=IORB,IORBEND)
         End Do
         KOCC=KOCC+NORB(ISYM)
      End Do

      if(iUHF.eq.1) then
      Write(Lu,'(A)') '#UONE'
      WRITE(LU,'(A)') '* Beta ONE ELECTRON ENERGIES'
      KOCC=0
      DO ISYM=1,NSYM
         DO IORB=1,NORB(ISYM),NDIV
           IORBEND=MIN(IORB+NDIV-1,NORB(ISYM))
         WRITE(LU,FMT) (EORB_ab(I+KOCC),I=IORB,IORBEND)
         End Do
         KOCC=KOCC+NORB(ISYM)
      End Do
      Endif  ! UHF
      Endif  ! iEne

* INDEX section. NOTE THIS SECTION SHOULD ALWAYS BE LAST (Gv constraint)
100   if(iInd.eq.1) then
       if(iAppend.eq.0.or.(iAppend.eq.1.and.isIndex.eq.0)) then
         Write(Lu,'(A)') '#INDEX'
       endif
      FMT='(A4)'
      iShift=0

      iBuff(0)=1
c       do i=1,7
       i=1
601    continue
       iBuff(i)=iBuff(i-1)+IndT(i+iShift)
       i=i+1
       if(i.le.7) goto 601

c       enddo

      DO ISYM=1,NSYM
      write(Lu,'(A)') '* 1234567890'
      iLab=0
      iBuff(0)=1
c       do i=1,7
       i=1
600    continue
       iBuff(i)=iBuff(i-1)+IndT(i+iShift)
       i=i+1
       if(i.le.7) goto 600

c       enddo
       Ip=1
         DO IORB=1,NORB(ISYM),10
         Buff='          '
          do i=1,10
          iBB=1
          do iB=1,7
          if(Ip.ge.iBuff(iB)) iBB=iBB+1
          enddo
          if(iBB.eq.8) then
           Buff(i:i)=' '
          else
           Buff(i:i)=Crypt(iBB:iBB)
          endif
          Ip=Ip+1
          enddo
            WRITE(LU,'(i1,A1,A10)') iLab,' ',Buff
            iLab=iLab+1
            if(iLab.gt.9) iLab=0
c            WRITE(*,FMT) Buff
         End Do
       iShift=iShift+7
      End Do

      Endif  ! iInd

      Close(Lu)
      RETURN
      END
