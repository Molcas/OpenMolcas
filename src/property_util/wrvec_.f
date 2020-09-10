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
*> Valid targets are: ``C``---CMO, ``O``---OCC, ``E``---EORB, ``I``---INDT, ``A``---Append Index, ``K``---Coordinates, ``B``---Basis section
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
      Character*12 lBuff
      Dimension iBuff(0:7)
      Character*7 Crypt
C-SVC: variable to hold birth certificate
      Character cDNA*256
      Character*120 inout
      Character*20 InpOrbVer
      Logical IsBorn
      Data Crypt/'fi123sd'/
      Integer, Save :: iVer=0
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
      iKoor=0
      iBasis=0
      if(index(Label,'C').ne.0) iCMO=1
      if(index(Label,'O').ne.0) iOcc=1
      if(index(Label,'E').ne.0) iEne=1
      if(index(Label,'I').ne.0) iInd=1
      if(index(Label,'T').ne.0) iTwoE=1
      if(index(Label,'A').ne.0) iAppend=1
      if(index(Label,'X').ne.0) iExtras=1
      if(index(Label,'K').ne.0) iKoor=1
      if(index(Label,'B').ne.0) iBasis=1
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

*
*  Get version
*
      iDefault=iVer22
      If (iVer.eq.0) Then
        Call getenvf('MOLCAS_INPORB_VERSION',InpOrbVer)
        If (InpOrbVer.eq.'') Then
          iVer=iDefault
        Else
          InpOrbVer='#INPORB '//Trim(AdjustL(InpOrbVer))
          Do jVer=1,mxVer
            if(Magic(jVer).eq.InpOrbVer) iVer=jVer
          End Do
          If (iVer.eq.0) Then
            Call WarningMessage(0,
     &           'Unknown INPORB version, using the default')
            iVer=iDefault
          End If
        End If
      End If

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
          WRITE(LU,'(E19.12)') EORB_ab(1)
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
      NDIV=nDivOcc(iVer)
      FMT=FmtOcc(iVer)
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

      NDIV=nDivOccHR(iVer)
      If (NDIV.gt.0) Then
         FMT=FmtOccHR(iVer)
         Write(Lu,'(A)') '#OCHR'
         WRITE(LU,'(A)') '* OCCUPATION NUMBERS (HUMAN-READABLE)'
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
         WRITE(LU,'(A)') '* Beta OCCUPATION NUMBERS (HUMAN-READABLE)'
         KOCC=0
         DO ISYM=1,NSYM
            DO IORB=1,NORB(ISYM),NDIV
              IORBEND=MIN(IORB+NDIV-1,NORB(ISYM))
            WRITE(LU,FMT) (OCC_ab(I+KOCC),I=IORB,IORBEND)
            EndDo
            KOCC=KOCC+NORB(ISYM)
         EndDo
         Endif  ! UHF
      End If
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
      Call getenvf('MOLCAS_SAGIT',inout)
      If (inout(1:1).eq.'y'.or.inout(1:1).eq.'Y') Then
        iKoor=iCMO
        iBasis=iCMO
      Else
        iKoor=0
        iBasis=0
      End If
      if(iKoor.eq.1.and.iBasis.eq.1) then
        in=16
        in=isfreeunit(in)
        call molcas_open(in,'ORB.std')
765     read(in,'(a)',end=766,err=766) InOut
        write(Lu,'(a)') InOut(1:mylen(InOut))
        goto 765
766     close(in)
      endif
* INDEX section. NOTE THIS SECTION SHOULD ALWAYS BE LAST (Gv constraint)
100   if(iInd.eq.1) then
       if(iAppend.eq.0.or.(iAppend.eq.1.and.isIndex.eq.0)) then
         Write(Lu,'(A)') '#INDEX'
       endif
      iShift=0
      nDiv=nDivInd(iVer)

      iBuff(0)=1
c       do i=1,7
       i=1
601    continue
       iBuff(i)=iBuff(i-1)+IndT(i+iShift)
       i=i+1
       if(i.le.7) goto 601

c       enddo

      DO ISYM=1,NSYM
      If (nSkpInd(iVer).gt.0) write(Lu,'(A)') '* 1234567890'
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
         DO IORB=1,NORB(ISYM),nDiv
         Buff='          '
          do i=1,nDiv
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
            WRITE(lBuff,FMTIND(iVer)) Buff
            If (Index(FMTIND(iVer),'X').gt.0)
     &         WRITE(lBuff(1:1),'(i1)') iLab
            WRITE(LU,'(A)') Trim(lBuff)
            iLab=iLab+1
            if(iLab.gt.9) iLab=0
         End Do
       iShift=iShift+7
      End Do

      Endif  ! iInd

      Close(Lu)
      RETURN
      END
      subroutine Koor2file(Lu)
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "real.fh"
      Dimension iOper(8)
      Dimension RotVec(3)
c      Character*16 FMT
      Character*(LENIN) AtomLbl(MxAtom)
      Character*(LENIN) Byte4
      Character*128 Line
cVV: the constant is used in all GV packages
#define R529 0.52917721067d0
        x529=R529
        write(Lu,'(A)') '#COORD'
*----------------------------------------------------------------------*
*     Read no.of symm. species                                         *
*----------------------------------------------------------------------*
      Call get_iScalar('nSym',nSym)
*----------------------------------------------------------------------*
*     Read symm. oper per symm. species                                *
*----------------------------------------------------------------------*
      Call Get_iArray('Symmetry operations',iOper,nSym)
*----------------------------------------------------------------------*
*     Read no. of unique atoms in the system                           *
*----------------------------------------------------------------------*
      Call Get_iScalar('Unique atoms',nAtoms)
*----------------------------------------------------------------------*
*     Read atom labels                                                 *
*----------------------------------------------------------------------*
      lw2=1
      Call Get_cArray('Unique Atom Names',AtomLbl,LENIN*nAtoms)
*----------------------------------------------------------------------*
*     Read coordinates of atoms                                        *
*----------------------------------------------------------------------*
      Call GETMEM('Coor','ALLO','REAL',ipCoor,3*nSym*nAtoms)
      Call Get_dArray('Unique Coordinates',Work(ipCoor),3*nAtoms)


      lw2=1
      nOper=0
      If ( nSym.eq.2 ) nOper=1
      If ( nSym.eq.4 ) nOper=2
      If ( nSym.eq.8 ) nOper=3
      nCenter=nAtoms
      Do i=1,nOper
        jOper=i+1
        If ( i.eq.3 ) jOper=5
        RotVec(1)=One
        If ( IAND(iOper(jOper),1).eq.1 ) RotVec(1)=-One
        RotVec(2)=One
        If ( IAND(iOper(jOper),2).eq.2 ) RotVec(2)=-One
        RotVec(3)=One
        If ( IAND(iOper(jOper),4).eq.4 ) RotVec(3)=-One
        newAt=0
        mCenter=nCenter
        Do iAt=0,mCenter-1
          Xold=WORK(ipCoor+iAt*3+0)
          Yold=WORK(ipCoor+iAt*3+1)
          Zold=WORK(ipCoor+iAt*3+2)
          Byte4=AtomLbl(lw2+iAt)
          Xnew=RotVec(1)*Xold
          Ynew=RotVec(2)*Yold
          Znew=RotVec(3)*Zold
          Do jAt=0,nCenter-1
             If (Byte4.eq.AtomLbl(Lw2+jAt)) Then
                Xold2=WORK(ipCoor+jAt*3+0)
                Yold2=WORK(ipCoor+jAt*3+1)
                Zold2=WORK(ipCoor+jAt*3+2)

          If ( Xnew.eq.Xold2.and.Ynew.eq.Yold2.and.Znew.eq.Zold2)
     &                   goto 999
             Endif
            Enddo
            WORK(ipCoor+nCenter*3+0)=Xnew
            WORK(ipCoor+nCenter*3+1)=Ynew
            WORK(ipCoor+nCenter*3+2)=Znew
            AtomLbl(lw2+nCenter)=Byte4
            nCenter=nCenter+1
            newAt=newAt+1
 999    Continue
        End Do
c        nCenter=nCenter+newAt
      End Do

      write (lu,*) Ncenter
      write (lu,*)
        DO IAT=0,NCENTER-1
          WRITE(LINE,'(A)') ATOMLBL(LW2+IAT)
          Byte4=ATOMLBL(LW2+IAT)(1:2)
          if(index ('0123456789',Byte4(2:2)).ne.0) Byte4(2:2)=' '
            WRITE(LINE,'(1X,A2,2X,3F15.8)') Byte4,
     &      WORK(IPCOOR+3*IAT)*x529,WORK(IPCOOR+3*IAT+1)*x529,
     &      WORK(IPCOOR+3*IAT+2)*x529
          write(lu,'(A50)') Line(1:50)
        ENDDO
      Call GETMEM('Coor','FREE','REAL',ipCoor,3*nSym*nAtoms)
      return
      end
c---------------------------------------------------
      subroutine Basi2file(Lu)
      use Basis_Info
      use Center_Info
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
      Integer ipSph(0:MxAng)
      Integer iix(2)
      Real*8 rix(2)
      Logical Found
*
      nbyte_i = iiloc(iix(2)) - iiloc(iix(1))
      nbyte_r = idloc(rix(2)) - idloc(rix(1))
*
*     Prologue
*
      write(Lu,'(A)') '#BASIS'
*
      Itemp=1+Mx_Shll*7+1+nShlls+1+nCnttp+1+64*Mx_mdc+
     +  1+ Mx_Unq +1 +8*Mx_AO+1+ Mx_Shll*3

      Icurr=0
      Call GETMEM('VVINfo','ALLO','INTE',iivv,Itemp)
      iWork(iivv+Icurr)=Mx_Shll
      Icurr=iCurr+1
c      Write(LU,'(i8)') Mx_Shll

      iWork(iivv+Icurr)=nCnttp
      Icurr=iCurr+1
      do i=1,nCnttp
        iWork(iivv+Icurr+i)=dbsc(i)%nOpt
      enddo
      Icurr=iCurr+nCnttp

c      Write(LU,'(i8)') nCnttp
c      Write(LU,'(8i8)') (nOpt(i),i=1,nCnttp)

      iWork(iivv+Icurr)=64*n_dc
      Icurr=iCurr+1
       do i=1,n_dc
          iA=iivv+iCurr
          iWork(iA+1:iA+8)=dc(i)%iCoSet(:,0)
          iCurr=iCurr+8
          iA=iivv+iCurr
          iWork(iA+1:iA+8)=dc(i)%iCoSet(:,1)
          iCurr=iCurr+8
          iA=iivv+iCurr
          iWork(iA+1:iA+8)=dc(i)%iCoSet(:,2)
          iCurr=iCurr+8
          iA=iivv+iCurr
          iWork(iA+1:iA+8)=dc(i)%iCoSet(:,3)
          iCurr=iCurr+8
          iA=iivv+iCurr
          iWork(iA+1:iA+8)=dc(i)%iCoSet(:,4)
          iCurr=iCurr+8
          iA=iivv+iCurr
          iWork(iA+1:iA+8)=dc(i)%iCoSet(:,5)
          iCurr=iCurr+8
          iA=iivv+iCurr
          iWork(iA+1:iA+8)=dc(i)%iCoSet(:,6)
          iCurr=iCurr+8
          iA=iivv+iCurr
          iWork(iA+1:iA+8)=dc(i)%iCoSet(:,7)
          iCurr=iCurr+8
       enddo

      Call ICopy(Mx_Unq,[1],0,IrrCmp,1)
      Call Get_iArray('IrrCmp',IrrCmp,Mx_Unq)

      iWork(iivv+Icurr)=Mx_Unq
      Icurr=iCurr+1
      call icopy(Mx_Unq,IrrCmp,1,iWork(iivv+Icurr),1)
c      do i=1,Mx_Unq
c        iWork(iivv+Icurr+i)=IrrCmp(i)
c      enddo
      Icurr=iCurr+Mx_Unq

c      Write(LU,'(i8)') Mx_Unq
c      call pack_me(LU,irrCmp,Mx_Unq)

*
*     Load the common LINFO
*
c      Len = iiLoc(lxEnd)-iiLoc(lxStrt)
c      Len = (Len+nbyte_i)/nbyte_i
c      Call Get_iArray('SewLInfo',lxStrt,Len)
      iWork(iivv+Icurr)=Mx_Shll
      Icurr=iCurr+1

      do i=1,Size(Shells)
        ii=0
        if(Shells(i)%Prjct) ii=1
        iWork(iivv+Icurr+i-1)=ii
      enddo
      Icurr=iCurr+Mx_Shll

      do i=1,Size(Shells)
        ii=0
        if(Shells(i)%Transf) ii=1
        iWork(iivv+Icurr+i-1)=ii
      enddo
      Icurr=iCurr+Mx_Shll

      do i=1,Size(Shells)
        ii=0
        if(Shells(i)%Aux) ii=1
        iWork(iivv+Icurr+i-1)=ii
      enddo
      Icurr=iCurr+Mx_Shll
cvv
c      do i=1,Mx_Shll
c        ii=0
c        if(Shells(i)%Frag) ii=1
c        iWork(iivv+Icurr+i)=ii
c      enddo
c      Icurr=iCurr+Mx_Shll

       write(Lu,'(i8)') iTemp

       call pack_me(Lu, iWork(iivv),iTemp)

      Call GETMEM('VVINfo','FREE','INTE',iivv,Itemp)
*
*     Load the common RINFO
*
      Len = idLoc(rxEnd)-idLoc(rxStrt)
      Len = (Len+nByte_r)/nByte_r
      Write(LU,'(i8)') Len

      Call GETMEM(' SewRInfo ','ALLO','REAL',irxStrt,Len)
      Call Get_dArray('SewRInfo',Work(irxStrt),Len)
      call pack_meR(LU,Work(irxStrt),Len)

c      Write(LU,'(4F16.8)') (Work(Info+i),i=0,Len-1)
      Call GETMEM(' SewRInfo ','FREE','REAL',irxStrt,Len)

c VVV /*
c      CLight=CLight_Info
*
*     Load the common CINFO
*
c      Len = icLoc(cxEnd)-icLoc(cxStrt)
c      Len = (Len+nByte_i)/nByte_i
c      Call Get_iArray('SewCInfo',cxStrt,Len)
c VVV */
*
       Call qpg_dArray('SewTInfo',Found,Len2)
      Write(LU,'(i8)') Len2

       Call GETMEM(' Sphere','ALLO','REAL',ipSph(0),Len2)

         Do 2 iAng = 0, iAngMx-1
            ipSph(iAng+1)= ipSph(iAng) + (iAng*(iAng+1)/2 + iAng + 1)**2
 2       Continue
      Call Get_dArray('SewTInfo',WORK(ipSph(0)),Len2)
      call pack_meR(LU,Work(ipSph(0)),Len2)

       Call GETMEM(' Sphere','FREE','REAL',ipSph(0),Len2)

      Return
      End
         subroutine pack_me(iout,iA,n)
         Dimension iA(n)
         Character*80 Str,StrC
         Character*120 StrB
#include "WrkSpc.fh"
         call getmem('Ind','ALLO','INTE',ipInd,n)
c         Allocatable:: Ind(:)
c         Allocate (Ind(n))
         iStrB=1
         ito=iA(1)
         iWork(ipInd)=1
         ii=2
         issame=1
         do i=2,n
          if(abs(iA(i)).gt.100000000) iA(i)=0
           if(ito.ne.iA(i)) then
             iWork(ipInd+ii-1)=i
             ii=ii+1
             ito=iA(i)
             issame=0
           endif
         enddo
         if(iA(n).eq.iA(n-1)) issame=1
         iWork(ipInd+ii-1)=n+issame
c Ind 1 2 3 4 5 100
          do jj=1,ii-1
           iD=iWork(ipInd+jj)-iWork(ipInd+jj-1)
          if(iD.eq.1) then
            write(Str,*) iA(iWork(ipInd+jj-1))
          else
            if(iD.eq.0) iWork(ipInd+jj)=iWork(ipInd+jj)+1
            if(iWork(ipInd+jj)-iWork(ipInd+jj-1).eq.1) then
            write(Str,*) iA(iWork(ipInd+jj-1))
            else
            write(Str,*) iWork(ipInd+jj)-iWork(ipInd+jj-1),':',
     *                   iA(iWork(ipInd+jj-1))
            endif
          endif
            icc=1
            do kk=1,80
             if(Str(kk:kk).ne.' ') then
               StrC(icc:icc)=Str(kk:kk)
               icc=icc+1
             endif
            enddo
            StrB(iStrB:iStrB+icc-1)=StrC(1:icc-1)
            iStrb=iStrB+icc
            StrB(iStrb:)=' '
            if(iStrb.ge.80) then
              write(iout,'(a)') StrB(1:iStrb)
              iStrb=1
            endif
          enddo
          if(iStrB.gt.1) write(iout,'(a)') StrB(1:iStrb)
         call getmem('Ind','FREE','INTE',ipInd,n)
          return
          end

         subroutine pack_meR(iout,A,n)
         IMPLICIT REAL*8 (A-H,O-Z)
         Dimension A(n)
         Character*80 Str,StrC
         Character*120 StrB
         Allocatable:: Ind(:)
         Allocate (Ind(n+1))
         iStrB=1
         to=A(1)
         Ind(1)=1
         ii=2
         issame=1
         do i=2,n
c          if(abs(A(i)).gt.100000000) iA(i)=0
           if(abs(to-A(i)).gt.10e-5) then
             Ind(ii)=i
             ii=ii+1
             to=A(i)
             issame=0
           endif
         enddo
         If (n.gt.1) Then
         if(abs(A(n)-A(n-1)).lt.10e-5) issame=1
         End If

         Ind(ii)=n+issame
c Ind 1 2 3 4 5 100
          do jj=1,ii-1
          if(Ind(jj+1)-Ind(jj).eq.1) then
             if(abs(A(Ind(jj))).lt.1E-6) then
              write(Str,'(A)') '0.0'
            else if (abs(A(Ind(jj))-1).lt.1E-6) then
              write(Str,'(A)') '1.0'
            else if (abs(A(Ind(jj))+1).lt.1E-6) then
              write(Str,'(A)') '-1.0'
            else
            write(Str,'(G18.4)') A(Ind(jj))
            endif
          else
            if(Ind(jj+1)-Ind(jj).eq.0) Ind(jj+1)=Ind(jj+1)+1
            if(Ind(jj+1)-Ind(jj).eq.1) then
                write(Str,'(G18.4)') A(Ind(jj))
                else
             if(abs(A(Ind(jj))).lt.1E-6) then
             write(Str,'(I8,A,A)') Ind(jj+1)-Ind(jj),':','0.0'
             else
             write(Str,'(I8,A,G18.4)') Ind(jj+1)-Ind(jj),':',A(Ind(jj))
             endif
                endif
          endif
            icc=1
            do kk=1,80
             if(Str(kk:kk).ne.' ') then
               StrC(icc:icc)=Str(kk:kk)
               icc=icc+1
             endif
            enddo
            StrB(iStrB:iStrB+icc-1)=StrC(1:icc-1)
            iStrb=iStrB+icc
            StrB(iStrb:)=' '
            if(iStrb.ge.80) then
              write(iout,'(a)') StrB(1:iStrb)
              iStrb=1
            endif
          enddo
          if(iStrB.gt.1) write(iout,'(a)') StrB(1:iStrb)
          Deallocate (Ind)
          return
          end


          subroutine unpack_me(iout,iA,n)
         Dimension iA(n)
         Character*120 Str, StrT
         ind=1
1        read(iout,'(a)') Str
c         write(*,'(a)') Str
10       if(Str.eq.' ') goto 200
20       j=index(Str,' ')
           if(j.eq.1) then
             StrT=Str(2:)
             Str=StrT
             goto 20
           endif
         StrT=Str(j+1:)
         Str(j:)=' '
         k=index(Str,':')
         if(k.eq.0) then
              read(Str,*) iA(ind)
              ind=ind+1
         else
              read(Str(k+1:),*) iA(ind)
              read(Str(1:k-1),*) l
               do ll=1,l
                 iA(ind+ll)=iA(ind)
               enddo
               ind=ind+l
         endif
         Str=StrT
         goto 10
200      continue
c         print *,'VV',ind-1,n
         if(ind-1.lt.n) goto 1
         return
         end


          subroutine unpack_meR(iout,A,n)
         IMPLICIT REAL*8 (A-H,O-Z)
         Dimension A(n)
         Character*120 Str, StrT
         ind=1
1        read(iout,'(a)') Str
c         write(*,'(a)') Str
10       if(Str.eq.' ') goto 200
20       j=index(Str,' ')
           if(j.eq.1) then
             StrT=Str(2:)
             Str=StrT
             goto 20
           endif
         StrT=Str(j+1:)
         Str(j:)=' '
         k=index(Str,':')
         if(k.eq.0) then
              read(Str,*) A(ind)
              ind=ind+1
         else
              read(Str(k+1:),*) A(ind)
              read(Str(1:k-1),*) l
               do ll=1,l
                 A(ind+ll)=A(ind)
               enddo
               ind=ind+l
         endif
         Str=StrT
         goto 10
200      continue
c         print *,'VV',ind-1,n
         if(ind-1.lt.n) goto 1
         return
         end
