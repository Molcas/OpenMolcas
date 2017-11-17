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
*  RDVEC
*
*> @brief
*>   A routine to read MO coefficients, occupation numbers, one-electron energies and type index information from ``INPORB`` file
*> @author V. Veryazov
*>
*> @details
*> New version of ::rdvec routine.
*> ::RDVEC is a wrapper to ::RDVEC_, which read UHF
*> information from ``INPORB`` file.
*>
*> \p Label defines the type of information to read from ``INPORB`` file
*> Valid targets are: ``C``---CMO, ``O``---OCC, ``E``---EORB, ``I``---INDT
*> ``A``---alpha values, ``B``---beta values
*>
*> ::RdVec checks that \p NBAS / \p NORB information is consistent,
*> and reacts according to \p iWarn. ``0``: No checks for \p NBAS / \p NORB;
*> ``1``: Print error message; ``2``: ::Abend.
*>
*> Example: Get CMO coeff. and OCC for RHF:
*>
*> \code
*> Call RdVec('INPORB',Lu,'CO',NSYM,NBAS,NBAS,CMO,OCC,Dummy,iDummy,Title,0,iErr)
*> \endcode
*>
*> @param[in]  Name  File name
*> @param[in]  LU_   Unit number
*> @param[in]  LABEL Task
*> @param[in]  NSYM  N symmetries
*> @param[in]  NBAS  N basis functions
*> @param[in]  NORB  N orbitals
*> @param[out] CMO   MO coefficients
*> @param[out] OCC   Occupations
*> @param[out] EORB  One electron energies
*> @param[out] INDT  Type Index information
*> @param[out] TITLE Title of orbitals
*> @param[in]  IWARN Warning level
*> @param[out] IERR  Return code
************************************************************************
      SUBROUTINE RDVEC(Name,LU_,LABEL,NSYM,NBAS,NORB,
     &   CMO, OCC, EORB, INDT,TITLE,iWarn,iErr)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NBAS(NSYM),NORB(NSYM),CMO(*),OCC(*), INDT(*), EORB(*)
      CHARACTER*(*) TITLE, Name, Label
      Dimension vDum(2)
      Call RdVec_(Name,LU_,LABEL,0,NSYM,NBAS,NORB,
     &   CMO,vDum, OCC, vDum, EORB, vDum,
     &   INDT,TITLE,iWarn,iErr,iWFtype)
       RETURN
       END

      SUBROUTINE RDVEC_(Name,LU_,LABEL,IUHF,NSYM,NBAS,NORB,
     &   CMO,CMO_ab, OCC, OCC_ab, EORB, EORB_ab,
     &   INDT,TITLE,iWarn,iErr,iWFtype)
* --------------------------------------------------------------------------------
*  Advanced RdVec (to remove all clones!)
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
#include "WrkSpc.fh"
      DIMENSION NBAS(NSYM),NORB(NSYM),CMO(*),OCC(*), INDT(*), EORB(*)
      DIMENSION CMO_ab(*),OCC_ab(*), EORB_ab(*)
      CHARACTER*(*) TITLE, Name, Label
      CHARACTER LINE*256,FMT*40
      LOGICAL Exist
      Character*7 Crypt, CryptUp
      Character*10 Buff
      Character*2 sDummy
*
* Note! the size of Magic must be exact (thanks to MS formatted inporb!)
*
      Character*8 Location
      data Crypt   /'fi123sd'/
      DATA CryptUP /'FIXXXSD'/
#include "inporbfmt.fh"
      Location='rdVec_'
      Line='not defined yet'
*
* Analyze Label
*
      iCMO=0
      iOCC=0
      iEne=0
      iInd=0
      iBeta=0
      If(index(Label,'C').ne.0) iCMO=1
      If(index(Label,'O').ne.0) iOcc=1
      If(index(Label,'E').ne.0) iEne=1
      If(index(Label,'I').ne.0) iInd=1
      If(index(Label,'A').ne.0) iBeta=-1
      If(index(Label,'B').ne.0) iBeta=1
*----------------------------------------------------------------------*
* Open file Name                                                       *
*----------------------------------------------------------------------*
      iErr=0
      Lu=Lu_
      Call OpnFl(Name,Lu,Exist)
      If (.Not.Exist) Then
        Write(6,*) 'RdVec: File ',Name(1:index(Name,' ')),' not found!'
        Call Abend()
      End If
      REWIND (LU)
*----------------------------------------------------------------------*
* Check version!                                                       *
*----------------------------------------------------------------------*
      READ(LU,'(A256)',END=999,ERR=999) Line
      iVer=0
      Do jVer=1,mxVer
        if(Magic(jVer).eq.Line(1:len(Magic(jVer)))) iVer=jVer
      End Do

      If(iVer.eq.0) Then
        Call SysAbendMsg(Location,'INPORB file in old format',' ')
      End If
*----------------------------------------------------------------------*
* INFO section, read it unconditionally                                *
*----------------------------------------------------------------------*
50    READ(LU,'(A256)',END=999,ERR=999) Line
      If(Line(1:5).ne.'#INFO') goto 50
      Read(Lu,'(a)',end=999,err=999) Title
      Read(Lu,'(a)',End=999,Err=999) Line
      Line(76:80)='0 0 0'
      Read(Line,*) myiUHF,myNSYM,iWFtype
*     Read(Lu,*,end=999,err=999) myiUHF,myNSYM
* In case of UHF mismatch:
      If(myiUHF.ne.iUHF) Then
* Stop if UHF requested, but the INPORB is not UHF
        If(myiUHF.eq.0)
     &    Call SysAbendFileMsg(Location,Name,'IUHF does not match',' ')
* With a UHF INPORB, only go on if alpha or beta orbitals
* explicitly requested
        If(iUHF.eq.0.and.iBeta.eq.0)
     &    Call SysAbendFileMsg(Location,Name,'IUHF does not match',' ')
      End If
      If(myNSYM.ne.NSYM) Call SysAbendFileMsg(Location,Name,
     &   'NSYM does not match',' ')
      Call GetMem('MYNBAS','Allo','Inte',imyNBAS,NSYM)
      Call GetMem('MYNORB','Allo','Inte',imyNORB,NSYM)
      Read(Lu,*,end=999,err=999) (iWork(imyNBAS+i-1),i=1,NSYM)
      Read(Lu,*,end=999,err=999) (iWork(imyNORB+i-1),i=1,NSYM)
*----------------------------------------------------------------------*
* Do checks                                                            *
*----------------------------------------------------------------------*
        Do i=1,NSYM
          If(iWork(imyNBAS+i-1).ne.NBAS(i)) Then
            Line='NBAS does not match'
            If(iWarn.eq.1) Then
              Call SysWarnMsg(Location,Line,' ')
            Else
              Call SysAbendFileMsg(Location,Name,Line,' ')
            End If
          End If
        End Do
      If(iWarn.gt.0) Then
        Do i=1,NSYM
          If(iWork(imyNORB+i-1).lt.NORB(i)) Then
            Line='NORB does not match'
            If(iWarn.eq.1) Then
              Call SysWarnMsg(Location,Line,' ')
            Else
              Call SysAbendFileMsg(Location,Name,Line,' ')
            End If
          End If
        End Do
      End If
      Call GetMem('MYNBAS','Free','Inte',imyNBAS,NSYM)
*----------------------------------------------------------------------*
* ORB section                                                          *
*----------------------------------------------------------------------*
      If(iCMO.eq.1) Then
        nDiv = nDivOrb(iVer)
        FMT = FmtOrb(iVer)
        Rewind(LU)
51      READ(LU,'(A256)',END=999,ERR=999) Line
        If(Line(1:4).ne.'#ORB') goto 51
        KCMO  = 0
        Do ISYM=1,NSYM
          Do IORB=1,iWork(imyNORB+ISYM-1)
            Do IBAS=1,NBAS(ISYM),NDIV
              IBASEND=MIN(IBAS+NDIV-1,NBAS(ISYM))
111           READ(LU,'(A256)',END=999,ERR=999) LINE
              If(LINE(1:1).EQ.'*') GOTO 111
              If(iOrb.le.nOrb(iSym)) Then
                 READ(LINE,FMT,err=888,end=888)
     &               (CMO(I+KCMO),I=IBAS,IBASEND)
              End If
            End Do
            If(iOrb.le.nOrb(iSym)) KCMO=KCMO+NBAS(ISYM)
          End Do
        End Do
        If(iUHF.eq.1.or.iBeta.eq.1) Then
52        READ(LU,'(A256)',END=999,ERR=999) Line
          If(Line(1:5).ne.'#UORB') goto 52
          KCMO  = 0
          Do ISYM=1,NSYM
            Do IORB=1,iWork(imyNORB+ISYM-1)
              Do IBAS=1,NBAS(ISYM),NDIV
                IBASEND=MIN(IBAS+NDIV-1,NBAS(ISYM))
112             READ(LU,'(A256)',END=999,ERR=999) LINE
                If(LINE(1:1).EQ.'*') GOTO 112
                If(iOrb.le.nOrb(iSym)) Then
                  If (iBeta.eq.1) Then
                    READ(LINE,FMT,err=888,end=888)
     &                  (CMO(I+KCMO),I=IBAS,IBASEND)
                  Else
                    READ(LINE,FMT,err=888,end=888)
     &                  (CMO_ab(I+KCMO),I=IBAS,IBASEND)
                  End If
                End If
              End Do
              If(iOrb.le.nOrb(iSym)) KCMO=KCMO+NBAS(ISYM)
            End Do
          End Do
        End If ! iUHF
      End If ! iCMO
*----------------------------------------------------------------------*
* OCC section                                                          *
*----------------------------------------------------------------------*
      If(iOcc.eq.1) Then
        nDiv = nDivOccMR(iVer)
        FMT = FmtOccMR(iVer)
        Rewind(LU)
53      READ(LU,'(A256)',END=999,ERR=999) Line
        If(Line(1:4).ne.'#OCC') goto 53
        KOCC  = 0
        Do ISYM=1,NSYM
          Do IORB=1,iWork(imyNORB+ISYM-1),NDIV
            IORBEND=MIN(IORB+NDIV-1,iWork(imyNORB+ISYM-1))
113         READ(LU,'(A256)',END=999,ERR=999) LINE
            If(LINE(1:1).EQ.'*') GOTO 113
            READ(LINE,FMT,err=888,end=888)
     &              (OCC(I+KOCC),I=IORB,IORBEND)
          End Do
*         KOCC=KOCC+iWork(imyNORB+ISYM-1)
          KOCC=KOCC+nOrb(iSym)
        End Do
        If(iUHF.eq.1.or.iBeta.eq.1) Then
54        READ(LU,'(A256)',END=999,ERR=999) Line
          If(Line(1:5).ne.'#UOCC') goto 54
          KOCC=0
          Do ISYM=1,NSYM
            Do IORB=1,iWork(imyNORB+ISYM-1),NDIV
              IORBEND=MIN(IORB+NDIV-1,iWork(imyNORB+ISYM-1))
114           READ(LU,'(A256)',END=999,ERR=999) LINE
              If(LINE(1:1).EQ.'*') GOTO 114
              If (iBeta.eq.1) Then
                READ(LINE,FMT,err=888,end=888)
     &              (OCC(I+KOCC),I=IORB,IORBEND)
              Else
                READ(LINE,FMT,err=888,end=888)
     &              (OCC_ab(I+KOCC),I=IORB,IORBEND)
              End If
            End Do
*           KOCC=KOCC+iWork(imyNORB+ISYM-1)
            KOCC=KOCC+nOrb(iSym)
          End Do
        End If ! iUHF
      End If ! iOCC
*----------------------------------------------------------------------*
* ONE section                                                          *
*----------------------------------------------------------------------*
      If(iEne.eq.1) Then
        nDiv = nDivEne(iVer)
        FMT = FmtEne(iVer)
        Rewind(LU)
55      READ(LU,'(A256)',END=666,ERR=666) Line
        If(Line(1:4).ne.'#ONE') goto 55
        KOCC  = 0
        Do ISYM=1,NSYM
          Do IORB=1,iWork(imyNORB+ISYM-1),NDIV
            IORBEND=MIN(IORB+NDIV-1,iWork(imyNORB+ISYM-1))
115         READ(LU,'(A256)',END=999,ERR=999) LINE
            If(LINE(1:1).EQ.'*') GOTO 115
            READ(LINE,FMT,err=888,end=888)
     &          (EORB(I+KOCC),I=IORB,IORBEND)
          End Do
*         KOCC=KOCC+iWork(imyNORB+ISYM-1)
          KOCC=KOCC+nOrb(iSym)
        End Do
        If(iUHF.eq.1.or.iBeta.eq.1) Then
56        READ(LU,'(A256)',END=999,ERR=999) Line
          If(Line(1:5).ne.'#UONE') goto 56
          KOCC=0
          Do ISYM=1,NSYM
            Do IORB=1,iWork(imyNORB+ISYM-1),NDIV
              IORBEND=MIN(IORB+NDIV-1,iWork(imyNORB+ISYM-1))
116           READ(LU,'(A256)',END=999,ERR=999) LINE
              If(LINE(1:1).EQ.'*') GOTO 116
              If (iBeta.eq.1) Then
                READ(LINE,FMT,err=888,end=888)
     &              (EORB(I+KOCC),I=IORB,IORBEND)
              Else
                READ(LINE,FMT,err=888,end=888)
     &              (EORB_ab(I+KOCC),I=IORB,IORBEND)
              End If
            End Do
*           KOCC=KOCC+iWork(imyNORB+ISYM-1)
            KOCC=KOCC+nOrb(iSym)
          End Do
        End If ! iUHF
      End If ! iOne
*----------------------------------------------------------------------*
* INDEX section                                                        *
*----------------------------------------------------------------------*
      If(iInd.eq.1) Then
        Rewind(LU)
57      READ(LU,'(A256)',END=666,ERR=666) Line
        If(Line(1:6).ne.'#INDEX') goto 57
        if(iVer.eq.iVer10) then
        FMT='(A4)'
        NDIV  = 4
        iShift=1
        Do ISYM=1,NSYM
c         iShift=(ISYM-1)*7
          Do IORB=1,iWork(imyNORB+ISYM-1),NDIV
            READ(LU,FMT,err=666,end=666) Buff
            Do i=1,4
              IND=index(Crypt,Buff(i:i))+index(CryptUp,Buff(i:i))
              If(Buff(i:i).ne.' ') Then
                If(IND.eq.0) Then
                  WRITE(6,*) '* ERROR IN RDVEC WHILE READING TypeIndex'
                  WRITE(6,'(3A)') '* Type=',Buff(i:i), ' is UNKNOWN'
                  WRITE(6,*) '* TypeIndex information is IGNORED'
                  iErr=1
                  Close(Lu)
                  goto 777
                End If
                IndT(iShift)=IND
                iShift=iShift+1
              End If
            End Do
          End Do
        End Do
        endif ! Ver10
c Ver 11
        if(iVer.ge.iVer11) then
        FMT='(A4)'
        NDIV  = 10
        iShift=1
        Do ISYM=1,NSYM
c         iShift=(ISYM-1)*7
        read(LU,*)
          Do IORB=1,iWork(imyNORB+ISYM-1),NDIV
            READ(LU,'(a2,A10)',err=666,end=666) sDummy,Buff
            Do i=1,10
              IND=index(Crypt,Buff(i:i))+index(CryptUp,Buff(i:i))
              If(Buff(i:i).ne.' ') Then
                If(IND.eq.0) Then
                  WRITE(6,*) '* ERROR IN RDVEC WHILE READING TypeIndex'
                  WRITE(6,'(3A)') '* Type=',Buff(i:i), ' is UNKNOWN'
                  WRITE(6,*) '* TypeIndex information is IGNORED'
                  iErr=1
                  Close(Lu)
                  goto 777
                End If
                IndT(iShift)=IND
                iShift=iShift+1
              End If
            End Do
          End Do
        End Do
        endif ! Ver10


        iA=1
        iB=1
        Do iSym=1,nSym
           Call iCopy(nOrb(iSym),IndT(iA),1,IndT(iB),1)
           iA=iA+nBas(iSym)
           iB=iB+nOrb(iSym)
        End Do
      End If  ! Index
      Close(Lu)
      Goto 777
*----------------------------------------------------------------------*
* a special case - INDEX information is not found                      *
*----------------------------------------------------------------------*
666   iErr=1
      WRITE(6,*) '* TypeIndex information is IGNORED *'
      Close(Lu)
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
777   Call GetMem('MYNORB','Free','Inte',imyNORB,NSYM)
      Return
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
999   Call SysAbendFileMsg(Location,Name,
     &   'Error during reading INPORB\n',Line)
888   Call SysAbendFileMsg(Location,Name,
     &   'Error during reading INPORB\n',Line)
      End
************************************************************************
*                                                                      *
************************************************************************
      SUBROUTINE VECSORT(NSYM,NBAS,NORB,CMO,OCC,INDT,NNWORD,NEWORD,iErr)
*
* Sorting routine: sort CMO, OCC according to INDT
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "WrkSpc.fh"

      DIMENSION NBAS(NSYM),NORB(NSYM),CMO(*),OCC(*), INDT(*)
* PAM 2012: Have VecSort return an array NEWORD with new orbital indices.
* If typedef info in an orbital file is used to change the order
* of orbitals, this is done by a call to VecSort. If user, as a result,
* needs to alter orbital indices ( e.g. in the supersymmetry array
* IXSYM) he needs to know how orbitals have changed order.
* The mapping is (New orbital index)=NEWORD(Old orbital index).
      DIMENSION NEWORD(*)
* PAM 2012: End of update

      MBAS=NBAS(1)
      Do i=2,nSym
         MBAS=Max(MBAS,NBAS(i))
      End Do
      Call GetMem('TCMO','Allo','Real',iTCMO,MBAS)

      kcmo=0
      kocc=0
      iii=0

* PAM 2012: NewOrd update
* If NNWORD is .gt. 0, this indicates the caller wish to get back
* a reindicing array. Then this must be large enough:
      If (NNWORD.gt.0) Then
       nw=0
       Do ISYM=1,NSYM
        Do i=1,NORB(ISYM)
         nw=nw+1
        End Do
       End Do
       If (nw.gt.NNWORD) Then
       End If

       Do nw = 1,NNWORD
        NEWORD(nw)=nw
       End Do
      End If
* PAM 2012: End of update

      Do ISYM=1,NSYM
*---- Check Do we need make sort?
      NeedSort=0
c      print *,'indt'
c      print *,(indt(i+iii),i=1,norb(isym))
        Do I=1,NORB(ISYM)
          If(IndT(i+iii).eq.0) Then
           iErr=1
          End If
          If(i.gt.1) Then
           If(IndT(i+iii).lt.IndT(i-1+iii)) NeedSort=1
          End If
        End Do
c       print *,'NeedSort=',NeedSort
       If(NeedSort.eq.1) Then
*---- Start sort
*---- we Do have only a few types of orbitals, so sorting is a simple...
        Do iType=1,7
          ip=0
          isfirst=0
          Do i=1,NORB(ISYM)
            If(isfirst.eq.0) Then
              If(IndT(i+iii).gt.iType) Then
               isfirst=1
              End If
              If(IndT(i+iii).le.iType) Then
               ip=i
              End If
            End If
            If(isfirst.eq.1.and.IndT(i+iii).eq.iType) Then
*---- We need to shift CMO, Occ
              m=IndT(i+iii)
              q=Occ(i+KOCC)
* PAM 2012: NewOrd update
              nw=NEWORD(i+KOCC)
* PAM 2012: End of update
              Do ii=1,NBAS(ISYM)
               Work(iTCMO+ii-1)=CMO((i-1)*NORB(ISYM)+KCMO+ii)
              End Do
              Do j=i,ip+2,-1
               IndT(j+iii)=IndT(j-1+iii)
               Occ(j+KOCC)=Occ(j-1+KOCC)
* PAM 2012: NewOrd update
               NEWORD(j+KOCC)=NEWORD(j-1+KOCC)
* PAM 2012: End of update
               Do ii=1,NBAS(ISYM)
                 CMO((j-1)*NORB(ISYM)+KCMO+ii)=
     &                 CMO((j-2)*NORB(ISYM)+KCMO+ii)
               End Do
              End Do

              IndT(ip+1+iii)=m
              Occ(ip+1+KOCC)=q
* PAM 2012: NewOrd update
               NEWORD(ip+1+KOCC)=nw
* PAM 2012: End of update
              Do ii=1,NBAS(ISYM)
               CMO((ip)*NORB(ISYM)+KCMO+ii)=Work(iTCMO+ii-1)
              End Do

              ip=ip+1
            End If
          End Do
      End Do
c      print *,'sorted:'
c      print '(10i2)',(IndT(i+iii), i=1,NORB(ISYM))
c      print '(4E18.12)',(Occ(i+KOCC), i=1,NORB(ISYM))

c      Do ii=1,NBAS(ISYM)
c      print *
c      print '(4E18.12)',(CMO(i+KOCC+(ii-1)*NBAS(ISYM)),i=1,NORB(ISYM))
c      End Do

*---- End sort
       End If
*---- Next symmetry
        KCMO=KCMO+NBAS(ISYM)*NORB(ISYM)
        KOCC=KOCC+NORB(ISYM)
        iii=iii+NORB(ISYM)

      End Do
      Call GetMem('TCMO','Free','Real',iTCMO,MBAS)

        Return
      End
*
      Subroutine Chk_vec_UHF(Name,Lu,isUHF)
c routine returns isUHF based on information in INPORB
      Character *(*) Name
      CHARACTER LINE*80,Location *11
      Logical Exist
#include "inporbfmt.fh"
      Location='Chk_vec_UHF'
      Line='not defined yet'

      Call OpnFl(Name,Lu,Exist)
      If (.Not.Exist) Then
       Write (6,*) 'RdVec: File ',Name(1:index(Name,' ')),' not found!'
       Call Abend()
      End If
      REWIND (LU)
* Check version!
      READ(LU,'(A80)',END=999,ERR=999) Line
      iVer=0
      Do jVer=1,mxVer
        if(Magic(jVer).eq.Line(1:len(Magic(jVer)))) iVer=jVer
      End Do

      If(iVer.eq.0) Then
          Call SysWarnMsg(Location,
     & 'INPORB file in old format',' ')
        Call SysPutsEnd()
       isUHF=0
       close(Lu)
       return
      End If
50    READ(LU,'(A80)',END=999,ERR=999) Line
      If(Line(1:5).ne.'#INFO') goto 50
* Now Do the real job
        Read (Lu,'(a)',end=999,err=999) Line
        Read(Lu,*,end=999,err=999) isUHF
        close(Lu)
      return
999      Call SysAbendFileMsg(Location,Name,
     &     'Error during reading INPORB\n',Line)
      end
