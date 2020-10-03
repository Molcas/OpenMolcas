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
* Copyright (C) 1991, Bjorn O. Roos                                    *
*               1991, Roland Lindh                                     *
*               Valera Veryazov                                        *
************************************************************************
      Subroutine Rdbsl(BasDir,BSLbl,Type,nCGTO,mCGTO,lAng,lCGTO,lUnit,
     &                 iAtmNr,BasisTypes,ExtBasDir)
************************************************************************
* Object: Decode the basis set label and read the basis set            *
*         from the library                                             *
*                                                                      *
* Called from: GetBS                                                   *
* Subroutines called: Decode                                           *
*                                                                      *
* Author: Bjoern Roos, Theoretical Chemistry, Chemical Centre          *
*         University of Lund, Lund, Sweden                             *
*         February 1991                                                *
* Patched: Valera Veryazov                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Logical Quit_On_Error
      common /getlnQOE/ Quit_On_Error
      Dimension nCGTO(0:lCGTO),mCGTO(0:lCGTO)
      Character*80 BSLBl,BSLB*180,string,atom,type,author,basis,blank,
     &            CGTO,CGTOm,atomb, aux,
     &            Get_Ln_Quit*180
      Character *(*) BasDir, ExtBasDir
* In case anyone would wonder: basis set path fixed to 100 char max,
* (but no problem in increasing this later on) JR 2005
      Character*256 BasLoc
#include "angtp.fh"
      Character*1 kAng(0:iTabMx),Str16*16, FileName*256, TmpString*256
      Integer StrnLn
      External StrnLn
      Logical lStop, Hit, IfTest, Exist,is_error
      Integer irecl
      Integer iLast_JR,iLast1
      Integer BasisTypes(4)
      Data IfTest/.False./
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      IfTest=.True.
#endif
      If (iTabMx.lt.lCGTO) Then
         Call WarningMessage(2,'RdBsL: iTabMx.lt.lCGTO;'//
     &                         'Update the code!')
         Call Abend()
      Else
         Do i=0,iTabMx
           kAng(i)=Angtp(i)
           Call UpCase(kAng(i))
         End Do
      End If
      If (IfTest) Then
         Write (6,*) ' Enter Rdbsl'
         Write (6,'(2a)') ' BsLbl=',BsLbl
         Write (6,'(2a)') ' BasDir=',BasDir
         Write (6,'(2a)') ' ExtDir=',ExtBasDir
      End If
      Call ResetErrorLine
      Call UpCase(BSLBl)

      iLast1=StrnLn(BasDir)
      BasLoc=BasDir
      if(ExtBasDir.ne.' ') then
      Hit=.True.
      Call Decode(BSLBl(1:80),type,2,Hit)
      endif
      call find_basis_set(BasLoc,ExtBasDir,type)
      iLast_JR=StrnLn(BasLoc)
      If (Index(BasDir,'c_Basis').ne.0) Then
         BasLoc(iLast_JR+1:iLast_JR+8) = '/c_Basis'
         iLast_JR=iLast_JR+8
      Else If (Index(BasDir,'j_Basis').ne.0) Then
         BasLoc(iLast_JR+1:iLast_JR+8) = '/j_Basis'
         iLast_JR=iLast_JR+8
      Else If (Index(BasDir,'jk_Basis').ne.0) Then
         BasLoc(iLast_JR+1:iLast_JR+9) = '/jk_Basis'
         iLast_JR=iLast_JR+9
      End If
      Call BasisTbl(BSLBl,BasLoc(1:iLast_JR))
*
      Hit=.True.
      Call Decode(BSLBl(1:80),atom,1,Hit)
      If (IfTest) write(6,'(1x,a,a)') 'Atom=',atom
      iAtmNr=Lbl2Nr(atom)
*
      Hit=.True.
      Call Decode(BSLBl(1:80),type,2,Hit)
      If (IfTest) write(6,'(1x,a,a)') 'Type=',type
*
      Hit=.True.
      Call Decode(BSLBl(1:80),author,3,Hit)
      If (IfTest) write(6,'(1x,a,a)') 'Author=',author
*
      Hit=.True.
      Call Decode(BSLBl(1:80),basis,4,Hit)
      If (IfTest) write(6,'(1x,a,a)') 'Basis=',basis
*
      Hit=.True.
      Call Decode(BSLBl(1:80),CGTO,5,Hit)
      If (IfTest) write(6,'(1x,a,a)') 'CGTO=',CGTO
*
      Hit=.False.
      Call Decode(BSLBl(1:80),Aux,6,Hit)
      If (.Not.Hit) Aux = ' '
      If (IfTest) write(6,'(1x,a,a)') 'Aux=',Aux
*
      blank=' '
      If(CGTO.eq.blank .and.
     &   (type.eq.'ECP'.or.
     &    type.eq.'PSD'.or.
     &    type.eq.'ANO')) Then
         If (IfTest) Write (6,*) ' Early exit'
         Call WarningMessage(2,
     &     'Abend in RdBsl:No CGTO basis set provided in basis label')
         Call Quit_OnUserError()
      End If
*
*     We do not need more than 10 dots
*
      j=0
      Do i=1,10
         j=j+Index(BsLbl(j+1:80),'.')
         If (j.eq.0) Exit
      End Do
      If (j.gt.0) BsLbl=BsLbl(1:j)
*
*     Open basis library
*
      TmpString=Type
      Call Upcase(TmpString)
      Call LeftAd(BasDir)
      Call LeftAd(TmpString)
      iLast1=StrnLn(BasDir)
      iLast2=StrnLn(TmpString)
      If (IfTest) Then
         Write (6,'(I3,A)') iLast1, BasDir
         Write (6,'(I3,A)') iLast2, TmpString
      End If
      Filename=BasLoc(1:iLast_JR)//'/'//TmpString(1:iLast2)
      TmpString=Author
      Call Upcase(TmpString)
      iLast4=StrnLn(Filename)
      iLast2=StrnLn(TmpString)
      TmpString=Filename(1:iLast4)//'.'//TmpString(1:iLast2)
      Call f_Inquire(TmpString,Exist)
      If (Exist) Filename=TmpString
      iLast3=StrnLn(Filename)
      If (IfTest) Write (6,'(A,A)') 'Filename=',Filename
      Call f_Inquire(Filename,Exist)
      If (.not.Exist) then
c Try to find name in trans.tbl file
        Call TransTbl(Filename)
        Call f_Inquire(Filename,Exist)
        If (.not.Exist) then
         iLast3=StrnLn(Filename)
         Call WarningMessage(2,'Basis set file '//
     &         Filename(1:iLast3)// ' does not exist!;;'//
     &         '(1) For a valence basis set: check the'//
     &         ' spelling of the basis set label'//
     &         ' and that the basis set file is present in the'//
     &         ' basis set library directory.;;'//
     &         '(2) For an external auxiliary basis set: check that '//
     &         'the basis set file is present in the appropiate '//
     &         'basis set library subdirectory.')
         Call Quit_OnUserError()
        End If
      endif
c check basistype
      Call BasisType(Filename,0,BasisTypes)
      call molcas_open_ext2(lUnit,FileName,'sequential','formatted',
     & iostat,.false.,irecl,'unknown',is_error)
c      Open(Unit=lUnit,File=Filename,
c     &        Form='FORMATTED',IOSTAT=IOStat)
      If (IOStat.ne.0) Then
         iLast3=StrnLn(Filename)
         Call WarningMessage(2,
     &         ' Problems opening basis set file '//
     &         Filename(1:iLast3))
         Call Quit_OnUserError()
       End If
       ReWind (lUnit)
*
*     loop over the basis set library to find the correct label
*
      If (IfTest) Write (6,*) ' Locate basis set label in library'
   10 BSLB = Get_Ln_Quit(lUnit,0)
      If (Quit_On_Error) Then
         iLast3=StrnLn(BsLbl)
         Call WarningMessage(2,
     &          'The requested basis set label: '//
     &          BsLbl(:iLast3)//';'//
     &          'was not found in basis library: '//Filename)
         Call Abend()
      End If
      Call UpCase(BSLB)
      If (BSLB(1:1).ne.'/') Go To 10
      n=Index(BSLB,' ')
      Do i=n,80
        BSLB(i:i)='.'
      End Do
      Hit=.True.
      Call Decode(BSLB(2:80),atomb,1,Hit)
      If (atomb.ne.atom) Go To 10
*
      If (type.ne.blank) then
         Hit=.True.
         Call Decode(BSLB(2:80),string,2,Hit)
         If (string.ne.type) Go To 10
      Endif
*
      If (author.ne.blank) then
         Hit=.True.
         Call Decode(BSLB(2:80),string,3,Hit)
         If (string.ne.author) Go To 10
      Endif
*
      If (basis.ne.blank) then
         Hit=.True.
         Call Decode(BSLB(2:80),string,4,Hit)
         If (string.ne.basis) Go To 10
      Endif
*
*     If a contraction sequence has been specified it must be identical
*     to what is in the library file if the basis set type does not
*     not allow any other contraction sequence.
*
      If (CGTO.ne.blank .and.
     &    type(1:3).ne.'ANO'.and.type.ne.'ECP'.and.type.ne.'PSD'.and.
     &    type.ne.'RYDBERG' ) Then
         Hit=.True.
         Call Decode(BSLB(2:80),string,5,Hit)
         If (string.ne.CGTO) Go To 10
      End If
*
      If (Aux.ne.blank) Then
         Hit=.True.
         Call Decode(BSLB(2:80),string,6,Hit)
         If (string.ne.aux) Go To 10
      End If
*
      If (IfTest) Write (6,*) ' Process library label'
      Hit=.True.
      Call Decode(BSLB(2:80),CGTOm,5,Hit)
      If (CGTO.eq.blank) CGTO=CGTOm
*
*     Here when the basis set label has been identified on the file
*     Now decode the CGTO label
*
      lAng=-1
      i1=1
      Do 20 i=1,80
         If(CGTO(i:i).eq.blank) go to 21
         Do k=0,lCGTO
            If (CGTO(i:i).eq.kAng(k)) then
*              Read(CGTO(i1:i-1),*) nCGTO(k)
               Str16 = CGTO(i1:i-1)
               Read(Str16,'(F16.0)') cg
               nCGTO(k)=nint(cg)
               i1=i+1
               If (k.ne.lAng+1) Then
                  Call WarningMessage(2,
     &                     'RdBsl: Error in contraction label')
                  Write (6,*) 'Conflict for ',kAng(k),' shell'
                  Write (6,*) 'Erroneous label:',CGTO
                  Call Abend()
               End If
               lAng=max(lAng,k)
               Go To 20
            End If
         End Do
   20 Continue
   21 Continue
#ifdef _DEMO_
       If (lAng.gt.1) Then
          Call WarningMessage(2,
     &             'Demo version can handle only '//
     &             's and p basis functions')
          Call Quit_OnUserError()
       End If
#endif
*     Write (*,*) ' lAng=',lAng
*
*     Check for size of contracted basis set
*
      lAngm=0
      i1=1
      Do 30 i=1,80
       If(CGTOm(i:i).eq.blank) go to 31
       Do 25 k=0,lCGTO
        If(CGTOm(i:i).eq.kAng(k)) then
*        Read(CGTOm(i1:i-1),*) mCGTO(k)
         Str16 = CGTOm(i1:i-1)
         Read(Str16,'(F16.0)') cg
         mCGTO(k)=nint(cg)
         i1=i+1
         lAngM=Max(lAngM,k)
         go to 30
        Endif
   25  Continue
   30 Continue
   31 Continue
      If (IfTest) Then
         Write (6,'(2a)') 'Type=',type
         Write (6,*) 'nCGTO=',(nCGTO(k),k=0,lCGTO)
         Write (6,*) 'mCGTO=',(mCGTO(k),k=0,lCGTO)
      End If
*     Write (*,*) ' lAngm=',lAngm
      If (lAngM.lt.lAng) Then
         Call WarningMessage(2,
     &    'Abend in RdBsl:Too high angular momentum in basis set input')
         Call Quit_OnUserError()
      End If
      lStop = .False.
      If(type.eq.'ANO') then
*       Gen.cont.: never more contracted functions than available
        Do 40 i=0,lang
           If (IfTest) Write (6,*) 'i,nCGTO(i),mCGTO(i)=',
     &                             i,nCGTO(i),mCGTO(i)
         If(nCGTO(i).gt.mCGTO(i)) then
          Write(6,4000) kAng(i),nCGTO(i),mCGTO(i)
 4000     Format(/1x,'Too many CGTOs of ',a,'-type ',I3,' Max=',I3)
          lStop = .True.
         Endif
   40   Continue
      Else If(type.eq.'ECP') then
*       Segmented basis set: never more functions than primitives
*       (to be checked later (getbs.f)) and never less functions
*       than available
        Do 41 i=0,lang
           If (IfTest) Write (6,*) 'i,nCGTO(i),mCGTO(i)=',
     &                             i,nCGTO(i),mCGTO(i)
         If(nCGTO(i).lt.mCGTO(i)) then
          Write(6,4100) kAng(i),nCGTO(i),mCGTO(i)
 4100     Format(/1x,'Too few segmented CGTOs of ',a,'-type ',I3
     &           ,' Min=',I3)
          lStop = .True.
         Endif
   41   Continue
      Endif
      If (lStop) Then
         Call WarningMessage(2,
     &     'Abend in RdBsl:Requested basis inconsistent with library')
         Call Quit_OnUserError()
      End If
      Return
      End
c
      Function Lbl2Nr(Atom)
      Integer Lbl2Nr, StrnLn
      Character*80 Atom
      Character*2 TmpLbl(2)
#include "periodic_table.fh"
*
      If (StrnLn(Atom).gt.2.or.StrnLn(Atom).le.0) Then
         Call WarningMessage(2,'The atom label;'//
     &                '-->'//Atom(1:4)//'<--;'//
     &                ' is not a proper string to define an element.')
         Call Quit_OnUserError()
      End If
*
      Lbl2Nr=-1
*
      lAtom=0
      Do i = 1, 80
         If (atom(i:i).ne.' ') Then
            lAtom=lAtom+1
         Else
            Go To 99
         End If
      End Do
 99   Continue
      If (lAtom.gt.2) lAtom=2
      TmpLbl(1)='  '
      If (lAtom.eq.1) Then
         TmpLbl(1)(2:2)=atom(1:1)
      Else If (lAtom.eq.2) Then
         TmpLbl(1)=atom(1:2)
      Else
         Go To 98
      End If
      Call UpCase(TmpLbl(1))
      Do i = 0, Num_Elem
         TmpLbl(2)=PTab(i)
        Call UpCase(TmpLbl(2))
         If (TmpLbl(2).eq.TmpLbl(1)) Then
            Lbl2Nr=i
            Go To 98
         End If
      End Do
98    Continue
      If (Lbl2Nr.eq.-1) Then
         Call WarningMessage(2,'The atom label;'//
     &               '-->'//Atom(1:4)//'<--;'//
     &               ' does not define an element.')
         Call Quit_OnUserError()
      End If
      Return
      End
c
