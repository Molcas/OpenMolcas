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
* Copyright (C) 1990,2020,  Roland Lindh                               *
*               1990, IBM                                              *
************************************************************************
************************************************************************
*                                                                      *
*    Object: to read basis set Exponents and Contraction Coefficients  *
*            from a library file.                                      *
*            The contraction coefficients will at this point be radial *
*            normalized.                                               *
*                                                                      *
* Called from: Input                                                   *
*                                                                      *
*                                                                      *
*            DDName is the path to the library directory               *
*                                                                      *
*                                                                      *
* Calling    : RecPrt, Rdbsl                                           *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
************************************************************************
#define _ACTUAL_
#include "getbs_interface.fh"
*     Local variables
      Character(LEN=80)  MPLbl*20, Filenm, Atom, Type
      Character(LEN=256) DirName
*
      Character Basis_Lib*256, Filename*263, DefNm*13
      Integer StrnLn
      External StrnLn
      Logical UnContracted
*
      Character*180 Line, Get_Ln
      External Get_Ln
      Character*24 Words(2)                     ! CGGn
      Logical inLn1, inLn2, inLn3, Hit, IfTest,
     &        isEorb,isFock
      Integer nCGTO(0:iTabMx),mCGTO(0:iTabMx)
      Logical Found
      Real*8, Allocatable:: ExpMerged(:),Temp(:,:)
      Data DefNm/'basis_library'/
*
#include "relmp.fh"
*     IRELMP =0  .... NOPAIR (DK2)
*     IRELMP =1  .... NOPAIR (DK1)
*     IRELMP =2  .... NOPAIR (DK2)
*     IRELMP =3  .... NOPAIR (DK3)
*     IRELMP =4  .... full NOPAIR (DK3)
*     IRELMP =11 .... RESC
*     IRELMP =21 .... ZORA
*     IRELMP =22 .... ZORA-FP
*     IRELMP =23 .... IORA
*                                                                      *
************************************************************************
*                                                                      *
      Interface
#include "getecp_interface.fh"
         Subroutine RecPrt(Title,FmtIn,A,nRow,nCol)
         Character*(*) Title
         Character*(*) FmtIn
         Integer nRow,nCol
         Real*8 A(nRow,nCol)
         End Subroutine RecPrt
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
!#define _DEBUG_
#ifdef _DEBUG_
      IfTest=.True.
      iPrint=99
#else
      IfTest=.False.
      iPrint=5
#endif
      If (IfTest) iPrint=99
      ip_Dummy=-1
      dbsc(nCnttp)%FOp = .True.
      nM1=0
      nM2=0
      lAng=0
*
      If (IfTest) Write (6,'(A,A)') 'DDName=',DDName
      Line=DDName
      call UpCase(Line)
      If (Index(Line,'INLINE').ne.0.and.Index(Line,'EXPERT').eq.0) Then
         if(.not.Expert) then
           Write(6,*) 'INLINE keyword is detected'
           Write(6,*) 'EXPERT keyword is not set...'
           Write(6,*) '       We better abort right now!'
           call abend()
         endif
      EndIf
*-----Locate the first word
      iFrst = 1
      Call NxtWrd(DDName,iFrst,iEnd)
      Filenm=DDName(iFrst:iEnd)
      DirName=DDName(iFrst:iEnd)
      Call UpCase(Filenm)
      If (Index(Filenm,'INLINE').ne.0) Then
         inLn1 = .true.
         inLn2 = .true.
         inLn3 = .true.
      Else If (Index(Filenm,'MM').ne.0) Then
         dbsc(nCnttp)%IsMM = 1
         inLn1 = .True.
         inLn2 = .false.
         inLn3 = .false.
      Else
         inLn1 = .false.
         inLn2 = .false.
         inLn3 = .false.
      End If
      iFrst = iEnd + 1
      iEnd = Len(DDName)
*-----Check the 2nd and the 3rd field
      If (Index(DDName(iFrst:iEnd),'INLINE').ne.0) Then
*        Here if T T or F T
         inLn3 = .true.
*--------Locate the 2nd field
         Call NxtWrd(DDName,iFrst,iEnd)
*        Check the 2nd field
         If (Index(DDName(iFrst:iEnd),'INLINE').ne.0) Then
            inLn2 = .true.
         Else
            inLn2 = .false.
         End If
      End If
      If (IfTest) Write (6,*) inLn1, inLn2, inLn3
*
      If (.Not.inLn1) Then
         lUnit=11
*
*        Find and decode basis set label
*
         Call Rdbsl(DirName,BSLbl,Type,nCGTO,mCGTO,lAng,Itabmx,lUnit,
     &              dbsc(nCnttp)%AtmNr,BasisTypes,ExtBasDir)
         Line=Get_Ln(lUnit)
         Ref(1)=Line(1:80)
         Line=Get_Ln(lUnit)
         Ref(2)=Line(1:80)
      Else
         Hit=.True.
         Call Decode(BSLbl,atom,1,Hit)
         dbsc(nCnttp)%AtmNr=Lbl2Nr(atom)
         lUnit=LuRd
         Ref(1) = ''
         Ref(2) = ''
         Hit=.True.
         Call Decode(BSLBl(1:80),type,2,Hit)
         Basis_Lib=' '
         i=1
         Basis_Lib=DefNm
         Call Find_Basis_Set(Basis_Lib,' ',' ')
         i = StrnLn(Basis_Lib)
         Call BasisType(Basis_Lib(1:i)//'/'//type,0,BasisTypes)
      End If
      Line(1:3)=Type(1:3)
      Call UpCase(Line(1:3))
      If (Line(1:3).eq.'AUX') dbsc(nCnttp)%Aux=.True.
      If (IfTest) Then
         Write (6,'(A,A)') 'Ref(1):',Ref(1)
         Write (6,'(A,A)') 'Ref(2):',Ref(2)
      End If
      Uncontracted = BasisTypes(1).eq.6
      If (dbsc(nCnttp)%IsMM .eq. 1) Then
         lAng = 0
         dbsc(nCnttp)%Charge = Zero
         Return
      End If
*
*--- begin parsing options
      isEorb=.false.
      isFock=.false.
      If (L_STDINP.AND.inLn1) then ! CGGn
        iSTDINP = iSTDINP + 1      ! CGGn
        Line = STDINP(iSTDINP)     ! CGGn
      else                         ! CGGn
        Line = Get_Ln(lUnit)
      EndIf                        ! CGGn
      If(Line.eq.'Options') Then
         Do
         Line=Get_Ln(lUnit)
         If(Line.ne.'EndOptions') Then
            If(Line.eq.'OrbitalEnergies') Then
               If (IfTest)
     &         Write(6,*) 'Orbital energies are included'
               isEorb=.true.
            Else If(Line.eq.'FockOperator') Then
               If (IfTest)
     &         Write(6,*) 'Fock operator is included'
               isEorb=.true.
               isFock=.true.
            Else
               Write(6,*) 'Illegal option: ',Line
               Call Abend()
            End If
         Else
            Exit
         End If
         End Do
         Line=Get_Ln(lUnit)
      End If
*--- end parsing options
      If (IfTest) Write (6,'(A,A)') 'Line=',Line
      If (L_STDINP.AND.inLn1) then              ! CGGn
        Call Pick_Words(Line,2,Nwords,Words)    ! CGGn
        If (Nwords.NE.2) Call Abend()           ! CGGn
        Call Get_dNumber(Words(1),dbsc(nCnttp)%Charge,iErr)  ! CGGn
        If (iErr.NE.0) Call Abend()             ! CGGn
        Call Get_iNumber(Words(2),lAng,iErr)    ! CGGn
        If (iErr.NE.0) Call Abend()             ! CGGn
      else                                      ! CGGn
        call get_f1(1,dbsc(nCnttp)%Charge)
        if (inLn1) call get_i1(2,lAng)
      EndIf                                     ! CGGn
      If (iPrint.ge.99) Then
         Write (6,*) 'lAng, Charge=',lAng, dbsc(nCnttp)%Charge
         Write (6,*) ' Start reading valence basis'
      End If
      If (lAng.gt.iTabMx) Then
         Write (6,*) 'GetBS: lAng.gt.iTabMx'
         Write (6,*) 'lAng,iTabMx=',lAng,iTabMx
         Call Abend()
      End If
*     Loop over each shell type (s,p,d,etc....)
      iValSh=iShll
      dbsc(nCnttp)%nVal=lAng+1
      mVal=0
      dbsc(nCnttp)%iVal=iShll+1
      Do 10 iAng = 0, lAng
         If (IfTest) Write (6,*) 'iAng=',iAng
         iShll = iShll + 1
         If (iShll.gt.MxShll) Then
            Write (6,*) 'GetBS: iShll.gt.MxShll'
            Write (6,*) 'iShll,MxShll=',iShll,MxShll
            Call Abend()
         End If
         If (IfTest) Then
            Write (6,'(A,A)') 'Line=',Line
            Write (6,*) L_STDINP,inLn1
         End If
         If (L_STDINP.AND.inLn1) then                   ! CGGn
           iSTDINP = iSTDINP + 1                        ! CGGn
           Line = STDINP(iSTDINP)                       ! CGGn
           Call Pick_Words(Line,2,Nwords,Words)         ! CGGn
           If (Nwords.NE.2) Call Abend()                ! CGGn
           Call Get_iNumber(Words(1),nPrim,iErr)        ! CGGn
           If (iErr.NE.0) Call Abend()                  ! CGGn
           Call Get_iNumber(Words(2),nCntrc,iErr)       ! CGGn
           If (iErr.NE.0) Call Abend()                  ! CGGn
           nCGTO(iAng)=0                                ! CGGn
           mCGTO(iAng)=0                                ! CGGn
         else                                           ! CGGn
           Line = Get_Ln(lUnit)
           If (IfTest) Write (6,'(A,A)') 'Line=',Line
           Call Get_i1(1,nPrim)
           If (inLn1) then
              Call Get_i1(2,nCntrc)
              nCGTO(iAng)=0
              mCGTO(iAng)=0
           Else
              nCntrc=nCGTO(iAng)
           Endif
         EndIf                                          ! CGGn
         If (IfTest) Write(6,*) ' nPrim, nCntrc=',nPrim, nCntrc
*
         Shells(iShll)%nExp=nPrim
         Shells(iShll)%nBasis_c = nCntrc
         Call mma_allocate(Shells(iShll)%Exp,nPrim,Label='Exp')
*        Read gaussian exponents
         If (nPrim.gt.0) then
            If (IfTest) Write(6,*) 'Read gaussian exponents'
            Call Read_v(lUnit,Shells(iShll)%Exp,1,nPrim,1,Ierr)
            If (Ierr.ne.0) Then
               Call WarningMessage(2,
     &                     'GetBS: Error while reading the exponents')
               Call Quit_OnUserError()
            End If
            If (IfTest) Write(6,*) 'Done with exponents'
         End If
         If (iPrint.ge.99)
     &      Call RecPrt(' Exponents',' ',Shells(iShll)%Exp,nPrim,1)
*
*        Storage of coefficients for both contracted and uncontracted case.
*
         Call mma_allocate(Shells(iShll)%Cff_c,nPrim,nCntrc,2,
     &                     Label='Cff_c')
         Shells(iShll)%Cff_c(:,:,1)=Zero
         Call mma_allocate(Shells(iShll)%pCff,nPrim,nCntrc,
     &                     Label='pCff')
         Shells(iShll)%nBasis=nCntrc
         Call mma_allocate(Shells(iShll)%Cff_p,nPrim,nPrim,2,
     &                     Label='Cff_p')
*        Read contraction coefficients
*        Observe that the matrix will have nPrim rows and
*        nCntrc columns
         If (IfTest) Then
            Write (6,'(2A)') ' Type=',Type
            Write (6,*) 'mCGTO(iAng)=',mCGTO(iAng)
            Write (6,*) 'nCGTO(iAng)=',nCGTO(iAng)
            Write (6,*) 'nCntrc=',nCntrc
         End If
         If (IfTest) Write (6,*) ' Read/Process coefficients'
*
         If ((inLn1 .or. mCGTO(iAng).eq.nCntrc).or. nCntrc.eq.0) Then
*           Read in coeffs. in GC format, as the standard case
            If (IfTest) Write (6,*) ' Standard case'
            Shells(iShll)%Cff_c(:,:,:)=Zero
            If (UnContracted) Then
               Do i=1,nPrim
                  Shells(iShll)%Cff_c(i,i,1)=One
               End Do
            Else
               Do iPrim = 1, nPrim
                  Call Read_v(lUnit,Shells(iShll)%Cff_c(1,1,1),
     &                        iPrim,nCntrc*nPrim,nPrim,Ierr)
                  If (Ierr.ne.0) Then
                     Call WarningMessage(2,
     &                      'GetBS: Error reading coeffs in GC format')
                     Call Quit_OnUserError()
                  End If
               End Do
            End If
*
         Else
*
*           Here, we want nCntrc generally contracted functions (gcf)
*           of nPrim primitives resulting from the mCGTO(iAng) gcfs
*           of the library plus the nCntrc-mCGTO(iAng) outermost
*           primitive functions of the unextended basis set, that is,
*           of the basis set which excludes any added polarization or
*           diffuse functions.
*           (Note that everything is in generally contraction format.)
*           Example:  library:.6s.3s.  + input:.6s.5s.   -->   result
*                            x x 0                          x x 0 0 0
*                            x x 0                          x x 0 0 0
*                            x x 0                          x x 0 0 0
*                            x x 0                          x x 1 0 0
*                            x x 0                          x x 0 1 0
*                            0 0 1                          0 0 0 0 1
*
            If (IfTest) Write (6,*) ' Initial GC + outermost primitives'
            If (nCntrc.gt.nPrim) Then
               Call WarningMessage(2,
     &                     'Number of contracted more than the number'
     &             //' of primitive: correct the basis set label!')
               Call Quit_OnUserError()
            End If
            Call mma_allocate(Temp,nPrim,Max(nCntrc,mCGTO(iAng)),
     &                        Label='Temp')
            Temp(:,:)=Zero
*           read the block in the library as it is
            Do iPrim = 1, nPrim
               Call Read_v(lUnit,Temp,
     &                     iPrim,nPrim*mCGTO(iAng),nPrim,Ierr)
               If (Ierr.ne.0) Then
                  Call WarningMessage(2,
     &                        'GetBS: Error reading the block')
                  Call Quit_OnUserError()
               End If
            End Do
*
*           Order the exponents
*
            Call OrdExp1(nPrim,Shells(iShll)%Exp,mCGTO(iAng),Temp)
*
*           identify the presence of added polarization and diffusse
*           functions;
            iAdded = 0
*           Examine all the contracted functions starting with the last
      Outer:Do jNow = mCGTO(iAng), 1, -1
               iFlgOne = 0
               ! Examine the primitives
               Do iNow = 1, nPrim
                  Coeff=Temp(iNow,jNow)
                  ! Stop if it is obvious that this is not a diffuse
                  ! or polarization function.
                  If (Coeff.ne.Zero.and.Coeff.ne.One) Exit Outer
                  If (Coeff.eq.One) Then
                     If (iFlgOne.eq.1)  Cycle Outer
                     iFlgOne = 1
                  End If
               End Do
               iAdded = iAdded + 1
               If (IfTest) Write (6,*)
     &                          'function',jNow,' is an added one'
            End Do Outer
*
            nAdded = iAdded
            If (nAdded.eq.mCGTO(iAng)) nAdded=0
            If (IfTest) Write (6,*) ' nAdded=',nAdded
            If (nAdded.gt.0) Then
*              shift the added polarization and diffuse functions to
*              the right
               Do jNow = 1, nAdded
                  j1= nCntrc     -jNow+1
                  j2= mCGTO(iAng)-jNow+1
                  Do iNow = 1, nPrim
                     Temp(iNow,j1)=Temp(iNow,j2)
                  End Do
               End Do
            End If
*           insert/append the outermost primitives (in GC format)
            Do jNow = mCGTO(iAng) + 1 - nAdded, nCntrc - nAdded
               If (IfTest) write (6,*) 'jNow=',jNow
               Temp(:,jNow)=Zero
               j = jNow - (mCGTO(iAng)-nAdded)
               iPrevNow = nPrim - nAdded - (nCntrc - mCGTO(iAng))
               iNow = iPrevNow + j
               Temp(iNow,jNow) = One
            End Do
            Shells(iShll)%Cff_c(:,:,1)=Temp(:,1:nCntrc)
            Call mma_deallocate(Temp)
         End If
*
         If (IfTest) Write (6,*) ' Done! Now Process.'
*
*        Order the exponents
*
         Call OrdExp(nPrim,Shells(iShll)%Exp,nCntrc,
     &                     Shells(iShll)%Cff_c(1,1,1))
         If (nPrim*nCntrc.ne.0) mVal = mVal + 1
*
*        Decontract if integrals required in the primitive basis
*
         If (nPrim.eq.0) Go To 777
         Shells(iShll)%Cff_p(:,:,1)=Zero
         Do i=1,nPrim
            Shells(iShll)%Cff_p(i,i,1)=One
         End Do
*
*------- Save the contraction coefficients once more after the coefficients.
*        The second set will not be normalized!
*
         Shells(iShll)%pCff(:,:) = Shells(iShll)%Cff_c(:,:,1)
         Shells(iShll)%Cff_c(:,:,2)=Shells(iShll)%Cff_c(:,:,1)
         Shells(iShll)%Cff_p(:,:,2)=Shells(iShll)%Cff_p(:,:,1)
*
*        The normalization coefficients are assumed to be for
*        normalized Gaussians. In Nrmlz the contraction coefficients are
*        multiplied with the normalization coefficient of each primitive
*        Gaussian. The contracted Gaussian are then normalized with respect
*        the radial overlap.
*
         If (.Not.UnNorm) Then
            Call Nrmlz(Shells(iShll)%Exp,nPrim,
     &                 Shells(iShll)%Cff_c(1,1,1),nCntrc,iAng)
            Call Nrmlz(Shells(iShll)%Exp,nPrim,
     &                 Shells(iShll)%Cff_p(1,1,1),nPrim,iAng)
         End If
*
         If (iPrint.ge.99) Then
            nPrim = Shells(iShll)%nExp
            nCntrc= Shells(iShll)%nBasis_C
            Call RecPrt(' Coefficients (normalized)',' ',
     &                  Shells(iShll)%Cff_c(1,1,1),nPrim,nCntrc)
            Call RecPrt(' Coefficients (unnormalized)',' ',
     &                  Shells(iShll)%Cff_c(1,1,2),nPrim,nCntrc)
         End If
 777     Continue
         If (nPrim.eq.0) Go To 778
*                                                                      *
************************************************************************
*                                                                      *
*        Begin read orbital energies
*
         Call mma_allocate(Shells(iShll)%FockOp,nCntrc,nCntrc,
     &                     Label='FockOp')
         Shells(iShll)%nFockOp=nCntrc
         Shells(iShll)%FockOp(:,:)=Zero

         If (isFock) Then
            dbsc(nCnttp)%FOp=dbsc(nCnttp)%FOp .and. .True.
            Line=Get_Ln(lUnit)
            Call Get_i1(1,nEorb)
            Call mma_allocate(Temp,nEorb,nEorb,Label='Temp')
            Do i=1,nEorb
               Line=Get_Ln(lUnit)
               Call Get_F(1,Temp(1,i),nEOrb)
            End Do
            Shells(iShll)%FockOp(1:Min(nEorb,nCntrc),
     &                           1:Min(nEorb,nCntrc))
     &                = Temp(1:Min(nEorb,nCntrc),
     &                       1:Min(nEorb,nCntrc))
            Call mma_deallocate(Temp)
#ifdef _DEBUG_
            Call RecPrt('Fock',' ',Shells(iShll)%FockOp,nCntrc,nCntrc)
#endif
         Else If(isEorb) Then
            dbsc(nCnttp)%FOp=dbsc(nCnttp)%FOp .and. .True.
            Line=Get_Ln(lUnit)
            Call Get_i1(1,nEorb)
            Call mma_allocate(Temp,nEorb,1,Label='Temp')
            If(nEorb.gt.0) Then
               Line=Get_Ln(lUnit)
               Call Get_F(1,Temp(1,1),nEorb)
            End If
            Do i=1,Min(nEOrb,nCntrc)
               Shells(iShll)%FockOp(i,i)=Temp(i,1)
            End Do
            Call mma_deallocate(Temp)
#ifdef _DEBUG_
            Call RecPrt('Eorb',' ',Shells(iShll)%FockOp,nCntrc,nCntrc)
#endif
         Else
            dbsc(nCnttp)%FOp=.False.
#ifdef _DEBUG_
            Call RecPrt('Empty',' ',Shells(iShll)%FockOp,nCntrc,nCntrc)
#endif
         End If
*
*        End read orbital energies
*                                                                      *
************************************************************************
*                                                                      *
*
 778     Continue
 10   Continue
      If (mVal.eq.0) dbsc(nCnttp)%nVal=0
***************************************************************************
*-----If PAM basis set read the potentials and coefficient!
*
      If (inLn2.and. .Not.inLn1) Then
         If (IfTest) Write (6,*) ' Close library and start to read from'
     &            //' standard input'
         Close(lUnit)
         lUnit = 5
      End If
      If ( Index(BSLBl,'.PAM.').ne.0) then
         If (IfTest) Write (6,*) ' Process PAM'
         dbsc(nCnttp)%lPAM2 = .True.
         If (iPrint.ge.99) Write (6,*) ' Start reading PAMs'
         Call GetPAM(lUnit,nCnttp)
*
         If (inLn3.and. .not.inLn2) Then
            Close(lUnit)
            lUnit=5
         End If
      End If
************************************************************************
*
*-----If FRAGMENT basis set read the fragment's coordinates, basis sets,
*     orbital energies and density matrix
*
      If (inLn2.and. .Not.inLn1) Then
         If (IfTest) Write (6,*) ' Close library and start to read from'
     &            //' standard input'
         Close(lUnit)
         lUnit = 5
      End If
      If ( Index(BSLBl,'.FRAGMENT.').ne.0) then
         If (IfTest) Write (6,*) ' Process FRAGMENT'
         If (iPrint.ge.99)
     &      Write (6,*) ' Start reading fragment data'
         Call GetFragment(lUnit,nCnttp)
*
         If (inLn3.and. .not.inLn2) Then
            Close(lUnit)
            lUnit=5
         End If
      End If
************************************************************************
*
*-----If ECP basis set read the rest!
*
      If (inLn2.and. .Not.inLn1) Then
         If (IfTest) Write (6,*) ' Close library and start to read from'
     &            //' standard input'
         Close(lUnit)
         lUnit = 5
      End If
*
      nProj = -1
      nAIMP = -1
      mSOC  = -1
      If ( Index(BSLBl,'.ECP.').ne.0  .or.
     &     Index(BSLBl,'.REL.').ne.0 ) then
         If (IfTest) Write (6,*) ' Process ECPs/RELs'
         dbsc(nCnttp)%ECP = .True.
         iPrSh=iShll
         If (iPrint.ge.99)
     &      Write (6,*) ' Start reading ECPs/RELs'
         dbsc(nCnttp)%iPrj=iShll+1
         Call GetECP(lUnit,iShll,nProj,UnNorm)
         dbsc(nCnttp)%nPrj=nProj+1
*
         If (inLn3.and. .not.inLn2) Then
            Close(lUnit)
            lUnit=5
         End If
*
*--------Now read the spectral resolvent basis set
*        Line=GetLn(lUnit)
         Line=Get_Ln(lUnit)
         Call UpCase(Line)
         If (Index(Line,'SPEC').eq.0) Then
            Call WarningMessage(2,
     &                 'ERROR: Keyword SPECTRAL expected,'
     &                 //' offending line : '//Line)
            Call Quit_OnUserError()
         Endif
*
         iMPShll=iShll
 9988    Continue
         Line=Get_Ln(lUnit)
         Call UpCase(Line)
         If (Line(1:4).eq.'END ') Go To 999
         If (Line(1:4).eq.'VALE') Go To 1001
         If (Line(1:4).eq.'CORE') Go To 1011
         If (Line(1:4).eq.'EXTE') Go To 1012
         If (Line(1:4).eq.'EXCH') Go To 1002
         If (Line(1:4).eq.'1STO') Go To 1003
         If (Line(1:4).eq.'NOPA') Go To 1005
*
         If (Line(1:4).eq.'NOP1') Go To 1006
         If (Line(1:4).eq.'NOP2') Go To 1007
         If (Line(1:4).eq.'NOP3') Go To 1008
         If (Line(1:4).eq.'NOPF') Go To 1010
         If (Line(1:4).eq.'RESC') Go To 1009
         If (Line(1:4).eq.'RA0H') Go To 9001
         If (Line(1:4).eq.'RA0F') Go To 9002
         If (Line(1:4).eq.'RAIH') Go To 9003
*
         If (Line(1:4).eq.'MIXE') Go To 1013
         If (Line(1:4).eq.'SOC ') Go To 1014
         If (Line(1:4).eq.'DKSO') Go To 1015
         Call WarningMessage(2,' Invalid keyword in GetBS;'//
     &                   Line)
         Call Quit_OnUserError()
*
*--------Valence basis set
*
 1001    Continue
         If (nAIMP.ne.-1) Then
            Call WarningMessage(2,
     &                  ' SR basis set is already defined!')
            Call Quit_OnUserError()
         End If
         nAIMP = lAng
         dbsc(nCnttp)%nSRO=nAIMP+1
         dbsc(nCnttp)%iSRO=iShll+1
         jValSh=iValSh
         Do iAIMP = 0, nAIMP
            iShll = iShll + 1
            If (iShll.gt.MxShll) Then
               Write (6,*) 'GetBS: iShll.gt.MxShll'
               Write (6,*) 'iShll,MxShll=',iShll,MxShll
               Call Abend()
            End If
            jValSh = jValSh + 1
            Call mma_allocate(Shells(iShll)%Exp,Shells(jValSh)%nExp,
     &                        Label='Exp')
            Shells(iShll)%Exp(:) = Shells(jValSh)%Exp(:)
            Shells(iShll)%nExp = Shells(jValSh)%nExp
            Shells(iShll)%nBasis  = 0
         End Do
         Go To 9988
*
*--------Core basis set
*
 1011    Continue
         If (nAIMP.ne.-1) Then
            Call WarningMessage(2,
     &                  ' SR basis set is already defined!')
            Call Quit_OnUserError()
         End If
         nAIMP = nProj
         dbsc(nCnttp)%nSRO=nAIMP+1
         dbsc(nCnttp)%iSRO=iShll+1
         jPrSh = iPrSh
         Do iAIMP = 0, nAIMP
            iShll = iShll + 1
            If (iShll.gt.MxShll) Then
               Write (6,*) 'GetBS: iShll.gt.MxShll'
               Write (6,*) 'iShll,MxShll=',iShll,MxShll
               Call Abend()
            End If
            jPrSh = jPrSh + 1
            Call mma_allocate(Shells(iShll)%Exp,Shells(jPrSh)%nExp,
     &                        Label='Exp')
            Shells(iShll)%Exp(:)=Shells(jPrSh)%Exp(:)
            Shells(iShll)%nExp=Shells(jPrSh)%nExp
            Shells(iShll)%nBasis  = 0
         End Do
         Go To 9988
*
*--------External basis set
*
 1012    Continue
         If (nAIMP.ne.-1) Then
            Call WarningMessage(2,
     &                  ' SR basis set is already defined!')
            Call Quit_OnUserError()
         End If
*        Line = GetLn(lUnit)
         Line = Get_Ln(lUnit)
         Call Get_i1(1,nAIMP)
         dbsc(nCnttp)%nSRO=nAIMP+1
         dbsc(nCnttp)%iSRO=iShll+1
         Do iAIMP = 0, nAIMP
            iShll = iShll + 1
            If (iShll.gt.MxShll) Then
               Write (6,*) 'GetBS: iShll.gt.MxShll'
               Write (6,*) 'iShll,MxShll=',iShll,MxShll
               Call Abend()
            End If
            Line = Get_Ln(lUnit)
            Call Get_i1(1,nPrim)
            Call mma_allocate(Shells(iShll)%Exp,nPrim,Label='Exp')
            Shells(iShll)%nExp=nPrim
            Shells(iShll)%nBasis=0
*
            If (nPrim.gt.0) then
               Call read_v(lUnit,Shells(iShll)%Exp,1,nPrim,1,Ierr)
               If (Ierr.ne.0) Then
                  Call WarningMessage(2,
     &                        'GetBS: Error reading SRO exponents')
                  Call Quit_OnUserError()
               End If
            End If
*
         End Do
         Go To 9988
*---------- Mixed basis set (valence + core), with dominance of the
*           valence basis set, i.e. the valence basis set plus those
*           exponents in the core primitives not too close to the
*           valence ones.

 1013    Continue
         If (nAIMP.ne.-1) Then
            Call WarningMessage(2,
     &                  ' SR basis set is already defined!')
            Call Quit_OnUserError()
         End If

*---------- Threshold for the ratio between consecutive exponents
*           (if the ratio between two consecutive exponents is lower
*           than the threshold, the exponent of the core basis set
*           will be removed)

         Line = Get_Ln(lUnit)
         Call Get_F1(1,RatioThres)
*
         nAIMP  = lAng
         dbsc(nCnttp)%iSRO=iShll+1
         dbsc(nCnttp)%nSRO=nAIMP+1
         jValSh = iValSh
         jPrSh  = iPrSh
*
         Do iAIMP = 0, nAIMP
            iShll = iShll + 1
            If (iShll.gt.MxShll) Then
               Write (6,*) 'GetBS: iShll.gt.MxShll'
               Write (6,*) 'iShll,MxShll=',iShll,MxShll
               Call Abend()
            End If
*
            jValSh = jValSh + 1
            nCntrc   = Shells(jValSh)%nBasis
*
*
            If (iAIMP.le.nProj) Then
               jPrSh = jPrSh + 1
*
               iDominantSet = 2
               Call mma_allocate(ExpMerged,
     &                           Shells(jPrSh)%nExp
     &                          +Shells(jValSh)%nExp,
     &                           Label='ExpMerged')
               Call MergeBS (Shells(jPrSh)%Exp,Shells(jPrSh)%nExp,
     &                       Shells(jValSh)%Exp,Shells(jValSh)%nExp,
     &                       ExpMerged,Shells(iShll)%nExp, RatioThres,
     &                       iDominantSet)
               Call mma_allocate(Shells(iShll)%Exp,Shells(iShll)%nExp,
     &                           Label='Exp')
               Shells(iShll)%Exp(:)=ExpMerged(1:Shells(iShll)%nExp)
               Call mma_deallocate(ExpMerged)

*
            Else
*
               Shells(iShll)%nExp=Shells(jValSh)%nExp
               Call mma_allocate(Shells(iShll)%Exp,Shells(iShll)%nExp,
     &                           Label='Exp')
               Shells(iShll)%Exp(:)=Shells(jValSh)%Exp(:)
*
            End If
*
            Shells(iShll)%nBasis=0
*
         End Do
         Go To 9988
*
*---------- SOC basis set
*
 1014    Continue
*
         If (mSOC.ne.-1) Then
            Call WarningMessage(2,
     &                  ' SOC basis set is already defined!')
            Call Quit_OnUserError()
         End If
*
      dbsc(nCnttp)%iSOC=iShll+1
      Line = Get_Ln(lUnit)
      Call Get_I1(1,mSOC)
      dbsc(nCnttp)%nSOC=mSOC+1
      IF (IfTest) Write(6,'(A,I4)') 'dbsc(nCnttp)%nSOC =',
     &                               dbsc(nCnttp)%nSOC
      If (mSOC.lt.0) Go To 990
      Do 12 iAng = 0, mSOC
         If (IfTest) Write (6,'(A,I4)') ' iAng=',iAng
         iShll = iShll + 1
         If (iShll.gt.MxShll) Then
            Call WarningMessage(2,
     &                      'Abend in GetBS: Increase MxShll')
            Call Quit_OnUserError()
         End If
         Line = Get_Ln(lUnit)
         Call Get_I1(1,nPrim)
         Call Get_I1(2,nCntrc)
         Call Get_I1(3,mDel)
         dbsc(nCnttp)%kDel(iAng)=mDel
         If (IfTest) Write(6,*) 'nPrim = ',nPrim,' nCntrc = ',nCntrc
         If (IfTest) Write(6,*) 'nDeleted = ', mDel
         Call mma_allocate(Shells(iShll)%Exp,nPrim,Label='Exp')
         Shells(iShll)%nExp=nPrim
         Shells(iShll)%nBasis = nCntrc
         If (IfTest) Write (6,*) 'getBS: ishll,nCntrc',ishll,nCntrc
         If (IfTest) Write (6,'(A)') ' Reading Exponents'
         If (nPrim.gt.0) Call Read_v(lUnit,Shells(iShll)%Exp,1,nPrim,1,
     &                               ierr)
         If (IfTest)
     &      Call RecPrt('Exponents',' ',Shells(iShll)%Exp,1,nPrim)
         Call mma_allocate(Shells(iShll)%Cff_c,nPrim,nCntrc,2,
     &                     Label='Cff_c')
         Call mma_allocate(Shells(iShll)%pCff,nPrim,nCntrc,
     &                     Label='pCff')
         Shells(iShll)%nBasis=nCntrc
         Call mma_allocate(Shells(iShll)%Cff_p,nPrim,nPrim,2,
     &                     Label='Cff_p')
         Shells(iShll)%Cff_p(:,:,:)=Zero
         If (IfTest) Write (6,'(A)') ' Reading coefficients'
         Do 20 iPrim = 1, nPrim
            Call Read_v(lUnit,Shells(iShll)%Cff_c(1,1,1),
     &                  iPrim,nPrim*nCntrc,nPrim,ierr)
 20      Continue
*
         Shells(iShll)%pCff(:,:) = Shells(iShll)%Cff_c(:,:,1)
         Shells(iShll)%Cff_c(:,:,2)=Shells(iShll)%Cff_c(:,:,1)
         Shells(iShll)%Cff_p(:,:,2)=Shells(iShll)%Cff_p(:,:,1)
*
         If (IfTest)
     &      Call RecPrt('Coefficients',' ',Shells(iShll)%Cff_c(1,1,1),
     &                  nPrim,nCntrc)
*
 12   Continue
 990     Continue
         Go To 9988
*
*------  Use DKSO on request
*
 1015    Continue
         dbsc(nCnttp)%SODK=.True.
         Go To 9988
*
*--------Exchange operator
*
 1002    Continue
         dbsc(nCnttp)%nOpt = iOr(dbsc(nCnttp)%nOpt,2**0)
         Go To 9988
*
*--------1st order relativistic correction
*
 1003    Continue
         dbsc(nCnttp)%nOpt = iOr(dbsc(nCnttp)%nOpt,2**1)
         dbsc(nCnttp)%nOpt = iOr(dbsc(nCnttp)%nOpt,2**2)
         Line=Get_Ln(lUnit)
         MPLbl=Line(1:20)
         Go To 9988
*
*
*        one-centre no-pair operators
*
 1005    Continue
         dbsc(nCnttp)%NoPair=.True.
         dbsc(nCnttp)%SODK=.True.
         IRELMP=0
         dbsc(nCnttp)%nOpt = iOr(dbsc(nCnttp)%nOpt,2**3)
         Go To 9988
*
*        one-centre no-pair operators (DK1)
*
 1006    Continue
         dbsc(nCnttp)%NoPair=.True.
         dbsc(nCnttp)%SODK=.True.
         IRELMP=1
         dbsc(nCnttp)%nOpt = iOr(dbsc(nCnttp)%nOpt,2**3)
         Go To 9988
*
*        one-centre no-pair operators (DK2)
*
 1007    Continue
         dbsc(nCnttp)%NoPair=.True.
         dbsc(nCnttp)%SODK=.True.
         IRELMP=2
         dbsc(nCnttp)%nOpt = iOr(dbsc(nCnttp)%nOpt,2**3)
         Go To 9988
*
*        one-centre no-pair operators (DK3)
*
 1008    Continue
         dbsc(nCnttp)%NoPair=.True.
         dbsc(nCnttp)%SODK=.True.
         IRELMP=3
         dbsc(nCnttp)%nOpt = iOr(dbsc(nCnttp)%nOpt,2**3)
         Go To 9988
*
*        one-centre no-pair operators (DK3)
*
 1010    Continue
         dbsc(nCnttp)%NoPair=.True.
         dbsc(nCnttp)%SODK=.True.
         IRELMP=4
         dbsc(nCnttp)%nOpt = iOr(dbsc(nCnttp)%nOpt,2**3)
         Go To 9988
*
*        one-centre RESC operators
*
 1009    Continue
         dbsc(nCnttp)%NoPair=.True.
         IRELMP=11
         dbsc(nCnttp)%nOpt = iOr(dbsc(nCnttp)%nOpt,2**3)
         Go To 9988
*
*        one-centre ZORA operators
*
 9001    Continue
         dbsc(nCnttp)%NoPair=.True.
         IRELMP=21
         dbsc(nCnttp)%nOpt = iOr(dbsc(nCnttp)%nOpt,2**3)
         Go To 9988
*
*        one-centre ZORA-FP operators
*
 9002    Continue
         dbsc(nCnttp)%NoPair=.True.
         dbsc(nCnttp)%SODK=.True.
         IRELMP=22
         dbsc(nCnttp)%nOpt = iOr(dbsc(nCnttp)%nOpt,2**3)
         Go To 9988
*
*        one-centre IORA operators
*
 9003    Continue
         dbsc(nCnttp)%NoPair=.True.
         dbsc(nCnttp)%SODK=.True.
         IRELMP=23
         dbsc(nCnttp)%nOpt = iOr(dbsc(nCnttp)%nOpt,2**3)
         Go To 9988
*
 999     Continue
         If (iAnd(dbsc(nCnttp)%nOpt,2**1).ne.0 .and.
     &       iAnd(dbsc(nCnttp)%nOpt,2**3).ne.0) Then
            Call WarningMessage(2,
     &                  ' 1st order relativistic correction and '
     &                //' no-pair approximation can not be used'
     &                //' simultaneously!')
            Call Quit_OnUserError()
         End If
         If (nAIMP.ge.0) Then
            Basis_Lib=DefNm
            Call Find_Basis_Set(Basis_Lib,' ',' ')
            Filename=Trim(Basis_Lib)//'/QRPLIB'
            Call f_Inquire(Filename,Found)
            If (.not. Found) Then
               Write(6,*) 'File '//Trim(Filename)//' not found'
               call abend()
            End If
            LUQRP=33
            call molcas_open(LUQRP,Filename)
c            Open(LUQRP,file='QRPLIB',form='formatted')
            Call CalcAMt(dbsc(nCnttp)%nOpt,LUQRP,MPLbl,nAIMP,iMPShll+1,
     &                   nProj,iPrSh+1,DBLE(dbsc(nCnttp)%AtmNr))
            Close (LUQRP)
         End If
      End If
*
      lAng = Max(lAng,nProj,nAIMP)
      If (.not.inLn3) Close(lUnit)
      Return
      End
      Subroutine Check_Info()
#include "itmax.fh"
#include "info.fh"
      Call Free_Work(LctInf)
      Call Abend()
      End Subroutine Check_Info
