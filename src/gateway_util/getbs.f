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
      SubRoutine GetBS(DDname,BSLbl,iBSLbl,
     &                 lAng,ipExp,ipCff,ipCff_Cntrct,ipCff_Prim,
     &                 nExp,nBasis,nBasis_Cntrct,MxShll,iShll,
     &                 MxAng, Charge,iAtmNr,BLine,Ref,
     &                 PAM2,FockOp, ECP,NoPairL,SODK,
     &                 CrRep,nProj,nAIMP,iOpt,
     &                 UnNorm,nDel,
     &                  nVal,  nPrj,  nSRO,  nSOC, nPP,
     &                 ipVal_,ipPrj_,ipSRO_,ipSOC_,ipPP_,LuRd,
     &                 BasisTypes,AuxCnttp, IsMM,
     &                 STDINP,iSTDINP,L_STDINP,Expert,ExtBasDir,
     &                 DInf,nDInf,nCnttp)
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
      Use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "itmax.fh"
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
      Real*8 DInf(nDInf)
      Character*80 BSLbl, BLine, Ref(2), MPLbl*20,
     &             Filenm, Atom, Type
      Character*256 DirName
*
      Character Basis_Lib*256, Filename*263, DefNm*13
      Integer StrnLn
      External StrnLn
      Logical UnContracted, L_STDINP
*
      Character*180 Line, Get_Ln, STDINP(MxAtom*2) ! CGGn
      External Get_Ln
      Character*(*) DDname
      Character*24 Words(2)                     ! CGGn
      Logical ECP, inLn1, inLn2, inLn3, Hit, IfTest,NoPairL,
     &        UnNorm, PAM2, SODK, AuxCnttp, FockOp,
     &        isEorb,isFock
      Integer ipExp(MxShll), ipCff(MxShll), ipCff_Cntrct(MxShll),
     &        ipCff_Prim(MxShll),
     &        nExp(MxShll), nBasis(MxShll), nCGTO(0:iTabMx),
     &        mCGTO(0:iTabMx), nDel(0:MxAng),
     &        nBasis_Cntrct(MxShll)
      Integer BasisTypes(4)
      Logical Expert, Found
      Character *(*) ExtBasDir
      Real*8, Allocatable:: ExpMerged(:)
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
*
*                                                                      *
************************************************************************
*                                                                      *
      Interface
         SubRoutine GetECP(lUnit,ipExp,ipCff,nExp,nBasis,MxShll,iShll,
     &                     BLine,CrRep,nProj,
     &                     ipPP,nPP,UnNorm,DInf,nDInf,nCnttp)
         Integer lUnit
         Integer ipExp(MxShll), ipCff(MxShll),
     &            nExp(MxShll), nBasis(MxShll)
         Integer MxShll,iShll
         Character*(*) BLine
         Real*8  CrRep
         Integer nProj
         Integer ipPP, nPP
         Logical UnNorm
         Real*8  DInf(nDInf)
         Integer nCnttp
         End SubRoutine GetECP
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
      iRout=6
      iPrint = nPrint(iRout)
      IfTest=.False.
!#define _DEBUG_
#ifdef _DEBUG_
      IfTest=.True.
      iPrint=99
#endif
      If (IfTest) iPrint=99
      ip_Dummy=-1
      PAM2 = .False.
      ECP = .False.
      FockOp = .True.
      NoPairL = .False.
      nVal=0
      nPrj=0
      nSRO=0
      nSOC=0
      nPP=0
      nM1=0
      nM2=0
      ipVal_=-1
      ipPrj_=-1
      ipSRO_=-1
      ipSOC_=-1
      ipPP_=-1
      IsMM = 0
      iOpt = 0
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
         IsMM = 1
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
     &              iAtmNr,BasisTypes,ExtBasDir)
         Line=Get_Ln(lUnit)
         Ref(1)=Line(1:80)
         Line=Get_Ln(lUnit)
         Ref(2)=Line(1:80)
      Else
         Hit=.True.
         Call Decode(BSLbl,atom,1,Hit)
         iAtmNr=Lbl2Nr(atom)
         lUnit=LuRd
         Ref(1) = BLine
         Ref(2) = BLine
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
      If (Line(1:3).eq.'AUX') AuxCnttp=.True.
      If (IfTest) Then
         Write (6,'(A,A)') 'Ref(1):',Ref(1)
         Write (6,'(A,A)') 'Ref(2):',Ref(2)
      End If
      Uncontracted = BasisTypes(1).eq.6
      If (IsMM .eq. 1) Then
         lAng = 0
         Charge = Zero
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
107      Continue
         Line=Get_Ln(lUnit)
         If(Line.ne.'EndOptions') Then
            If(Line.eq.'OrbitalEnergies') Then
C              Write(6,*) 'Orbital energies are included'
               isEorb=.true.
            Else If(Line.eq.'FockOperator') Then
C              Write(6,*) 'Fock operator is included'
               isEorb=.true.
               isFock=.true.
            Else
               Write(6,*) 'Illegal option: ',Line
            End If
            Go To 107
         End If
         Line=Get_Ln(lUnit)
      End If
*--- end parsing options
      If (IfTest) Write (6,'(A,A)') 'Line=',Line
      If (L_STDINP.AND.inLn1) then              ! CGGn
        Call Pick_Words(Line,2,Nwords,Words)    ! CGGn
        If (Nwords.NE.2) Call Abend()           ! CGGn
        Call Get_dNumber(Words(1),Charge,iErr)  ! CGGn
        If (iErr.NE.0) Call Abend()             ! CGGn
        Call Get_iNumber(Words(2),lAng,iErr)    ! CGGn
        If (iErr.NE.0) Call Abend()             ! CGGn
      else                                      ! CGGn
        call get_f1(1,Charge)
        if (inLn1) call get_i1(2,lAng)
      EndIf                                     ! CGGn
      If (iPrint.ge.99) Then
         Write (6,*) 'lAng, Charge=',lAng, Charge
         Write (6,*) ' Start reading valence basis'
      End If
      If (lAng.gt.MxAng) Then
         Write (6,*) 'GetBS: lAng.gt.MxAng'
         Write (6,*) 'lAng,MxAng=',lAng,MxAng
         Call Abend()
      End If
*     Loop over each shell type (s,p,d,etc....)
      iValSh=iShll
      iStrt1=ipExp(iShll+1)
      nVal=lAng+1
      mVal=0
      ipVal_=iShll+1
      Do 10 iAng = 0, lAng
         If (IfTest) Write (6,*) 'iAng=',iAng
         iShll = iShll + 1
         If (iShll.gt.MxShll) Then
            Write (6,*) 'GetBS: iShll.gt.MxShll'
            Write (6,*) 'iShll,MxShll=',iShll,MxShll
            Call Abend()
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
           Call Get_i1(1,nPrim)
           If (inLn1) then
              Call Get_i1(2,nCntrc)
              nCGTO(iAng)=0
              mCGTO(iAng)=0
           Else
              nCntrc=nCGTO(iAng)
           Endif
         EndIf                                          ! CGGn
         nDntrc=nCntrc
         If (IfTest) Write(6,*) ' nPrim, nCntrc=',nPrim, nCntrc
*
         iStrt=ipExp(iShll)
         nExp(iShll) = nPrim
         nBasis_Cntrct(iShll) = nCntrc
         iEnd = iStrt - 1
         Call mma_allocate(Shells(iShll)%Exp,nPrim,Label='Exp')
         Shells(iShll)%nExp=nPrim
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
         iStrt = iEnd + 1
*
*        Storage of coefficients for both contracted and uncontracted case.
*
         ipCff_c = iStrt
         ipCff_Cntrct(iShll)=iStrt
         iEnds= iEnd + 2*nPrim*nCntrc
         ipCff_Prim(iShll)= iEnds + 1
         ipCff_p = ipCff_Prim(iShll)
         iEnds= iEnds+ 2*nPrim**2
culf
         iEndc = iStrt + nPrim*nDntrc - 1
         iEnd  = iStrt + nPrim*nCntrc - 1
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
         If ((inLn1 .or. mCGTO(iAng).eq.nCntrc).or.
     &       nCntrc.eq.0) Then
*           Read in coeffs. in GC format, as the standard case
            If (IfTest) Write (6,*) ' Standard case'
            Call FZero(DInf(iStrt),nPrim*nDntrc)
            If (UnContracted) Then
               Call DCopy_(nPrim,[One],0,DInf(iStrt),nPrim+1)
            Else
               Do iPrim = 0, nPrim-1
                  Call Read_v(lUnit,DInf,iStrt+iPrim,iEndc,nPrim,Ierr)
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
*           read the block in the library as it is
            iEndNow = iStrt + nPrim*mCGTO(iAng) - 1
            Do iPrim = 0, nPrim-1
               Call Read_v(lUnit,DInf,iStrt+iPrim,iEndNow,nPrim,Ierr)
               If (Ierr.ne.0) Then
                  Call WarningMessage(2,
     &                        'GetBS: Error reading the block')
                  Call Quit_OnUserError()
               End If
            End Do
*
*           Order the exponents
*
            Call OrdExp1(nPrim,Shells(iShll)%Exp,
     &                   mCGTO(iAng),DInf(ipCff_c))
*
*           identify the presence of added polarization and diffusse
*           functions;
            iAdded = 0
            Do 212 jNow = mCGTO(iAng), 1, -1
               iFlgOne = 0
               Do iNow = 1, nPrim
                  i = iStrt - 1 + nPrim*(jNow-1)+iNow
                  If (DInf(i).ne.Zero.and.DInf(i).ne.One) go to 2129
                  If (DInf(i).eq.One) Then
                     If (iFlgOne.eq.1)  go to 212
                     iFlgOne = 1
                  End If
               End Do
               iAdded = iAdded + 1
               If (IfTest) Write (6,*)
     &                          'function',jNow,' is an added one'
212         Continue
2129        Continue
            nAdded = iAdded
            If (nAdded.eq.mCGTO(iAng)) nAdded=0
            If (IfTest) Write (6,*) ' nAdded=',nAdded
            If (nAdded.gt.0) Then
*              shift the added polarization and diffuse functions to
*              the right
               Do jNow = 1, nAdded
                  Do iNow = 1, nPrim
                     j = nCntrc     -jNow+1
                     iNew = iStrt - 1 + nPrim*(j-1)+iNow
                     j = mCGTO(iAng)-jNow+1
                     iOld = iStrt - 1 + nPrim*(j-1)+iNow
                     DInf(iNew) = DInf(iOld)
                  End Do
               End Do
            End If
*           insert/append the outermost primitives (in GC format)
            Do jNow = mCGTO(iAng) + 1 - nAdded, nCntrc - nAdded
               If (IfTest) write (6,*) 'jNow=',jNow
               iPrevNow = iStrt - 1 + nPrim*(jNow-1)
*              write (*,*) '0:',nPrim*(jNow-1)+1,nPrim*(jNow-1)+nPrim
               Do i = iPrevNow + 1, iPrevNow + nPrim
                  DInf(i) = Zero
               End Do
               j = jNow - (mCGTO(iAng)-nAdded)
               iPrevNow = nPrim - nAdded - (nCntrc - mCGTO(iAng))
               iNow = iPrevNow + j
               i = iStrt - 1 + nPrim*(jNow-1)+iNow
*              write (*,*) 'j=',j,'i=',nPrim*(jNow-1)+iNow
               DInf(i) = One
            End Do
         End If
*
         If (IfTest) Write (6,*) ' Done! Now Process.'
*
*        Order the exponents
*
         Call OrdExp(nPrim,Shells(iShll)%Exp,nCntrc,DInf(ipCff_c))
         If (nPrim*nCntrc.ne.0) mVal = mVal + 1
*
*        Decontract if integrals required in the primitive basis
*
         If (nPrim.eq.0) Go To 777
         Call DCopy_(nPrim*nPrim,[Zero],0,DInf(ipCff_p),1)
         Call DCopy_(nPrim,[One],0,DInf(ipCff_p),nPrim+1)
*
*------- Save the contraction coefficients once more after the coefficients.
*        The second set will not be normalized!
*
         iOff = nPrim*nCntrc
         Call DCopy_(nPrim*nCntrc,DInf(ipCff_c),1,DInf(ipCff_c+iOff),1)
         iOff = nPrim*nPrim
         Call DCopy_(nPrim*nPrim, DInf(ipCff_p),1,DInf(ipCff_p+iOff),1)
*
*        The normalization coefficients are assumed to be for
*        normalized Gaussians. In Nrmlz the contraction coefficients are
*        multiplied with the normalization coefficient of each primitive
*        Gaussian. The contracted Gaussian are then normalized with respect
*        the radial overlap.
*
         If (.Not.UnNorm) Then
            Call Nrmlz(Shells(iShll)%Exp,nPrim,
     &                 DInf(ipCff_c),nCntrc,iAng)
            Call Nrmlz(Shells(iShll)%Exp,nPrim,
     &                 DInf(ipCff_p),nPrim,iAng)
         End If
*
         If (iPrint.ge.99) Then
            ipCff_x = ipCff_Cntrct(iShll)
            nPrim = nExp(iShll)
            nCntrc= nBasis_Cntrct(iShll)
            Call RecPrt(' Coefficients (normalized)',' ',
     &                  DInf(ipCff_x),nPrim,nCntrc)
            Call RecPrt(' Coefficients (unnormalized)',' ',
     &                  DInf(ipCff_x+nPrim*nCntrc),nPrim,nCntrc)
         End If
 777     Continue
         iEnd = iEnds
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
            FockOp=FockOp .and. .True.
            Line=Get_Ln(lUnit)
            Call Get_i1(1,nEorb)
            Do i=1,nEorb
               Line=Get_Ln(lUnit)
               Call Get_F(1,Shells(iShll)%FockOp(1,i),nEOrb)
            End Do
#ifdef _DEBUG_
            Call RecPrt('Fock',' ',Shells(iShll)%FockOp,nCntrc,nCntrc)
#endif
         Else If(isEorb) Then
            FockOp=FockOp .and. .True.
            Line=Get_Ln(lUnit)
            Call Get_i1(1,nEorb)
            If(nEorb.gt.0) Then
               Line=Get_Ln(lUnit)
               Call Get_F(1,Shells(iShll)%FockOp,nEorb)
            End If
            Do i=2,nCntrc
               Shells(iShll)%FockOp(i,i)=Shells(iShll)%FockOp(i,1)
               Shells(iShll)%FockOp(i,1)=Zero
            End Do
#ifdef _DEBUG_
            Call RecPrt('Eorb',' ',Shells(iShll)%FockOp,nCntrc,nCntrc)
#endif
         Else
            FockOp=.False.
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
         If (iShll.lt.MxShll) ipExp(iShll+1) = iEnd + 1
 10   Continue
      If (mVal.eq.0) nVal=0
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
         PAM2 = .True.
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
         ECP = .True.
         iPrSh=iShll
         If (iPrint.ge.99)
     &      Write (6,*) ' Start reading ECPs/RELs'
         ipPrj_=iShll+1
         Call GetECP(lUnit,ipExp,ipCff,nExp,nBasis,MxShll,iShll,Bline,
     &               CrRep,nProj,ipPP_,nPP,UnNorm,
     &               DInf,nDinf,nCnttp)
         nPrj=nProj+1
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
         nSRO=nAIMP+1
         ipSRO_=iShll+1
         jValSh=iValSh
         Do iAIMP = 0, nAIMP
            iShll = iShll + 1
            iStrt= ipExp(iShll)
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
            nExp(iShll)  = nExp(jValSh)
            nBasis(iShll)  = 0
            ipCff(iShll)   = ip_Dummy
            iEnd = iStrt - 1
            If (iShll.lt.MxShll) ipExp(iShll+1) = iEnd + 1
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
         nSRO=nAIMP+1
         ipSRO_=iShll+1
         jPrSh = iPrSh
         Do iAIMP = 0, nAIMP
            iShll = iShll + 1
            iStrt= ipExp(iShll)
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
            nExp(iShll)  = nExp(jPrSh)
            nBasis(iShll)  = 0
            ipCff(iShll)   = ip_Dummy
            iEnd = iStrt - 1
            If (iShll.lt.MxShll) ipExp(iShll+1) = iEnd + 1
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
         nSRO=nAIMP+1
         ipSRO_=iShll+1
         Do iAIMP = 0, nAIMP
            iShll = iShll + 1
            iStrt= ipExp(iShll)
            If (iShll.gt.MxShll) Then
               Write (6,*) 'GetBS: iShll.gt.MxShll'
               Write (6,*) 'iShll,MxShll=',iShll,MxShll
               Call Abend()
            End If
            Line = Get_Ln(lUnit)
            Call Get_i1(1,nPrim)
            iStrt = ipExp(iShll)
            Call mma_allocate(Shells(iShll)%Exp,nPrim,Label='Exp')
            Shells(iShll)%nExp=nPrim
            nExp(iShll) = nPrim
            nBasis(iShll) = 0
            ipCff(iShll)   = ip_Dummy
            iEnd = iStrt - 1
*
            If (nPrim.gt.0) then
               Call read_v(lUnit,Shells(iShll)%Exp,1,nPrim,1,Ierr)
               If (Ierr.ne.0) Then
                  Call WarningMessage(2,
     &                        'GetBS: Error reading SRO exponents')
                  Call Quit_OnUserError()
               End If
            End If
            iStrt = iEnd + 1
            iEnd = iStrt
*
            If (iShll.lt.MxShll) ipExp(iShll+1) = iEnd + 1
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
         ipSRO_=iShll+1
         nSRO=nAIMP+1
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
            nCntrc   = nBasis(jValSh)
*
            iStrt = ipExp(iShll)
*
            If (iAIMP.le.nProj) Then
               jPrSh = jPrSh + 1
*
               iDominantSet = 2
               Call mma_allocate(ExpMerged,nExp(jPrSh)+nExp(jValSh),
     &                           Label='ExpMerged')
               Call MergeBS (Shells(jPrSh)%Exp,nExp(jPrSh),
     &                       Shells(jValSh)%Exp,nExp(jValSh),
     &                       ExpMerged,nExp(iShll), RatioThres,
     &                       iDominantSet)
               Call mma_allocate(Shells(iShll)%Exp,nExp(iShll),
     &                           Label='Exp')
               Shells(iShll)%Exp(:)=ExpMerged(1:nExp(iShll))
               Shells(iShll)%nExp=nExp(iShll)
               Call mma_deallocate(ExpMerged)

*
            Else
*
               nExp(iShll) = nExp(jValSh)
               Call mma_allocate(Shells(iShll)%Exp,nExp(iShll),
     &                           Label='Exp')
               Shells(iShll)%Exp(:)=Shells(jValSh)%Exp(:)
               Shells(iShll)%nExp=nExp(iShll)
*
            End If
*
            iEnd = iStrt - 1
            nBasis(iShll) = 0
            ipCff(iShll)   = ip_Dummy
*
            iStrt = iEnd + 1
            iEnd = iStrt
*
            If (iShll.lt.MxShll) ipExp(iShll+1) = iEnd + 1
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
      ipSoc_=iShll+1
      Line = Get_Ln(lUnit)
      Call Get_I1(1,mSOC)
      nSOC=mSOC+1
      IF (IfTest) Write(6,'(A,I4)') 'nSOC =',nSOC
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
         nDel(iAng)=mDel
         If (IfTest) Write(6,*) 'nPrim = ',nPrim,' nCntrc = ',nCntrc
         If (IfTest) Write(6,*) 'nDeleted = ', mDel
         iStrt = ipExp(iShll)
         Call mma_allocate(Shells(iShll)%Exp,nPrim,Label='Exp')
         Shells(iShll)%nExp=nPrim
         nExp(iShll) = nPrim
         nBasis(iShll) = nCntrc
         If (IfTest) Write (6,*) 'getBS: ishll,nCntrc',ishll,nCntrc
         If (IfTest) Write (6,'(A)') ' Reading Exponents'
         iEnd = iStrt + nPrim - 1
         If (nPrim.gt.0) Call Read_v(lUnit,Shells(iShll)%Exp,1,nPrim,1,
     &                               ierr)
         If (IfTest)
     &      Call RecPrt('Exponents',' ',Shells(iShll)%Exp,1,nPrim)
         iStrt = iEnd + 1
         ipCff(iShll) = iStrt
         iEnd = iStrt + nPrim*nCntrc - 1
         If (IfTest) Write (6,'(A)') ' Reading coefficients'
         Do 20 iPrim = 0, nPrim-1
            Call Read_v(lUnit,DInf,iStrt+iPrim,iEnd,nPrim,ierr)
 20      Continue
         If (IfTest) Call RecPrt('Exponents',' ',DInf(ipCff(iShll)),
     &                           nPrim,nCntrc)
         iStrt = iEnd + 1
         If (iShll.lt.MxShll) ipExp(iShll+1) = iEnd + 1
*
 12   Continue
 990     Continue
         Go To 9988
*
*------  Use DKSO on request
*
 1015    Continue
         SODK=.True.
         Go To 9988
*
*--------Exchange operator
*
 1002    Continue
         iOpt = iOr(iOpt,2**0)
         Go To 9988
*
*--------1st order relativistic correction
*
 1003    Continue
         iOpt = iOr(iOpt,2**1)
         iOpt = iOr(iOpt,2**2)
         Line=Get_Ln(lUnit)
         MPLbl=Line(1:20)
         Go To 9988
*
*
*        one-centre no-pair operators
*
 1005    Continue
         NoPairL=.True.
         SODK=.True.
         IRELMP=0
         iOpt = iOr(iOpt,2**3)
         Go To 9988
*
*        one-centre no-pair operators (DK1)
*
 1006    Continue
         NoPairL=.True.
         SODK=.True.
         IRELMP=1
         iOpt = iOr(iOpt,2**3)
         Go To 9988
*
*        one-centre no-pair operators (DK2)
*
 1007    Continue
         NoPairL=.True.
         SODK=.True.
         IRELMP=2
         iOpt = iOr(iOpt,2**3)
         Go To 9988
*
*        one-centre no-pair operators (DK3)
*
 1008    Continue
         NoPairL=.True.
         SODK=.True.
         IRELMP=3
         iOpt = iOr(iOpt,2**3)
         Go To 9988
*
*        one-centre no-pair operators (DK3)
*
 1010    Continue
         NoPairL=.True.
         SODK=.True.
         IRELMP=4
         iOpt = iOr(iOpt,2**3)
         Go To 9988
*
*        one-centre RESC operators
*
 1009    Continue
         NoPairL=.True.
         IRELMP=11
         iOpt = iOr(iOpt,2**3)
         Go To 9988
*
*        one-centre ZORA operators
*
 9001    Continue
         NoPairL=.True.
         IRELMP=21
         iOpt = iOr(iOpt,2**3)
         Go To 9988
*
*        one-centre ZORA-FP operators
*
 9002    Continue
         NoPairL=.True.
         SODK=.True.
         IRELMP=22
         iOpt = iOr(iOpt,2**3)
         Go To 9988
*
*        one-centre IORA operators
*
 9003    Continue
         NoPairL=.True.
         SODK=.True.
         IRELMP=23
         iOpt = iOr(iOpt,2**3)
         Go To 9988
*
 999     Continue
         If (iAnd(iOpt,2**1).ne.0 .and. iAnd(iOpt,2**3).ne.0) Then
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
            Call CalcAMt(iOpt,LUQRP,MPLbl,nAIMP,iMPShll+1,nProj,
     &                   iPrSh+1,ipCff,nExp,nBasis,MxShll,
     &                   DBLE(iAtmNr),DInf,nDInf)
            Close (LUQRP)
         End If
      End If
*
      lAng = Max(lAng,nProj,nAIMP)
      If (.not.inLn3) Close(lUnit)
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iBSLbl)
      End
