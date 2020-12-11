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
* Copyright (C) 1991,1997, Roland Lindh                                *
************************************************************************
      SubRoutine DefInt2(BVct,dBVct,nBvct,Labels,BMtrx,mInt,nAtom,
     &                   nLines,Value,rInt,rInt0,Lbl,
     &                   lWrite,rMult,dBMtrx,Value0,lIter,iFlip)
************************************************************************
*                                                                      *
* Object: to generate the B matrix for the constraints                 *
*                                                                      *
*     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
*             University of Lund, SWEDEN                               *
*             May '91                                                  *
*                                                                      *
*             Modified to constraints, June '97 by R. Lindh            *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "Molcas.fh"
      Real*8 BVct(3*nAtom,nBVct), dBVct(3*nAtom,3*nAtom,nBVct),
     &       Value(nBVct), BMtrx(3*nAtom,mInt), rInt(mInt), rInt0(mInt),
     &       rMult(nBVct,nBVct),
     &       dBMtrx(3*nAtom,3*nAtom,mInt), Value0(nBVct), MaxErr
      Character Line*120, Labels(nBVct)*8, Type*6, Format*8,
     &          Temp*120, Lbl(mInt)*8, filnam*16
      Logical lWrite, Start, rInt0_on_file,rInt0_in_memory, InSlapaf
      Integer, Parameter:: Flip=1, NoFlip=0
      Integer, External:: StrnLn
      Integer iFlip(nBVct)
      Character(LEN=100), External:: Get_SuperName
      Integer, Allocatable:: Ind(:), tpc(:)
      Real*8, Allocatable:: Hess(:), Mass(:), Grad(:), xyz(:), r0(:)
#include "angstr.fh"

      Lu=6
*
*     Initiate some stuff for automatic setting of
*     the constraints to be those that the structure
*     actually has initially.
*
      rInt0_on_file=.FALSE.
      InSlapaf = (Get_SuperName().eq.'slapaf')
      If (InSlapaf) Call qpg_dArray('rInt0',rInt0_on_file,nrInt0)
      If (.not.rInt0_on_File) nrInt0 = mInt
      rInt0_in_memory=.FALSE.
*
      iRout = 30
      iPrint = nPrint(iRout)
      Start=lIter.eq.1
      Call ICopy(nBVct,[Flip],0,iFlip,1)
*
      nTemp=Len(Temp)
      Write (Format,'(A,I3.3,A)') '(F',nTemp,'.0)'
*
      Lu_UDC=91
      filnam='UDC'
      Call molcas_open(Lu_UDC,filnam)
c      Open(Lu_UDC,File=filnam,Form='FORMATTED',Status='OLD')
      Rewind (Lu_UDC)
*
      call dcopy_(nBVct**2,[Zero],0,rMult,1)
      If (iPrint.ge.99) lWrite = .True.
      If (iPrint.ge.99 .or. lWrite) Then
          Write (Lu,*)
          Call CollapseOutput(1,'Constraints section')
          Write (Lu,'(34X,A)') 'CONSTRAINTS'
         Write (Lu,*)
         Write (Lu,'(80A)') ('*',i=1,80)
         Do iLines = 1, nLines
            Read(Lu_UDC,'(A)') Line
            Write (Lu,'(A)') Line(:StrnLn(Line))
         End Do
         Write (Lu,'(80A)') ('*',i=1,80)
         Write (Lu,*)
         Write (Lu,*)
         Write (Lu,*) '*********************************************'//
     &                '****************'
         Write (Lu,*) '* Values of the primitive constraints        '//
     &                '               *'
         Write (Lu,*) '*********************************************'//
     &                '****************'
         Rewind (Lu_UDC)
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Step 1. Set up the b vectors from which we will define the
*             constraints.
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      iBVct = 0
      Call mma_allocate(tpc ,nBVct,Label='tpc')
      Do 10 iLines = 1, nLines
         Read(Lu_UDC,'(A)') Line
         Temp = Line
         Call UpCase(Temp)
         If (Temp(1:4).Eq.'VALU') Go To 100
         iBVct = iBVct + 1
*
*        Move the label of the internal coordinate
*
         neq = Index(Line,'=')
         If (neq.Eq.0) Then
            Call WarningMessage(2,'Error in DefInt2')
            Write (Lu,'(A)') ' Syntax error in line:'
            Write (Lu,'(A)') Line
            Call Quit_OnUserError()
         Else
            iFrst = 1
            Call NxtWrd(Line,iFrst,iEnd)
            jEnd = iEnd
            If (Line(iEnd:iEnd).eq.'=') jEnd = jEnd - 1
            If (jEnd-iFrst+1.gt.8) Then
               Call WarningMessage(2,'Error in DefInt2')
               Write (Lu,'(A,A)') Line(iFrst:jEnd),
     &               ' has more than 8 character, syntax error!'
               Call Quit_OnUserError()
            End If
            Call ChkLbl(Line(iFrst:jEnd),Labels,iBVct-1)
            Labels(iBVct) = Line(iFrst:jEnd)
         End If
*
*--------Construct the corresponding transformation vector
*
         mCntr = 0
         iType = 0
         If (Index(Temp,'CART').Ne.0) Then
            nCntr=1
            nGo = Index(Temp,'CART')
            nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
            If (Index(Temp(nGo:nTemp),'X').Ne.0) Then
               nGo = nGo-1+Index(Temp(nGo:nTemp),'X')
               nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
               Type='X     '
               iType=1
            Else If (Index(Temp(nGo:nTemp),'Y').Ne.0) Then
               nGo = nGo-1+Index(Temp(nGo:nTemp),'Y')
               nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
               Type='Y     '
               iType=2
            Else If (Index(Temp(nGo:nTemp),'Z').Ne.0) Then
               nGo = nGo-1+Index(Temp(nGo:nTemp),'Z')
               nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
               Type='Z     '
               iType=3
            Else
              nGo=-1
              Call WarningMessage(2,'Error in DefInt2')
              Write (Lu,*) 'DefInt2: wrong cartesian type'
              Write (Lu,'(A,A)') 'Temp=',Temp
              Call Quit_OnUserError()
            End If
         Else If (Index(Temp,'BOND').Ne.0) Then
            nCntr=2
            nGo = Index(Temp,'BOND')
            nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
            Type='STRTCH'
            iType=4
         Else If (Index(Temp,'LANGLE(2)').Ne.0) Then
            nCntr=3
            nGo = Index(Temp,'LANGLE(2)')
            nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
            Type ='LBEND2'
            iType=5
         Else If (Index(Temp,'LANGLE(1)').Ne.0) Then
            nCntr=3
            nGo = Index(Temp,'LANGLE(1)')
            nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
            Type ='LBEND1'
            iType=6
         Else If (Index(Temp,'ANGL').Ne.0) Then
            nCntr=3
            nGo = Index(Temp,'ANGL')
            nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
            Type ='BEND  '
            iType=7
         Else If (Index(Temp,'DIHE').Ne.0) Then
            nCntr=4
            nGo = Index(Temp,'DIHE')
            nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
            Type ='TRSN  '
* AOM!! Remove flip for Torsions!!!
            iFlip(iBVct)=NoFlip
            iType=8
         Else If (Index(Temp,'OUTO').Ne.0) Then
            nCntr=4
            nGo = Index(Temp,'OUTO')
            nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
            Type ='OUTOFP'
            iFlip(iBVct)=NoFlip
            iType=9
         Else If (Index(Temp,'EDIF').Ne.0) Then
            nCntr=nAtom
            nGo = Index(Temp,'EDIF')
            nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
            Type ='EDIFF '
            iFlip(iBVct)=NoFlip
            iType=10
         Else If (Index(Temp,'SPHE').Ne.0) Then
            nCntr=nAtom
            nGo = Index(Temp,'SPHE')
            nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
            Type ='SPHERE'
            iType=11
         Else If (Index(Temp,'NAC ').Ne.0) Then
            nCntr=nAtom
            nGo = Index(Temp,'NAC ')
            nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
            Type ='NAC   '
            iType=12
         Else If (Index(Temp,'TRAN').Ne.0) Then
            nCntr=nAtom
            nGo = Index(Temp,'TRAN')
            nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
            Type ='TRANSV'
            iFlip(iBVct)=NoFlip
            iType=13
         Else If (Index(Temp,'DISS').Ne.0) Then
            i1 = Index(Line,'(')
            i2 = Index(Line,'+')
            i3 = Index(Line,')')
            If (i1.ge.i2 .or. i2.ge.i3) Then
               Call WarningMessage(2,'Error in DefInt2')
               Write (6,*)
               Write (6,*) '********** ERROR ************'
               Write (6,*)' Line contains syntax error !'
               Write (6,'(A)') Line
               Write (6,*) i1,i2,i3
               Write (6,*) '*****************************'
               Call Quit_OnUserError()
            End If
            nGo = i3+1
            Temp = Line(i1+1:i2-1)
            Read (Temp,Format) Tmp
            nCntr=NInt(Tmp)
            Temp=Line(i2+1:i3-1)
            Read (Temp,Format) Tmp
            mCntr=NInt(Tmp)
            Type ='DISSOC'
            iType=14
         Else
            nGo=-1
            Call WarningMessage(2,'Error in DefInt2')
            Write (Lu,*)' Line contains syntax error!'
            Write (Lu,'(A)') Line
            Call Quit_OnUserError()
         End If
         tpc(iBVct)=iType
*
         msAtom=nCntr+mCntr
         Call mma_allocate(xyz ,3*msAtom,Label='xyz')
         Call mma_allocate(Grad,3*msAtom,Label='Grad')
         Call mma_allocate(Ind ,2*msAtom,Label='Ind')
         Call mma_allocate(Mass,msAtom,Label='Mass')
         Call mma_allocate(Hess,(3*msAtom)**2,Label='Hess')
*
         Call Cllct2(Line(nGo:nTemp),BVct(1,iBVct),dBVct(1,1,iBVct),
     &               Value(iBVct),nAtom,nCntr,mCntr,xyz,Grad,
     &               Ind,Type,Mass,Labels(iBVct),lWrite,
     &               rMult(iBVct,iBVct),Hess,lIter)
*
         If (Type.eq.'TRSN  ' .and.
     &       Abs(Value(iBVct)).lt.Pi*Half) iFlip(iBVct)=NoFlip
*
         Call mma_deallocate(Hess)
         Call mma_deallocate(Mass)
         Call mma_deallocate(Ind)
         Call mma_deallocate(Grad)
         Call mma_deallocate(xyz)
*
 10   Continue
      Call WarningMessage(2,'Error in DefInt2')
      Write (Lu,*) 'DefInt2: No internal coordinates are defined!'
      Call Quit_OnUserError()
*
 100  Continue

*
*---- Process primitive value to correct for flips
*
      If (.Not.Start) Then
         Do iBVct = 1, nBVct
*
*---------- Test if we have a flip in the sign of the value
*
            If (Value(iBVct)*Value0(iBVct).lt.Zero .and.
     &          iFlip(iBVct).eq.Flip) Then
C              Write (Lu,*) 'Flip Sign for ',Labels(iBVct)
               Value(iBVct)=-Value(iBVct)
            End If
         End Do
      End If
      call dcopy_(nBVct,Value,1,Value0,1)
*

      If (iPrint.ge.59) Call RecPrt(' The B-vectors',' ',
     &                              BVct,3*nAtom,nBVct)
      If (iPrint.ge.19) Then
         Call RecPrt(
     &        ' Value of primitive internal coordinates / au or rad',
     &                             ' ',Value,nBVct,1)
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Step 2. Define constraints as linear combinations of
*             the previously defined primitive internal coordinates.
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      iBMtrx = 0
      jLines=iLines
 201  Continue
         jLines=jLines+1
         If (jLines.gt.nLines) Go To 200
*
         Read (Lu_UDC,'(A)') Line
         Temp = Line
         Call UpCase(Temp)
*
         iBMtrx = iBMtrx + 1
         rInt(iBMtrx) = Zero
         Write(Lbl(iBMtrx),'(A,I3.3)') 'Cns',iBMtrx
         RR=Zero
*
*------- Find the label of the first primitive
*
         iFrst = 1
         Call NxtWrd(Line,iFrst,iEnd)
         jEnd = iEnd
         If (Line(iEnd:iEnd).eq.'=') jEnd=iEnd-1
*
         nPlus = Index(Line,' + ')
         nMinus= Index(Line,' - ')
*                                                                      *
************************************************************************
*                                                                      *
         If (nPlus.Eq.0.and.nMinus.Eq.0.and.Index(Line,'&').eq.0) Then
*                                                                      *
************************************************************************
*                                                                      *
*           a single vector (this will only extend over one line)
*
            If (Index(Line,'&').ne.0) Then
               Call WarningMessage(2,'Error in DefInt2')
               Write (Lu,*) 'Single vector lines should not extend'
               Write (Lu,*) 'over more than one line!'
               Write (Lu,'(A)') Line
               Call Quit_OnUserError()
            End If
            iBVct = 0
            Do jBVct = 1, nBVct
               If (Line(iFrst:jEnd).eq.Labels(jBVct))
     &            iBVct=jBVct
            End Do
            If (iBVct.eq.0) Then
               Call WarningMessage(2,'Error in DefInt2')
               Write (Lu,*) ' A single vector'
               Write (Lu,*) ' Undefined internal coordinate'
               Write (Lu,'(A,A)') Line
               Write (Lu,'(A,A)') Line(iFrst:jEnd)
               Call ErrTra
               Call Quit_OnUserError()
            End If
*
            call dcopy_(3*nAtom,BVct(1,iBVct),1,BMtrx(1,iBMtrx),1)
            call dcopy_((3*nAtom)**2,dBVct(1,1,iBVct),1,
     &                              dBMtrx(1,1,iBMtrx),1)
            rInt(iBMtrx) = Value(iBVct)
*
            iFrst=iEnd+1
            Call NxtWrd(Line,iFrst,iEnd)
            If (Line(iEnd:iEnd).eq.'=') Then
               iFrst=iEnd+1
               Call NxtWrd(Line,iFrst,iEnd)
            End If
            Temp=Line(iFrst:iEnd)
*
            If (Index(Temp,'FIX').ne.0) Then
*
*              Pick up values from the runfile. Written there on
*              the first iteration.
*
               If (.Not.rInt0_in_memory) Then
                  rInt0_in_memory=.TRUE.
                  Call mma_allocate(r0,nrInt0,Label='r0')
                  If (rInt0_on_file) Then
                     Call Get_dArray('rInt0',r0,nrInt0)
                  Else
                    r0(:)=Zero
                  End If
               End If
*
               If (rInt0_on_file) Then
                  rInt0(iBMtrx)=r0(iBMtrx)
               Else
                  r0(iBmtrx)   =rInt(iBMtrx)
                  rInt0(iBMtrx)=rInt(iBMtrx)
               End If
*
            Else
*
*              Read value from input file.
*
               Read (Temp,Format) rInt0(iBMtrx)
               Temp = Line
               Call UpCase(Temp)
               If (Index(Temp,'ANGSTROM').ne.0) rInt0(iBMtrx) =
     &                                       rInt0(iBMtrx)/angstr
               If (Index(Temp,'DEGREE').ne.0) rInt0(iBMtrx) =
     &                                     rInt0(iBMtrx)*Pi/1.800D+02
            End If
*
* AOM: trying to correct torsion...
            If (tpc(iBVct).eq.8) Then
              n0 = Int((rInt(iBMtrx)-rInt0(iBMtrx))/(Two*Pi))
              If (Abs(rInt(IBMtrx)-rInt0(iBMtrx)-Dble(n0)*Two*Pi) .gt.
     &           Abs(rInt(IBMtrx)-rInt0(iBMtrx)-Dble(n0+1)*Two*Pi)) Then
                n0=n0+1
              Else If (Abs(rInt(IBMtrx)-rInt0(iBMtrx)-Dble(n0)*Two*Pi)
     &                 .gt.
     &           Abs(rInt(IBMtrx)-rInt0(iBMtrx)-Dble(n0-1)*Two*Pi)) Then
                n0=n0-1
              Endif
              rInt(iBMtrx)=rInt(iBMtrx)-Dble(n0)*Two*Pi
            Endif
*                                                                      *
************************************************************************
*                                                                      *
         Else
*                                                                      *
************************************************************************
*                                                                      *
*
*-----------A linear combination of vectors
*
            call dcopy_(3*nAtom,[Zero],0,BMtrx(1,iBMtrx),1)
            call dcopy_((3*nAtom)**2,[Zero],0,dBMtrx(1,1,iBMtrx),1)
            iFrst = 1
            Sgn=One
*
*-----------Process the rest of the line and possible extension lines
*
 22         Continue
*
*            ***********************************************************
*            *                                                         *
*------------> Get the factor
               Call NxtWrd(Line,iFrst,iEnd)
               Temp=Line(iFrst:iEnd)
               Read (Temp,Format) Fact
               Fact = Fact * Sgn
               iFrst = iEnd + 1
*------------> Get the label
               Call NxtWrd(Line,iFrst,iEnd)
               If (Line(iEnd:iEnd).eq.'=') iEnd=iEnd-1
               iBVct = 0
               Do jBVct = 1, nBVct
                  If (Line(iFrst:iEnd).eq.Labels(jBVct))
     &               iBVct=jBVct
               End Do
               If (iBVct.eq.0) Then
                  Call WarningMessage(2,'Error in DefInt2')
                  Write (Lu,*) ' Linear combinations of vectors'
                  Write (Lu,*) ' Undefined internal coordinate'
                  Write (Lu,'(A,A)') Line
                  Write (Lu,'(A,A)') Line(iFrst:iEnd)
                  Call Quit_OnUserError()
               End If
*
               Call DaXpY_(3*nAtom,Fact,BVct(1,iBvct),1,
     &                                 BMtrx(1,iBMtrx),1)
               Call DaXpY_((3*nAtom)**2,Fact,dBVct(1,1,iBvct),1,
     &                                      dBMtrx(1,1,iBMtrx),1)
               rInt(iBMtrx) = rInt(iBMtrx)+ Fact * Value(iBVct)
*            *                                                         *
*            ***********************************************************
*
               iFrst = iEnd + 1
 25            Continue
               Temp=Line(iFrst:nTemp)
               nEq   = Index(Temp,'=')
               nPlus = Index(Temp,'+')
               nMinus= Index(Temp,'-')
*
            If (nEq.ne.0 .and.
     &          (nEq.lt.nMinus.eqv.nMinus.gt.0) .and.
     &          (nEq.lt.nPlus .eqv.nPlus .gt.0)) Then
               If (Index(Line,'&').ne.0) Then
                  Call WarningMessage(2,'Error in DefInt2')
                  Write (Lu,*) 'This line should not be extended'
                  Write (Lu,'(A)') Line
                  Call Quit_OnUserError()
               End If
               iFrst = iFrst + nEq + 1
               Call NxtWrd(Line,iFrst,iEnd)
               Temp=Line(iFrst:iEnd)
*
               If (Index(Temp,'FIX').ne.0) Then
*
*                 Pick up values from the runfile. Written there on
*                 the first iteration.
*
                  If (.Not.rInt0_in_memory) Then
                     rInt0_in_memory=.TRUE.
                     Call mma_allocate(r0,nrInt0,Label='r0')
                     If (rInt0_on_file) Then
                        Call Get_dArray('rInt0',r0,nrInt0)
                     Else
                        r0(:)=Zero
                     End If
                  End If
*
                  If (rInt0_on_file) Then
                     rInt0(iBMtrx)=r0(iBMtrx)
                  Else
                     r0(iBmtrx)=rInt(iBMtrx)
                  End If
*
               Else
*
*                 Read value from input file.
*
                  Read (Temp,Format) rInt0(iBMtrx)
                  Temp = Line
                  Call UpCase(Temp)
                  If (Index(Temp,'ANGSTROM').ne.0) rInt0(iBMtrx) =
     &                                          rInt0(iBMtrx)/angstr
                  If (Index(Temp,'DEGREE').ne.0) rInt0(iBMtrx) =
     &                                        rInt0(iBMtrx)*Pi/1.800D+02
               End If
               Go To 24
            End If
*
            If (nPlus.ne.0.and.(nPlus.lt.nMinus.eqv.nMinus.gt.0)) Then
               Sgn=One
               iFrst = iFrst + nPlus
               Go To 22
            End If
*
            If (nMinus.ne.0.and.(nMinus.lt.nPlus.eqv.nPlus.gt.0)) Then
               Sgn=-One
               iFrst = iFrst + nMinus
               Go To 22
            End If
*
*---------- Here if all statements processed of this line
*
            If (Index(Line,'&').ne.0) Then
               jLines=jLines+1
               If (jLines.gt.nLines) Then
                  Call WarningMessage(2,'Error in DefInt2')
                  Write(Lu,*)'DefInt2: jLines.gt.nLines'
                  Call Quit_OnUserError()
               End If
               Read (Lu_UDC,'(A)') Line
               iFrst = 1
               Call NxtWrd(Line,iFrst,iEnd)
               If (Line(iFrst:iEnd).eq.'+') Then
                  iFrst = iEnd + 1
                  Sgn=One
                  Go To 22
               Else If (Line(iFrst:iEnd).eq.'-') Then
                  iFrst = iEnd + 1
                  Sgn=-One
                  Go To 22
               Else If (Line(iFrst:iEnd).eq.'=') Then
                  iFrst = 1
                  Go To 25
               Else
                  Call WarningMessage(2,'Error in DefInt2')
                  Write (Lu,*) ' Syntax Error: first character in '
     &                      //' extension line is not + or -'
                  Write (Lu,'(A)') Line
                  Write (Lu,'(3A)') '-->',Line(iFrst:iEnd),'<--'
                  Call Quit_OnUserError()
               End If
            End If
*
*---------- At the end of this line
*
 24         Continue
*                                                                      *
************************************************************************
*                                                                      *
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Go To 201
*
 200  Continue
      If (iPrint.ge.99) Then
         Call RecPrt(' The B-matrix',' ',BMtrx,3*nAtom,mInt)
         Do iInt = 1, mInt
            Call RecPrt(' The dB-matrix',' ',dBMtrx(1,1,iInt),3*nAtom,
     &                                                        3*nAtom)
         End Do
      End If
      Close(Lu_UDC)
      Call mma_deallocate(tpc)
      If (rint0_in_memory) Then
         If (InSlapaf) Call Put_dArray('rInt0',r0,mInt)
         Call mma_deallocate(r0)
      End If
*
*     Compute the maximum error in the constraints
*
      MaxErr = Zero
      Do iInt = 1, mInt
         MaxErr = Max(MaxErr, Abs(rInt(iInt)-rInt0(iInt)))
      End Do
      Call Put_dScalar('Max error',MaxErr)
*
      If (iPrint.ge.99 .or. lWrite) Then
         Write (Lu,*)
         Write (Lu,*)
         Write (Lu,*) '*******************************************'
         Write (Lu,*) '* Values of the constraints   / au or rad *'
         Write (Lu,*) '*******************************************'
         Write (Lu,*) '  Label        C         C0'
         Write (Lu,'(1X,A,2X,F10.6,F10.6)')
     &         (Lbl(iInt),rInt(iInt),rInt0(iInt),iInt=1,mInt)
         Write (Lu,*)
         Call CollapseOutput(0,'Constraints section')
         Write (Lu,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
