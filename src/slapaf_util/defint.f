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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      SubRoutine DefInt(BVct,nBVct,Labels,BMtrx,nQQ,nAtom,nLines,Value,
     &                  rInt,Lbl,Name,Coor,dMass,
     &                  jStab,nStab,mxdc,rMult,nDim,Redundant)
************************************************************************
*                                                                      *
* Object: to generate the B matrix which is the transformation matrix  *
*         between a infinitesimal displacement in the symmetry adapted *
*         internal coordinates and the symmetry unique cartesian       *
*         coordinates.                                                 *
*                                                                      *
*     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
*             University of Lund, SWEDEN                               *
*             May '91                                                  *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "Molcas.fh"
      Real*8 BVct(3*nAtom,nBVct), Value(nBVct), BMtrx(3*nAtom,nQQ),
     &       rInt(nQQ), Coor(3,nAtom), dMass(nAtom),
     &       rMult(nBVct)
      Character Labels(nBVct)*8, Type*6, Temp*120, Lbl(nQQ)*8, cNum*4,
     &          Name(nAtom)*(LENIN), Line*120, Format*8, filnam*16
      Logical First, lWrite, Flip, lPIC(6*nAtom), lAtom(nAtom), lErr,
     &        Redundant
      Integer nStab(mxdc), jStab(0:7,mxdc)
      Save First
      Data First/.True./
*
      lWrite = .False.
      lErr = .False.
      iRout = 30
      iPrint = nPrint(iRout)
      If (iPrint.ge.6) lWrite = .True.
      Do i = 1, 6*nAtom
         lPIC(i)  = .True.
      EndDo
      Do i = 1, nAtom
         lAtom(i) = .True.
      EndDo
      Do jBVct = 1, nBVct
         Labels(jBVct)=' '
      End Do
*
      Lu=6
      nTemp=Len(Temp)
      Write (Format,'(A,I3.3,A)') '(F',nTemp,'.0)'
*
      Lu_UDIC=91
      filnam='UDIC'
      Call molcas_open(Lu_UDIC,filnam)
c      Open(Lu_UDIC,File=filnam,Form='Formatted',Status='OLD')
      Rewind(Lu_UDIC)
*
      call dcopy_(nBVct,[Zero],0,rMult,1)
      If (iPrint.eq.99) First = .True.
      If (lWrite) Then
         Write (6,*)
         Write (6,'(80A)') ('*',i=1,60)
         Write (6,*) ' User-Defined Internal coordinates'
         Write (6,'(80A)') ('-',i=1,60)
         Do iLines = 1, nLines
            Read (Lu_UDIC,'(A)') Temp
            Write (6,'(A)') Temp
         End Do
         Write (6,'(80A)') ('*',i=1,60)
         Write (6,*)
         Write (6,*)
         Write (6,*) '***********************************************'//
     &               '**************'
         Write (6,*) '* Values of primitive internal coordinates     '//
     &               '             *'
         Write (6,*) '-----------------------------------------------'//
     &               '--------------'
         Rewind(Lu_UDIC)
      End If
*
*     Step 1. Set up the b vectors from which we will define the
*     internal coordinates.
*
      iBVct = 0
      Do 10 iLines = 1, nLines
         Flip=.False.
         Read (Lu_UDIC,'(A)') Line
         Temp=Line
         Call UpCase(Temp)

         If (Temp(1:4).Eq.'VARY') Go To 100
         iBVct = iBVct + 1
*
*        Move the label of the internal coordinate
*
         neq = Index(Line,'=')
         If (neq.Eq.0) Then
            Call WarningMessage(2,'Error in DefInt')
            Write (6,*)
            Write (6,'(A)') '***********************************'
            Write (6,'(A)') ' Syntax error in line :            '
            Write (6,'(A)') Line(1:33),'...'
            Write (6,'(A)') '***********************************'
            Call Quit_OnUserError()
         Else
            iFrst = 1
            Call NxtWrd(Line,iFrst,iEnd)
            jEnd = iEnd
            If (Line(iEnd:iEnd).eq.'=') jEnd = jEnd - 1
            If (jEnd-iFrst+1.gt.8) Then
               Call WarningMessage(2,'Error in DefInt')
               Write (6,*)
               Write (6,'(A)') '***********************************'
               Write (6,'(A)') ' Syntax error in line :            '
               Write (6,'(A)') Line(1:33),'...'
               Write (6,'(A,A)') Line(iFrst:jEnd),
     &               ' has more than 8 character'
               Write (6,'(A)') '***********************************'
               Call Quit_OnUserError()
            End If
            Call ChkLbl(Line(iFrst:jEnd),Labels,iBVct-1)
            Labels(iBVct) = Line(iFrst:jEnd)
         End If
*
*--------Construct the corresponding transformation vector
*
         mCntr = 0
         If (Index(Temp,'CART').Ne.0) Then
            nCntr=1
            nGo = Index(Temp,'CART')
            nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
            If (Index(Temp(nGo:nTemp),'X').Ne.0) Then
               nGo = nGo-1+Index(Temp(nGo:nTemp),'X')
               nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
               Type='X     '
            Else If (Index(Temp(nGo:nTemp),'Y').Ne.0) Then
               nGo = nGo-1+Index(Temp(nGo:nTemp),'Y')
               nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
               Type='Y     '
            Else If (Index(Temp(nGo:nTemp),'Z').Ne.0) Then
               nGo = nGo-1+Index(Temp(nGo:nTemp),'Z')
               nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
               Type='Z     '
            Else
               nGo=-1
               Call WarningMessage(2,'Error in DefInt')
               Write (6,*)
               Write (6,*) '*********** ERROR ************'
               Write (6,*) ' DefInt: wrong cartesian type '
               Write (6,'(A,A)') ' Temp=',Temp
               Write (6,*) '******************************'
               Call Quit_OnUserError()
            End If
         Else If (Index(Temp,'BOND').Ne.0) Then
            nCntr=2
            nGo = Index(Temp,'BOND')
            nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
            Type='STRTCH'
         Else If (Index(Temp,'LANGLE(2)').Ne.0) Then
            nCntr=3
            nGo = Index(Temp,'LANGLE(2)')
            nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
            Type ='LBEND2'
         Else If (Index(Temp,'LANGLE(1)').Ne.0) Then
            nCntr=3
            nGo = Index(Temp,'LANGLE(1)')
            nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
            Type ='LBEND1'
         Else If (Index(Temp,'ANGL').Ne.0) Then
            nCntr=3
            nGo = Index(Temp,'ANGL')
            nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
            Type ='BEND  '
         Else If (Index(Temp,'DIHE').Ne.0) Then
            nCntr=4
            nGo = Index(Temp,'DIHE')
            nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
            Type ='TRSN  '
            Flip=.True. .and. .Not. First
         Else If (Index(Temp,'OUTO').Ne.0) Then
            nCntr=4
            nGo = Index(Temp,'OUTO')
            nGo = nGo-1+Index(Temp(nGo:nTemp),' ')
            Type ='OUTOFP'
         Else If (Index(Temp,'DISS').Ne.0) Then
            i1 = Index(Line,'(')
            i2 = Index(Line,'+')
            i3 = Index(Line,')')
            If (i1.ge.i2 .or. i2.ge.i3) Then
               Call WarningMessage(2,'Error in DefInt')
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
         Else
            nGo=-1
            Call WarningMessage(2,'Error in DefInt')
            Write (6,*)
            Write (6,*) '*********** ERROR ***********'
            Write (6,*)' Line contains syntax error !'
            Write (6,'(A)') Line
            Write (6,*) '*****************************'
            Call Quit_OnUserError()
         End If
*
         msAtom = nCntr + mCntr
         Call GetMem('xyz ','Allo','Real',ipxyz ,3*msAtom)
         Call GetMem('Temp','Allo','Real',ipTemp,3*msAtom)
         Call GetMem('Ind ','Allo','Inte',ipInd ,2*msAtom)
         Call GetMem('Mass','Allo','Real',ipMass,2*msAtom)
         Call GetMem('TMtrx','Allo','Real',ipTM,9*nAtom*(nCntr+mCntr))
*
         Call Cllct(Line(nGo:nTemp),BVct(1,iBVct),Value_Temp,
     &              Name,nAtom,Coor,nCntr,mCntr,
     &              Work(ipxyz),Work(ipTemp),iWork(ipInd),Type,
     &              dMass,Work(ipMass),Work(ipTM),lWrite,
     &              Labels(iBVct),lWrite,jStab,nStab,mxdc,
     &              rMult(iBVct),lAtom)
*
         If (.Not.First .and.
     &       Type.eq.'TRSN' .and.
     &       Abs(Value_Temp).lt.Pi*Half) Flip=.False.
         If (Flip .and.
     &       Value(iBVct)*Value_Temp.lt.Zero ) Then
*           Write (Lu,*) 'Flip Sign for ',Labels(iBVct)
            If (Value(iBVct).lt.Zero) Then
               Value(iBVct)= -Pi - (Pi-Value_Temp)
            Else
               Value(iBVct)=  Pi + (Pi+Value_Temp)
            End If

         Else
            Value(iBVct)= Value_Temp
         End If

*
         Call GetMem('TMtrx','Free','Real',ipTM,9*nAtom*(nCntr+mCntr))
         Call GetMem('Mass','Free','Real',ipMass,2*msAtom)
         Call GetMem('Ind ','Free','Inte',ipInd ,2*msAtom)
         Call GetMem('Temp','Free','Real',ipTemp,3*msAtom)
         Call GetMem('xyz ','Free','Real',ipxyz ,3*msAtom)
*
 10   Continue
      Call WarningMessage(2,'Error in DefInt')
      Write (6,*) '**********************************************'
      Write (6,*) ' ERROR: No internal coordinates are defined ! '
      Write (6,*) '**********************************************'
      Call Quit_OnUserError()
*
 100  Continue
      nDefPICs = iBVct
      If (iPrint.ge.59) Call RecPrt(' The B-vectors',' ',
     &                              BVct,3*nAtom,nBVct)
      If (iPrint.ge.19) Call RecPrt(
     &        ' Values of primitive internal coordinates / au or rad',
     &                             ' ',Value,nBVct,1)
*
*     Step 2. Define internal coordinates as linear combinations of
*     the previously defined primitive internal coordinates.
*
      iBMtrx = 0
      jLines=iLines
 201  Continue
         jLines=jLines+1
         If (jLines.gt.nLines) Go To 200
*
*------- Jump over the FIX keyword if present
*
         Read (Lu_UDIC,'(A)') Line
         Temp = Line
         Call UpCase(Temp)
         If (Temp(1:3).Eq.'FIX') Go To 20
         If (Temp(1:4).Eq.'ROWH') Go To 30
*
         iBMtrx = iBMtrx + 1
         rInt(iBMtrx) = Zero
         RR=Zero
*
         iFrst = 1
         Call NxtWrd(Line,iFrst,iEnd)
         jEnd = iEnd
         If (Line(iEnd:iEnd).eq.'=') jEnd=iEnd-1
         Lbl(iBMtrx) = Line(iFrst:jEnd)
         neq = Index(Line,'=')
         If (neq.Eq.0) Then
*
*           a single vector (this will only extend over one line)
*
            iBVct = 0
            Do 21 jBVct = 1, nBVct
               If (Line(iFrst:jEnd).eq.Labels(jBVct))
     &            iBVct=jBVct
 21         Continue
            If (iBVct.eq.0) Then
               Call WarningMessage(2,'Error in DefInt')
               Write (6,*)
               Write (6,*) '*******************************'
               Write (6,*) ' ERROR: A single vector        '
               Write (6,*) ' Undefined internal coordinate '
               Write (6,'(A,A)') Line
               Write (6,'(A,A)') Line(iFrst:jEnd)
               Write (6,*) '*******************************'
               Call ErrTra
               Call Quit_OnUserError()
            End If
*
            lPIC(iBVct) = .False.
            call dcopy_(3*nAtom,BVct(1,iBVct),1,BMtrx(1,iBMtrx),1)
            Call DScal_(3*nAtom,rMult(iBVct)**2,
     &                 BMtrx(1,iBMtrx),1)
            rInt(iBMtrx) = rMult(iBVct)**2*Value(iBVct)
            RR=RR+rMult(iBVct)**2
*
         Else
*
*-----------A linear combination of vectors
*
            call dcopy_(3*nAtom,[Zero],0,BMtrx(1,iBMtrx),1)
            iFrst = neq + 1
            Sgn=One
*
*-----------Process the rest of the line and possible extension lines
*
 22         Continue
*              Get the factor
               Call NxtWrd(Line,iFrst,iEnd)
               Temp=Line(iFrst:iEnd)
               Read (Temp,Format) Fact
               Fact = Fact * Sgn
               iFrst = iEnd + 1
*              Get the label
               Call NxtWrd(Line,iFrst,iEnd)
               iBVct = 0
               Do jBVct = 1, nBVct
                  If (Line(iFrst:iEnd).eq.Labels(jBVct)) iBVct=jBVct
               End Do
               If (iBVct.eq.0) Then
                  Call WarningMessage(2,'Error in DefInt')
                  Write (6,*)
                  Write (6,*) '************ ERROR *************'
                  Write (6,*) ' Linear combinations of vectors '
                  Write (6,*) ' Undefined internal coordinate  '
                  Write (6,'(A,A)') Line
                  Write (6,'(A,A)') Line(iFrst:iEnd)
                  Write (6,*) '********************************'
                  Call ErrTra
                  Call Quit_OnUserError()
               End If
*
               lPIC(iBVct) = .False.
               Call DaXpY_(3*nAtom,Fact*rMult(iBVct)**2,
     &                    BVct(1,iBvct),1,
     &                    BMtrx(1,iBMtrx),1)
               rInt(iBMtrx) = rInt(iBMtrx)
     &                      + Fact * rMult(iBVct)**2
     &                      * Value(iBVct)
               RR=RR+rMult(iBVct)**2*Fact**2
*
               iFrst = iEnd + 1
               Temp=Line(iFrst:nTemp)
               nPlus = Index(Temp,'+')
               nMinus= Index(Temp,'-')
            If (nPlus.ne.0.and.(nPlus.lt.nMinus.eqv.nMinus.gt.0)) Then
               Sgn=One
               iFrst = iFrst + nPlus
               Go To 22
            End If
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
                  Call WarningMessage(2,'Error in DefInt')
                  Write(6,*)
                  Write(6,*) '********** ERROR *********'
                  Write(6,*) ' DefInt: jLines.gt.nLines '
                  Write(6,*) '**************************'
                  Call Quit_OnUserError()
               End If
               Read (Lu_UDIC,'(A)') Line
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
               Else
                  Call WarningMessage(2,'Error in DefInt')
                  Write (6,*)
                  Write (6,*) '************** ERROR *************'
                  Write (6,*) ' Syntax Error: first character in '
     &                      //' extension line is not + or -'
                  Write (6,'(A)') Line
                  Write (6,'(3A)') '-->',Line(iFrst:iEnd),'<--'
                  Write (6,*) '**********************************'
                  Call Quit_OnUserError()
               End If
            End If
*
         End If
*
       rInt(iBMtrx)=rInt(iBMtrx)/Sqrt(RR)
       Call DScal_(3*nAtom,One/Sqrt(RR),BMtrx(1,iBMtrx),1)
*
 20   Continue
      Go To 201
*
* --- Skip the  RowH  section of input ---
*
 30   jLines=jLines+1
      If (jLines.gt.nLines) Go To 200
      Read (Lu_UDIC,'(A)') Line
      Temp = Line
      Call UpCase(Temp)
      Go To 30
*-----------------------------------------
 200  Continue
      If (iPrint.ge.99) Call RecPrt(' The B-matrix',' ',
     &                              BMtrx,3*nAtom,nQQ)
*
      If (lWrite) Then
         Write (6,*)
         Write (6,*)
         Write (6,*) '*********************************************'
         Write (6,*) '* Value of internal coordinates / au or rad *'
         Write (6,*) '---------------------------------------------'
         Write (6,'(1X,A,2X,F10.4)')
     &         (Lbl(iInt),rInt(iInt),iInt=1,nQQ)
         Write (6,*)
         Call XFlush(6)
      End If
*
* --- Some checks: Warnings & Errors ---
*
      iDiff = 3*nAtom-nDim
      If (iDiff.eq.0) Then
        cNum = '3N'
      Else
        Write(cNum,'(A,I1)') '3N-',iDiff
      End If

      If (iBMtrx.LT.nDim) then
         Write (6,*) '**************** ERROR **********************'
         Write (6,*) ' N.r of Internal Coordinates lower than ',cNum
         Write (6,*) '*********************************************'
         Write (6,*)
         lErr = .True.
      EndIf

      If (iBMtrx.GT.nDim.and..NOT.Redundant) then
         Write (6,*) '***************** ERROR ***********************'
         Write (6,*) ' N.r of Internal Coordinates greater than ',cNum
         Write (6,*) '***********************************************'
         Write (6,*)
         lErr = .True.
      EndIf

      If (nDefPICs.LT.iBMtrx) then
         Write (6,*)
         Write (6,*) '****************** ERROR ******************'
         Write (6,*) ' Too many Internal Coordinates !           '
         Write (6,'(A,I4)') ' N.r of Primitive Internal Coordinates:',
     &                                                    nDefPICs
         Write (6,'(A,I4)') ' N.r of Internal Coordinates:          ',
     &                                                      iBMtrx
         Write (6,*) '*******************************************'
         Write (6,*)
         lErr = .True.
      EndIf

      Do i = 1, nDefPICs
         If (lPIC(i)) then
            Write (6,*)
            Write (6,*) '*************** ERROR *******************'
            Write (6,*) ' Primitive Internal Coordinate not used: '
            Write (6,*) ' ',Labels(i)
            Write (6,*) '*****************************************'
            Write (6,*)
            lErr = .True.
         EndIf
      EndDo

      Do i = 1, nAtom
         If (lAtom(i)) then
            Call WarningMessage(2,'Error in DefInt')
            Write (6,*)
            Write (6,*) '*********** ERROR ****************'
            Write (6,*) ' No Coordinate defined for atom:  '
            Write (6,*) ' ',Name(i)
            Write (6,*) '**********************************'
            Write (6,*)
            lErr = .True.
         EndIf
      End Do

      If (lErr) Call Quit_OnUserError()

      If (iBMtrx.ne.nQQ) Then
         Call WarningMessage(2,'Error in DefInt')
         Write (6,*)
         Write (6,*) '******************* ERROR *****************'
         If (iBMtrx.gt.nQQ)
     &      Write (6,*) ' Too many internal coordinates are defined'
         If (iBMtrx.lt.nQQ)
     &      Write (6,*) ' Too few internal coordinates are defined'
         Write (6,*) ' You have defined', iBMtrx
         Write (6,*) ' There should be ', nQQ
         Write (6,*) '*******************************************'
         Call Quit_OnUserError()
      End If

      Close(Lu_UDIC)
      First = .False.
      Return
      End
