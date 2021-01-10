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
      Subroutine Fix_UDC(iRow_c,nLambda,nsAtom,nStab,Remove)
      use Slapaf_Info, only: AtomLbl
      Implicit Real*8 (A-H,O-Z)
************************************************************************
*                                                                      *
*     Object: Replace the fragment specifications with actual          *
*             constraint specifications.                               *
*             Remove constraints to be ignored.                        *
*                                                                      *
************************************************************************
#include "Molcas.fh"
      External Get_Ln
      Character(LEN=180) Get_Ln, Line, Line2, Lines(iRow_c)
      Character(LEN=16) FilNam
      Character(LEN=180) Label1, Label2, Label3, Label4,
     &                   FragLabels(iRow_c),
     &                   SoftLabels(iRow_c)
      Integer nStab(nsAtom), FragZMat(iRow_c,2),ZMatOffset
      Logical Ignore, Values, Remove
      Character(LEN=100), External:: Get_SuperName
      Character(LEN=100) SuperName
*                                                                      *
************************************************************************
*                                                                      *
      If (iRow_c.le.1) Return
      iRow_c_Save=iRow_c
      nLambda_save=nLambda
      SuperName=Get_SuperName()
*                                                                      *
************************************************************************
*                                                                      *
*     Open files for user defined internal constraints.
*
      Lu_UDC=91
      FilNam='UDC'
      Lu_UDC=IsFreeUnit(Lu_UDC)
      Call Molcas_Open(Lu_UDC,FilNam)
      ReWind(Lu_UDC)
*
      Lu_TMP=Lu_UDC+1
      FilNam='UDCTMP'
      Lu_TMP=IsFreeUnit(Lu_TMP)
      Call Molcas_Open(Lu_TMP,FilNam)
      ReWind(Lu_TMP)
*                                                                      *
************************************************************************
*                                                                      *
*     Process the content of LU_UDC and put it onto LU_TMP
*
      iRow_c=0 ! Number of lines in the section
      nZMat=0  ! Number of additional lines with constraints
      nFrag=0  ! Number of fragment constraints
      nSoft=0  ! Number of default soft constraints
      Values=.False.
*
      iZMat=0
      ZMatOffset=0
 100  Continue
         Line=Get_Ln(LU_UDC)   ! Read the line
         Call UpCase(Line)
*
         If (Line(1:4).eq.'VALU') Values=.True.
         If (Line(1:4).ne.'END') Then
*
*           Read all possible continuation lines
*
            nLines=1
            Lines(nLines)=Line
            Do While (Index(Lines(nLines),'&').gt.0)
               nLines=nLines+1
               Lines(nLines)=Get_Ln(LU_UDC)
               Call UpCase(Lines(nLines))
            End Do
            iLine=1
            Line=Lines(iLine)
            iEq=Index(Line,'=')
*
*           If there are already some ZMAT constraints in the file,
*           start numbering from the next one
*
            If (Line(1:4) .eq. 'ZMAT') Then
               Read(Line(5:7),*) Num
               ZMatOffset=Max(ZMatOffset,Num)
            End If
*
            If (.Not.Values) Then
*
*           Find out if the primitive defaults to soft
*
            iFrst=iEq+1
            Call NxtWrd(Line,iFrst,iEnd)
            If ((Line(iFrst:iFrst+3).eq.'SPHE') .or.
     &          (Line(iFrst:iFrst+3).eq.'TRAN') .or.
     &          (Line(iFrst:iFrst+3).eq.'EDIF') .or.
     &          (Line(iFrst:iFrst+3).eq.'NAC ')) Then
               nSoft=nSoft+1
               iFrst=1
               Call NxtWrd(Line,iFrst,iEnd)
               SoftLabels(nSoft)=Line(iFrst:iEnd)
            End If
*
*           Modify the lines with Fragments to Z-matrix format
*
            iFrst=iEq+1
            Call NxtWrd(Line,iFrst,iEnd)
            If (Line(iFrst:iFrst+7).eq.'FRAGMENT') Then
*
*              Here we do the expansion
*
               nFrag=nFrag+1
               FragZMat(nFrag,1)=0
               FragZMat(nFrag,2)=0
               iFrst=1
               Call NxtWrd(Line,iFrst,iEnd)
               FragLabels(nFrag)=Line(iFrst:iEnd)
               iFrst=Index(Line,'FRAGMENT')
               Call NxtWrd(Line,iFrst,iEnd)
*
*              Pick up the first atom label
*
               iFrst=iEnd+1
               Call NxtWrd(Line,iFrst,iEnd)
               If (Line(iFrst:iEnd).eq.'&') Then
                  iLine=iLine+1
                  Line=Lines(iLine)     ! Read the continuation line
                  iFrst=1
                  Call NxtWrd(Line,iFrst,iEnd)
               End If
               Label1=Line(iFrst:iEnd)
               nLabel1=iEnd-iFrst+1
*
*              Pick up the second atom label
*
               iFrst=iEnd+1
               Call NxtWrd(Line,iFrst,iEnd)
               If (Line(iFrst:iEnd).eq.'&') Then
                  iLine=iLine+1
                  Line=Lines(iLine)     ! Read the continuation line
                  iFrst=1
                  Call NxtWrd(Line,iFrst,iEnd)
               End If
               Label2=Line(iFrst:iEnd)
               nLabel2=iEnd-iFrst+1
*
*              Case: Bond
*
               nZMat=nZMat+1
               Write (Line2,'(A,I3.3,A,2(1X,A))')
     &         'ZMAT',nZMat+ZMatOffset,' = Bond',Label2(1:nLabel2),
     &                                           Label1(1:nLabel1)

               Write(Lu_TMP,'(A)') Line2
               iRow_c=iRow_c+1
               iZMat=iZMat+1
               FragZMat(nFrag,1)=iZMat
               FragZMat(nFrag,2)=iZMat
*
*              Pick up the third atom label
*
               iFrst=iEnd+1
               Call NxtWrd(Line,iFrst,iEnd)
               If (iEnd.eq.-1) Go To 300
               If (Line(iFrst:iEnd).eq.'&') Then
                  iLine=iLine+1
                  Line=Lines(iLine)     ! Read the continuation line
                  iFrst=1
                  Call NxtWrd(Line,iFrst,iEnd)
               End If
               Label3=Line(iFrst:iEnd)
               nLabel3=iEnd-iFrst+1
               jAtom=0
               Do iAtom = 1, nsAtom
                  If (Label3.eq.AtomLbl(iAtom)) jAtom=iAtom
               End Do
               If (nStab(jAtom).eq.1) Then
                  nDeg=3
               Else If (nStab(jAtom).eq.2) Then
                  nDeg=2
               Else If (nStab(jAtom).eq.1) Then
                  nDeg=1
               Else
                  nDeg=0
               End If
*
*              Case: Bond-Angle
               If (nDeg.le.1) Go To 200
               nZMat=nZMat+1
               Write (Line2,'(A,I3.3,A,2(1X,A))')
     &         'ZMAT',nZMat+ZMatOffset,' = Bond',Label3(1:nLabel3),
     &                                           Label2(1:nLabel2)
               Write(Lu_TMP,'(A)') Line2
               iRow_c=iRow_c+1
*
               If (nDeg.le.2) Go To 200
               nZMat=nZMat+1
               Write (Line2,'(A,I3.3,A,3(1X,A))')
     &         'ZMAT',nZMat+ZMatOffset,' = Angle',Label3(1:nLabel3),
     &                                            Label2(1:nLabel2),
     &                                            Label1(1:nLabel1)
               Write(Lu_TMP,'(A)') Line2
               iRow_c=iRow_c+1
               iZMat=iZMat+2
               FragZMat(nFrag,2)=iZMat
 200           Continue
*
 400           Continue
*
*              Pick up next atom (fourth)
*
               iFrst=iEnd+1
               Call NxtWrd(Line,iFrst,iEnd)
*
*              Branch out if no more labels. Look out for line
*              continuations
*
               If (iEnd.eq.-1) Go To 300
               If (Line(iFrst:iEnd).eq.'&') Then
                  iLine=iLine+1
                  Line=Lines(iLine)     ! Read the continuation line
                  iFrst=1
                  Call NxtWrd(Line,iFrst,iEnd)
               End If
*
               Label4=Line(iFrst:iEnd)
               nLabel4=iEnd-iFrst+1
               jAtom=0
               Do iAtom = 1, nsAtom
                  If (Label4.eq.AtomLbl(iAtom)) jAtom=iAtom
               End Do
               If (nStab(jAtom).eq.1) Then
                  nDeg=3
               Else If (nStab(jAtom).eq.2) Then
                  nDeg=2
               Else If (nStab(jAtom).eq.1) Then
                  nDeg=1
               Else
                  nDeg=0
               End If
*
*              Case: Bond-Angle-Dihedral
               If (nDeg.eq.0) Go To 500
               nZMat=nZMat+1
               Write (Line2,'(A,I3.3,A,2(1X,A))')
     &         'ZMAT',nZMat+ZMatOffset,' = Bond',Label4(1:nLabel4),
     &                                           Label3(1:nLabel3)
               Write(Lu_TMP,'(A)') Line2
               iRow_c=iRow_c+1
*
               If (nDeg.eq.1) Go To 500
               nZMat=nZMat+1
               Write (Line2,'(A,I3.3,A,3(1X,A))')
     &         'ZMAT',nZMat+ZMatOffset,' = Angle',Label4(1:nLabel4),
     &                                            Label3(1:nLabel3),
     &                                            Label2(1:nLabel2)
               Write(Lu_TMP,'(A)') Line2
               iRow_c=iRow_c+1
*
               If (nDeg.eq.2) Go To 500
               nZMat=nZMat+1
               Write (Line2,'(A,I3.3,A,4(1X,A))')
     &         'ZMAT',nZMat+ZMatOffset,' = Dihedral',Label4(1:nLabel4),
     &                                               Label3(1:nLabel3),
     &                                               Label2(1:nLabel2),
     &                                               Label1(1:nLabel1)
               Write(Lu_TMP,'(A)') Line2
               iRow_c=iRow_c+1
               iZMat=iZMat+3
               FragZMat(nFrag,2)=iZMat
 500           Continue
*
*              Push the labels
*
               Label1=Label2
               nLabel1=nLabel2
               Label2=Label3
               nLabel2=nLabel3
               Label3=Label4
               nLabel3=nLabel4
*
*              Get the next label
*
               Go To 400
*
 300           Continue
*
            Else
*
*              No processing just write to the file.
*
               Write(Lu_TMP,'(A)') Line
               iRow_c=iRow_c+1
            End If
*
            Else ! Values
*
*              Remove constraints to be ignored
*
               Ignore=.False.
               iFrst=1
               Call NxtWrd(Lines(1),iFrst,iEnd)
               iFrag=0
               Do i=1,nFrag
                  If (Line(iFrst:iEnd).eq.Trim(FragLabels(i))) Then
                     iFrag=i
                  End If
               End Do
*              Some constraints default to soft unless marked as hard
               iSoft=0
               If (Index(Lines(nLines),'HARD').le.iEq) Then
                  Do i=1,nSoft
                     If (Line(iFrst:iEnd).eq.Trim(SoftLabels(i))) Then
                        iSoft=i
                     End If
                  End Do
               End If
               iEq=Index(Lines(nLines),'=')
*              Ignore soft in num. grad., phantom in slapaf
               If ((SuperName.eq.'numerical_gradient') .and.
     &             ((Index(Lines(nLines),'SOFT').gt.iEq) .or.
     &              (iSoft.gt.0))) Ignore=.True..And.Remove
               If ((SuperName.eq.'slapaf') .and.
     &             (Index(Lines(nLines),'PHANTOM').gt.iEq))
     &               Ignore=.True..And.Remove
*              Ignore if this is a fragment constraint
               If (iFrag.gt.0) Then
*                 if it is already ignored, ignore all replacement constraints
                  If (Ignore) FragZMat(iFrag,2)=0
                  Ignore=.True.
               End If
               If (.Not.Ignore) Then
                  Do i=1,nLines
                     Write(Lu_TMP,'(A)') Lines(i)
                  End Do
                  iRow_c=iRow_c+nLines
               Else
                  nLambda=nLambda-1
               End If
*
            End If ! Values
*
            Go To 100 ! Get the next line
*
         End If
*
*     Put in the additional constraints here.
*
      Do iFrag=1,nFrag
*        Note that FragZMat(iFrag,2)=0 if the fragment is ignored
         Do iZMat=FragZMat(iFrag,1),FragZMat(iFrag,2)
            Write (Lu_TMP,'(A,I3.3,A)') 'ZMAT',iZMat+ZMatOffset,' = FIX'
            iRow_c=iRow_c+1
            nLambda=nLambda+1
         End Do
      End Do
*
*                                                                      *
************************************************************************
*                                                                      *
*     Write out the last line again
*
      Line='END OF CONSTRAINTS'
      Write(Lu_TMP,'(A)') Line
*                                                                      *
************************************************************************
*                                                                      *
*     Copy the stuff from LU_TMP back to LU_UDC
*
      REWIND(Lu_UDC)
      REWIND(Lu_TMP)
 600  Continue
      Line=Get_Ln(Lu_TMP)
*     Write (*,*) Line
*
      Write (Lu_UDC,'(A)') Trim(Line)
      If (Line(1:4).ne.'END ') Go To 600
*     Write (*,*) 'iRow_C,nLambda=',iRow_C,nLambda
*
*                                                                      *
************************************************************************
*                                                                      *
*     Close the files.
*
      Close(Lu_TMP)
      Close(Lu_UDC)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
