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
      Subroutine Get_Mpprop_input(nAtoms,iPol,LNearestAtom,LAllCenters
     &,AveOrb,LLumOrb,Diffuse,dLimmo,Thrs1,Thrs2,nThrs,ThrsMul,iPrint)
*
      Implicit Real*8 (a-h,o-z)
*
#include "MolProp.fh"
#include "warnings.fh"

      Integer      nBondCount(mxAtomMP)
      Character*3  EndKey
*Jose Character*4  TestLabe(0:nAtoms), KWord
      Character*4   KWord
      Character*6  TestLabe(0:nAtoms)
      Character*180 Key, BLine
      Character*180 Get_Ln
      Logical Debug, lTtl, LNearestAtom
      Logical LAllCenters, AveOrb, Diffuse(3)
      Logical LLumorb
      Dimension dLimmo(2)

      External Get_Ln

      Data Debug/.False./

      iStdOut = 6

      Do i = 1, 180
         BLine(i:i) = ' '
      End Do
      Title=BLine

*
*     Call qEnter('RdCtl')
*
*

      LuRd=21
      Call SpoolInp(LuRd)

      Rewind(LuRd)
      Call RdNLst(LuRd,'MPPROP')
*
*     KeyWord directed input
*
 998  lTtl = .False.
!EB 9988 Continue
      Key = Get_Ln(LuRd)
      If (Debug) Write (iStdOut,*) ' Processing:',Key
      KWord = Trim(Key)
      Call UpCase(KWord)
      If (KWord(1:1).eq.'*')    Go To 998
      If (KWord.eq.BLine)       Go To 998
      If (KWord(1:4).eq.'BOND') Go To 981
!      If (KWord(1:4).eq.'METH') Go To 982
      If (KWord(1:4).eq.'TITL') Go To 983
      If (KWord(1:4).eq.'TYPE') Go To 984
      If (KWord(1:4).eq.'POLA') Go To 985
      If (KWord(1:4).eq.'NONE') Go To 986
      If (KWord(1:4).eq.'ALLC') Go To 987
      If (KWord(1:4).eq.'LUMO') Go To 988
      If (KWord(1:4).eq.'DIFF') Go To 989
      If (KWord(1:4).eq.'PRIN') Go To 991
! Keyword added to handle average orbitals
      If (KWord(1:4).eq.'AVER') Go To 7912
*
      If (KWord(1:4).eq.'END ') Go To 997
      iChrct=Len(KWord)
      Last=iCLast(KWord,iChrct)
      Write (6,*)
      Write (6,'(1X,A,A)') KWord(1:Last),' is not a keyword!'
      Write (6,*) ' Error in keyword.'
      Call ErrTra
      Call Quit(_RC_INPUT_ERROR_)
      Write (6,*) ' Premature end of input file.'
      Call Quit(_RC_INPUT_ERROR_)
      Write (6,*) ' Error while reading input file.'
 990  Call ErrTra
      Call Quit(_RC_INPUT_ERROR_)
*                                                                      *
************************************************************************
*                                                                      *
* Get the input for bonds
981   LAllCenters=.True.
      Do i=1,MxAtomMP
         nBondCount(i)=0
         NUB(i) = 0
         Do j=1,MxAtomMP
            NBI(i,j) = 0
            BondMat(i,j)=.False.
         EndDo
      EndDo
      Do i=1,nAtoms
         nBonds = 0
         Key = Get_Ln(LuRd)
         m = 1
         Do j=1,nAtoms
            TestLabe(j) = ' '
         EndDo
         Do j=1,180
            EndKey=Key(j:j+2)
            Call UpCase(EndKey)
            If((Key(j:j).eq.' ').or.(Key(j:j).eq.',').or.
     &      (Key(j:j).eq.';')) Then
               If((j.eq.1).and.
     &         ((Key(j:j).eq.';').or.(Key(j:j).eq.','))) Then
                  Write(iStdOut,*)
     &            'Error in input, breaker in first position'
                  Goto 990
               ElseIf(m.lt.j) Then
                  TestLabe(nBonds)=Key(m:j-1)
                  nBonds = nBonds+1
                  m = j+1
               Else
                  m = j+1
               EndIf
            ElseIf (EndKey.eq.'END') Then
               Goto 9811
            EndIf
         EndDo
         Do j=1,nAtoms
            If(TestLabe(0).eq.Labe(j)) Then
               Do k=1,nAtoms
                  Do l=1,nAtoms
                     If(TestLabe(k).eq.Labe(l)) Then
                        BondMat(j,l) = .True.
                        BondMat(l,j) = .True.
                     Else

                     EndIf
                  EndDo
               EndDo
            EndIf
         EndDo
      EndDo
9811  Continue

      Do i=1,nAtoms
         Do j=1,nAtoms
            If(BondMat(i,j)) Then
               NuB(i)=NuB(i)+1
               NBI(i,NuB(i))=j
            EndIf
         EndDo
      EndDo

      Write(iStdOut,*)
      Write(iStdOut,'(10X,A)') '************************'
      Write(iStdOut,'(10X,A)') '**** Bonding matrix ****'
      Write(iStdOut,'(10X,A)') '************************'
      Write(iStdOut,*)
      Write(iStdOut,'(A8,A,A)') 'Atom', '  No bonds', '   Bonding with'
      Do i=1,nAtoms
         Write(iStdOut,'(A8,I6,A11,1000A8)') LABE(I),NUB(I),
     &   (LABE(NBI(I,J)),J=1,NUB(I))
      EndDo
      Write(iStdOut,*)
      Write(iStdOut,*)
      Goto 998

* Set Method level
!982   Key = Get_Ln(LuRd)
!      If (Debug) Write (*,*) ' Processing:',Key
!      Method = Key
!      Call UpCase(Method)
!      If( (Method=='SCF') .or. (Method=='RASSCF') ) Then
!
!      Else
!         Write(*,*) 'The Method Label is not correct'
!         Write(*,*) Method
!         Call Quit(20)
!      EndIf
!      Goto 998
* Set the Title
983   Key = Get_Ln(LuRd)
      If (Debug) Write (iStdOut,*) ' Processing:',Key
      Call UpCase(Key)
      Title=Key
      Goto 998
* Set the specifik atom types
984   Key = Get_Ln(LuRd)
      Call UpCase(Key)
      If (Key(1:3).eq.'END') Then
         Go To 998
      Else
         Read(Key,*) j,iAtomPar(j)
      EndIf
      Goto 984
* Set that local polarizabilities should be calculated
985   Key = Get_Ln(LuRd)
      Read(Key,*) iPol
      Goto 998
* Set the nearest atom calculation to .False.
986   LNearestAtom=.False.
      Goto 998
* Set the average orbital option to .True.
7912  AveOrb=.True.
      iPol=0
      Goto 998
* Set the logical variable that defines if all ceters or not
987   LAllCenters=.True.
      LNearestAtom=.False.
*                                                                      *
************************************************************************
*                                                                      *
*---- Set the bonds to all sites
*---- Get information from input
!      Do i=1,mxAtomMP
!         NuB(i)=0
!         Do j=1,mxAtomMP
!            NBI(i,j)=0
!         End Do
!      End Do
      Do i=1,nAtoms
         NuB(i)=nAtoms-1
         Do j=1,nAtoms-1
            If(j.ge.i) Then
               NBI(i,j)=j+1
               BondMat(i,j+1) = .True.
            Else
               NBI(i,j)=j
               BondMat(i,j) = .True.
            EndIf
         End Do
      End Do
      Goto 998
* Get the vectors from the INPORB file
988   LLumOrb=.True.
      Goto 998

* Obtain diffuse stuff to the MpProp decomposition.
989   Continue
      Key = Get_Ln(LuRd)
      Call UpCase(Key)
      If(Key(1:4).eq.'NUME') then
        Diffuse(1)=.true.
        Diffuse(2)=.true.
9891    Continue
        Key = Get_Ln(LuRd)
        Call UpCase(Key)
        If(Key(1:4).eq.'LIMI') then
          Key = Get_Ln(LuRd)
          Call Get_F(1,dLimmo,2)
        Elseif(Key(1:4).eq.'THRE') then
          Key=Get_Ln(LuRd)
          Call Get_F(1,Thrs1,1)
          Call Get_F(2,Thrs2,1)
          Call Get_I(3,nThrs,1)
          Call Get_F(4,ThrsMul,1)
        Elseif(Key(1:4).eq.'END ') then
          GoTo 987 !Not an error, Diffuse implies AllCenters
        Else
          Write(6,*) 'Undefined option for ''DIFFuse'':',Key
          Call FindErrorLine
          Call Quit_OnUserError()
        Endif
        GoTo 9891
      Elseif(Key(1:4).eq.'REXT') then
        Diffuse(1)=.true.
        Diffuse(3)=.true.
      Else
        Write(6,*) 'Undefined option for ''DIFFuse'':',Key
        Call FindErrorLine
        Call Quit_OnUserError()
      Endif
      GoTo 987  !Not an error, Diffuse implies AllCenters

*PrintLevel.
991   Continue
      Key = Get_Ln(LuRd)
      Call UpCase(Key)
      Call Get_I(1,iPrint,1)
      GoTo 998

997   Continue
      If(Title.eq.Bline) Then
         Write(iStdOut,*)
         Write(iStdOut,*) ' !!WARNING!! The molecule do not have a name'
         Write(iStdOut,*)
      EndIf

      Return
      End
