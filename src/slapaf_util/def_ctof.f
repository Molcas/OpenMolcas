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
*               2008, Giovanni Ghigo                                   *
************************************************************************
      SubRoutine Def_CtoF(lNew,nAtom,Coor)
************************************************************************
*                                                                      *
*     Author: Giovanni Ghigo, Dep. of General and Organic Chemistry    *
*             University of Torino, ITALY                              *
*             July 2008                                                *
*     Adapted from  DefInt by                                          *
*             Roland Lindh, Dep. of Theoretical Chemistry,             *
*             University of Lund, SWEDEN                               *
************************************************************************
      use Slapaf_Info, only: dMass
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "Molcas.fh"
      Real*8    Coor(3,nAtom)
      Character Labels*8, Type*6, Temp*120,
     &          Line*120, Format*8, filnam*16
      Logical lWrite, lNew
      Integer, Allocatable:: Ind(:,:)
      Real*8, Allocatable:: xyz(:,:), Temp2(:,:), Mass(:,:)
*
      lWrite = .True.
      Lu=6
      nTemp=Len(Temp)
      Write (Format,'(A,I3.3,A)') '(F',nTemp,'.0)'
*
      Lu_UDIC=91
      filnam='UDIC'
      call molcas_open(Lu_UDIC,filnam)
c      Open(Lu_UDIC,File=filnam,Form='Formatted',Status='OLD')
      Rewind(Lu_UDIC)
      Write (6,*)
      Write (6,*)
     &'****************************************************************'
      If (lNew) then
         Write (6,*)
     &'* New value of the internal coordinate to follow               *'
      else
         Write (6,*)
     &'* Original value of the internal coordinate to follow          *'
      EndIf
      Write (6,*)
     &'****************************************************************'
*
*     Step 1. BSet up the b vectors from which we will define the
*     internal coordinates.
*
c      iBVct = 0
      Read (Lu_UDIC,'(A)') Line
      Temp=Line
      Call UpCase(Temp)
*
*     Move the label of the internal coordinate
*
      neq = Index(Line,'=')
      If (neq.Eq.0) Then
         Call WarningMessage(2,'Error in Def_CTOF')
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
            Call WarningMessage(2,'Error in Def_CTOF')
            Write (6,'(A)') '***********************************'
            Write (6,'(A)') ' Syntax error in line :            '
            Write (6,'(A)') Line(1:33),'...'
            Write (6,'(A,A)') Line(iFrst:jEnd),
     &            ' has more than 8 character'
            Write (6,'(A)') '***********************************'
            Call Quit_OnUserError()
         End If
         Labels = Line(iFrst:jEnd)
      End If
*
*-----Construct the corresponding transformation vector
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
            Call WarningMessage(2,'Error in Def_CTOF')
            Write (6,*) 'DefInt: wrong cartesian type'
            Write (6,'(A,A)') 'Temp=',Temp
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
            Call WarningMessage(2,'Error in Def_CTOF')
            Write (6,*)' Line contains syntax error!'
            Write (6,'(A)') Line
            Write (6,*) i1,i2,i3
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
         Call WarningMessage(2,'Error in Def_CTOF')
         Write (6,*)' Line contains syntax error!'
         Write (6,'(A)') Line
         Call Quit_OnUserError()
      End If
*
      msAtom = nCntr + mCntr
      Call mma_allocate(xyz ,3,msAtom,Label='xyz')
      Call mma_allocate(Temp2,3,msAtom,Label='Temp2')
      Call mma_allocate(Ind ,2,msAtom,Label='Ind')
      Call mma_allocate(Mass,2,msAtom,Label='Mass')
*
      Call CllCtoF(Line(nGo:nTemp),nAtom,Coor,nCntr,mCntr,xyz,
     &             Temp2,Ind,Type,dMass,Mass,Labels)
*
      Call mma_deallocate(Mass)
      Call mma_deallocate(Ind)
      Call mma_deallocate(Temp2)
      Call mma_deallocate(xyz)
*
      Close(Lu_UDIC)
      Return
      End
