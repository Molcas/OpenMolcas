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
      SubRoutine ChTab(iOper,nIrrep,lIrrep,lBsFnc,iSigma)
************************************************************************
*                                                                      *
* Object: to generate the character table of a point group within      *
*         D2h.                                                         *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             September '91                                            *
************************************************************************
      use Symmetry_Info, only: Symmetry_Info_Set
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Integer iOper(nIrrep), iChTbl(1:8,1:8)
      Integer iTest(8)
      Character*80 lIrrep(8)*3, lBsFnc(8), Tmp
      Character*6 xyz(0:7), SymLab*3
      Common /SymLab/SymLab
      Logical Inv, Rot
      Data xyz/'      ','x     ','y     ','xy, Rz',
     &         'z     ','xz, Ry','yz, Rx','I     '/
*                                                                      *
************************************************************************
*                                                                      *
*     Call qEnter('ChTab')
*
      Do 1 i = 1, 8
         lIrrep(i) = ' '
         lBsFnc(i) = ' '
 1    Continue
      If (nIrrep.eq.1) Then
         SymLab='C1 '
         iSigma=1
      Else If (nIrrep.eq.2) Then
         If (iOper(2).eq.7) Then
            SymLab='Ci '
            iSigma=1
         Else If (iOper(2).eq.1.or.iOper(2).eq.2.or.iOper(2).eq.4) Then
            SymLab='Cs'
            iSigma=1
         Else
            SymLab='C2'
            iSigma=2
         End If
      Else If (nIrrep.eq.4) Then
         If (iOper(2).eq.7.or.iOper(3).eq.7.or.iOper(4).eq.7) Then
            SymLab='C2h'
            iSigma=2
         Else
            Rot = .True.
            Do 600 i = 1, nIrrep
               If (iOper(i).eq.1.or.iOper(i).eq.2.or.iOper(i).eq.4)
     &            Rot = .False.
 600        Continue
            If (Rot) Then
               SymLab='D2 '
               iSigma=2
            Else
               SymLab='C2v'
               iSigma=2
            End If
         End If
      Else If (nIrrep.eq.8) Then
         SymLab='D2h'
         iSigma=2
      Else
         Call WarningMessage(2,'ChTab: Illegal value of nIrrep')
         Write (6,*) 'nIrrep=',nIrrep
         Call Abend()
      End If
      Call ICopy(8**2,[0],0,iChTbl,1)
*
*     Go through the functions x, y, and z, and the dyadic functions.
*
      iSymX = 0
      iSymY = 0
      iSymZ = 0
      Do 104 i = 1, nIrrep
         If (iAnd(iOper(i),1).ne.0) iSymX = 1
         If (iAnd(iOper(i),2).ne.0) iSymY = 2
         If (iAnd(iOper(i),4).ne.0) iSymZ = 4
 104  Continue
*
*-----Loop over basis functions (a' la Malmqvist)
*
      Do 10 iFnc = 0, 7
         Tmp=xyz(iFnc)
*
*        Generate a row in the character table of this function
*
         ix = iAnd(iFnc,iSymX)
         iy = iAnd(iFnc,iSymY)/2
         iz = iAnd(iFnc,iSymZ)/4
*--------Loop over all operators
         Do 20 i = 1, nIrrep
            jx = iAnd(iOper(i),iSymX)
            jy = iAnd(iOper(i),iSymY)/2
            jz = iAnd(iOper(i),iSymZ)/4
            iCh = 1
            If (ix.ne.0 .and. jx.ne.0) iCh = -iCh
            If (iy.ne.0 .and. jy.ne.0) iCh = -iCh
            If (iz.ne.0 .and. jz.ne.0) iCh = -iCh
            iTest(i) = iCh
 20      Continue
*
*--------Compute place of Irrep
*
         If (nIrrep.eq.1) Then
            jIrrep=1
         Else If (nIrrep.eq.2) Then
            jIrrep=1+(1-iTest(2))/2
         Else If (nIrrep.eq.4) Then
            jIrrep=1+( (1-iTest(2))+2*(1-iTest(3)) )/2
         Else If (nIrrep.eq.8) Then
            jIrrep=1+( (1-iTest(2))+2*(1-iTest(3))+4*(1-iTest(5)) )/2
         Else
            jIrrep=-1
            Call WarningMessage(2,'ChTab: Illegal nIrrep value!')
            Write (6,*) 'nIrrep=',nIrrep
            Call Abend()
         End If
         If (lBsFnc(jIrrep)(1:1).eq.' ') Then
            lBsFnc(jIrrep) = Tmp
            Call ICopy(nIrrep,iTest,1,iChTbl(jIrrep,1),8)
         Else
            LenlBs=Len(lBsFnc(jIrrep))
            LenTmp=Len(Tmp)
            i1 = iCLast(lBsFnc(jIrrep),LenlBs)
            i2 = iCLast(Tmp,LenTmp)
            lBsFnc(jIrrep) = lBsFnc(jIrrep)(1:i1)//', '//Tmp(1:i2)
         End If
 10   Continue
*
*     Set up some Mulliken symbols for the irreps
*
      Do 100 iIrrep = 1, nIrrep
         lIrrep(iIrrep)='a'
         Do 110 i = 1, nIrrep
*           Write (*,*) ' iIrrep,i=',iIrrep,i
*
*           If the character of an rotation in an irreps is -1 then
*           the irreps is assigned the character B, otherwise A.
*
*           Write (*,*) iOper(i),iChTbl(iIrrep,i)
            If ((iOper(i).eq.3 .or. iOper(i).eq.5 .or.
     &           iOper(i).eq.6) .and. iChTbl(iIrrep,i).eq.-1)
     &          lIrrep(iIrrep)='b'
*
 110     Continue
 100  Continue
      iSub = 0
*
*     Subscript according to C2 operations
*
      Rot = .False.
      Do 300 i = 1, nIrrep
         If (iOper(i).eq.3 .or.
     &       iOper(i).eq.5 .or.
     &       iOper(i).eq.6) Rot = .True.
 300  Continue
      If (Rot.and.SymLab.ne.'C2v') Then
         iSub = iSub + 1
*
*        Find the number of A's and B's
*
         ia = 0
         ib = 0
         Do 310 i = 1, nIrrep
            If (lIrrep(i)(1:1).eq.'a') ia = ia + 1
            If (lIrrep(i)(1:1).eq.'b') ib = ib + 1
 310     Continue
         If (nIrrep.eq.8) Then
            ia = ia/2
            ib = ib/2
         End If
         If (SymLab.eq.'C2h') Then
            ia = ia/2
            ib = ib/2
         End If
*
*        Find the rotations
*
         iRot = 0
         Do 320 i = 1, nIrrep
            If ( iOper(i).eq.3.or.iOper(i).eq.5.or.iOper(i).eq.6) Then
               iRot = iRot + 1
               Write (Tmp,'(I1)') iRot
               If (ia.gt.1) Then
                  Do 321 j = 1, nIrrep
                     If (lIrrep(j)(1:1).eq.'a'.and.iChTbl(j,i).eq.1)
     &                  lIrrep(j)=lIrrep(j)(1:1)//Tmp(1:1)
 321              Continue
               End If
               If (ib.gt.1) Then
                  Do 322 j = 1, nIrrep
                     If (lIrrep(j)(1:1).eq.'b'.and.iChTbl(j,i).eq.1)
     &                  lIrrep(j)=lIrrep(j)(1:1)//Tmp(1:1)
 322              Continue
               End If
            End If
 320     Continue
      Else If (Rot.and.SymLab.eq.'C2v') Then
*
*        Find the Rotation
*
         iRot = -1
         Do 350 i = 1, nIrrep
            If (iOper(i).eq.3.or.iOper(i).eq.5.or.iOper(i).eq.6)
     &         iRot = iOper(i)
 350     Continue
*
*        Find the first vertical mirror plane to this axis
*
         Do 351 i = 1, nIrrep
            If (iOper(i).ne.3.and.iOper(i).ne.5.and.iOper(i).ne.6
     &         .and.iOper(i).ne.7.and.iAnd(iOper(i),iRot).ne.1)
     &         iRot = i
 351     Continue
         Do 352 i = 1, nIrrep
            If (iChTbl(i,iRot).eq.1) Then
               j = 1
            Else
               j = 2
            End If
            Write (Tmp,'(I1)') j
            lIrrep(i)=lIrrep(i)(1:1)//Tmp(1:1)
 352     Continue
      End If
*
*     Subscript according to inversion if present
*
      Inv=.False.
      Do 200 i = 1, nIrrep
         Inv = iOper(i).eq.7 .or. Inv
 200  Continue
      If (Inv) Then
         iSub = iSub + 1
*
*------- Loop over each Irrep
*
         Do 210 iIrrep = 1, nIrrep
            LenlIrr=Len(lIrrep(iIrrep))
            i1 = 1 + iCLast(lIrrep(iIrrep),LenlIrr)
*
*---------- Loop over operators
*
            Do 211 i = 1, nIrrep
               If (iOper(i).eq.7) Then
                  If (iChTbl(iIrrep,i).eq.1) Then
                     lIrrep(iIrrep)(i1:i1)='g'
                  Else If (iChTbl(iIrrep,i).eq.-1) Then
                      lIrrep(iIrrep)(i1:i1)='u'
                  End If
               End If
 211        Continue
 210     Continue
      End If
*
*     Fix labels for Cs
*
      If (SymLab(1:2).eq.'Cs') Then
         lIrrep(1) = 'a'''
         lIrrep(2) = 'a"'
      End If
*
      Call Symmetry_Info_Set(nIrrep,iOper,iChTbl)
*                                                                      *
************************************************************************
*                                                                      *
*     Call qExit('ChTab')
      Return
      End
