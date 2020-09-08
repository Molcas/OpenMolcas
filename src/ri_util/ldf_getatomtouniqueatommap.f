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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_GetAtomToUniqueAtomMap(A2UA,nA)
      Use Basis_Info
C
C     Thomas Bondo Pedersen, June 2010.
C     - based on Print_Geometry by Roland Lindh.
C
C     Compute map from atom to unique atom.
C
C     A2UA(i): returns unique atom corresponding to atom i
C              (i is an index in the LDF Atom Info list)
C
      Implicit Real*8 (A-H,O-Z)
      Integer A2UA(nA)
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"

#if defined (_DEBUG_)
      Character*62 Msg, Msg2, Msg3, Msg4
      Parameter (Msg=
     & 'LDF_GetAtomToUniqueAtomMap: LDF_AtomWithCoordinates returned 0')
      Parameter (Msg2=
     & 'LDF_GetAtomToUniqueAtomMap: LDF_AtomWithCoordinates too large ')
      Parameter (Msg3=
     & 'LDF_GetAtomToUniqueAtomMap: Internal error detected...        ')
      Parameter (Msg4=
     & 'LDF_GetAtomToUniqueAtomMap: Some atoms were not found...      ')
#endif
      Integer  LDF_AtomWithCoordinates
      External LDF_AtomWithCoordinates

#if defined (_DEBUG_)
      Call iZero(A2UA,nA)
      nCount=0
#endif

      ndc=0
      l_UAR=3
      Call GetMem('LDFUAR','Allo','Real',ip_UAR,l_UAR)
      Do jCnttp=1,nCnttp
         mCnt=dbsc(jCnttp)%nCntr
         If (dbsc(jCnttp)%pChrg .or.
     &       dbsc(jCnttp)%Aux .or.
     &       dbsc(jCnttp)%Frag) Then
            ndc=ndc+mCnt
         Else
            Do i=0,2
               Work(ip_UAR+i)=dbsc(jCnttp)%Coor(i+1,1)*
     &                               dble(iPhase(i+1,iCoset(0,0,ndc+1)))
            End Do
            ndc=ndc+1
            jxyz=jxyz+3
            iAtom=LDF_AtomWithCoordinates(Work(ip_UAR))
#if defined (_DEBUG_)
            nCount=nCount+1
            If (iAtom.lt.1) Then
               Call WarningMessage(2,Msg)
               Call LDF_Quit(1)
            End If
            If (iAtom.gt.nA) Then
               Call WarningMessage(2,Msg2)
               Call LDF_Quit(1)
            End If
            If (A2UA(iAtom).ne.0) Then
               Call WarningMessage(2,Msg3)
               Call LDF_Quit(1)
            End If
#endif
            iAtom_Unique=iAtom
            A2UA(iAtom)=iAtom_Unique
            Do jCnt=2,mCnt
               Do i=0,2
                  Work(ip_UAR+i)=dbsc(jCnttp)%Coor(i+1,jCnt)*
     &                               dble(iPhase(i+1,iCoset(0,0,ndc+1)))
               End Do
               ndc=ndc+1
               jxyz=jxyz+3
               iAtom=LDF_AtomWithCoordinates(Work(ip_UAR))
#if defined (_DEBUG_)
               nCount=nCount+1
               If (iAtom.lt.1) Then
                  Call WarningMessage(2,Msg)
                  Call LDF_Quit(1)
               End If
               If (iAtom.gt.nA) Then
                  Call WarningMessage(2,Msg2)
                  Call LDF_Quit(1)
               End If
               If (A2UA(iAtom).ne.0) Then
                  Call WarningMessage(2,Msg3)
                  Call LDF_Quit(1)
               End If
#endif
               A2UA(iAtom)=iAtom_Unique
            End Do
         End If
      End Do
      Call GetMem('LDFUAR','Free','Real',ip_UAR,l_UAR)

#if defined (_DEBUG_)
      If (nCount.ne.nA) Then
         Call WarningMessage(2,Msg4)
         Call LDF_Quit(1)
      End If
      nErr=0
      Do iA=1,nA
         If (A2UA(iA).lt.1 .or. A2UA(iA).gt.nA) Then
            nErr=nErr+1
         End If
      End Do
      If (nErr.ne.0) Then
         Call WarningMessage(2,Msg3)
         Call LDF_Quit(1)
      End If
#endif

      End
