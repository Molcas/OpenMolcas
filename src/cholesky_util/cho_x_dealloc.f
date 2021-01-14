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
* Copyright (C) 2004, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine Cho_X_Dealloc(irc)
      use ChoArr, only: iSOShl, iBasSh, nBasSh, nBstSh, iSP2F, iAtomShl,
     &                  iShlSO
C
C     T.B. Pedersen, July 2004.
C
C     Purpose: deallocate ALL index arrays for the Cholesky utility.
C              On exit, irc=0 signals sucessful completion.
C
      Implicit None
      Integer irc
#include "choptr.fh"
#include "chosew.fh"
#include "cholq.fh"
#include "chopar.fh"
#include "stdalloc.fh"

      Character*13 SecNam
      Parameter (SecNam = 'Cho_X_Dealloc')

      Integer nAlloc

C     Initialize allocation counter.
C     ------------------------------

      nAlloc = 0

C     Deallocate.
C     -----------

      If (l_InfRed .ne. 0) Then
         Call GetMem('InfRed','Free','Inte',ip_InfRed,l_InfRed)
      End If
      nAlloc = nAlloc + 1

      If (l_InfVec .ne. 0) Then
         Call GetMem('InfVec','Free','Inte',ip_InfVec,l_InfVec)
      End If
      nAlloc = nAlloc + 1

      If (l_IndRed .ne. 0) Then
         Call GetMem('IndRed','Free','Inte',ip_IndRed,l_IndRed)
      End If
      nAlloc = nAlloc + 1

      If (l_IndRSh .ne. 0) Then
         Call GetMem('IndRSh','Free','Inte',ip_IndRSh,l_IndRSh)
      End If
      nAlloc = nAlloc + 1

      If (l_iScr .ne. 0) Then
         Call GetMem('iScr','Free','Inte',ip_iScr,l_iScr)
      End If
      nAlloc = nAlloc + 1

      If (l_iiBstRSh .ne. 0) Then
         Call GetMem('iiBstRSh','Free','Inte',ip_iiBstRSh,l_iiBstRSh)
      End If
      nAlloc = nAlloc + 1

      If (l_nnBstRSh .ne. 0) Then
         Call GetMem('nnBstRSh','Free','Inte',ip_nnBstRSh,l_nnBstRSh)
      End If
      nAlloc = nAlloc + 1

      If (l_IntMap .ne. 0) Then
         Call GetMem('IntMap','Free','Inte',ip_IntMap,l_IntMap)
      End If
      nAlloc = nAlloc + 1

      If (l_nDimRS .ne. 0) Then
         Call GetMem('nDimRS','Free','Inte',ip_nDimRS,l_nDimRS)
      End If
      nAlloc = nAlloc + 1

      If (l_iRS2F .ne. 0) Then
         Call GetMem('iRS2F','Free','Inte',ip_iRS2F,l_iRS2F)
      End If
      nAlloc = nAlloc + 1

      If (Allocated(iSOShl)) Call mma_deallocate(iSOShl)

      If (Allocated(iShlSO)) Call mma_deallocate(iShlSO)

      If (l_iQuab .ne. 0) Then
         Call GetMem('iQuab','Free','Inte',ip_iQuab,l_iQuab)
      End If
      nAlloc = nAlloc + 1

      If (Allocated(iBasSh)) Call mma_deallocate(iBasSh)

      If (Allocated(nBasSh)) Call mma_deallocate(nBasSh)

      If (Allocated(nBstSh)) Call mma_deallocate(nBstSh)

      If (Allocated(iAtomShl)) Call mma_deallocate(iAtomShl)

      If (Allocated(iSP2F)) Call mma_deallocate(iSP2F)

C     Check that #allocations agrees with choptr.fh.
C     -----------------------------------------------

      irc = CHO_NALLOC - nAlloc
      If (irc .ne. 0) Then
         Write(6,*) SecNam,' is out of sync with choptr.fh !!!'
         Write(6,*) '(Note that this is due to a programming error...)'
         Return
      End If

C     Zero entire common block.
C     -------------------------

      Call Cho_PtrIni(irc)
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': Cho_PtrIni is out of sync ',
     &              'with choptr.fh !!!'
         Write(6,*) '(Note that this is due to a programming error...)'
         Return
      End If

C     Deallocate any used pointer in chosew.fh
C     -----------------------------------------

      If (l_iShP2RS .ne. 0) Then
         Call GetMem('SHP2RS','Free','Inte',ip_iShP2RS,l_iShP2RS)
         ip_iShP2RS=0
         l_iShP2RS=0
      End If

      If (l_iShP2Q .ne. 0) Then
         Call GetMem('SHP2Q','Free','Inte',ip_iShP2Q,l_iShP2Q)
         ip_iSHP2Q=0
         l_iSHP2Q=0
      End If

C     Deallocate any used pointer in cholq.fh
C     ----------------------------------------

      If (l_iQuAB_L .ne. 0) Then
         Call GetMem('IQUAB_L','Free','Inte',ip_iQuAB_L,l_iQuAB_L)
         ip_iQuAB_L=0
         l_iQuAB_L=0
      End If

      If (l_iQL2G .ne. 0) Then
         Call GetMem('IQL2G','Free','Inte',ip_iQL2G,l_iQL2G)
         ip_iQL2G=0
         l_iQL2G=0
      End If

      If (l_LQ .ne. 0) Then
         Call GetMem('LQ','Free','Real',ip_LQ,l_LQ)
         ip_LQ=0
         l_LQ=0
      End If

C     Deallocate any used pointer in chopar.fh
C     -----------------------------------------

      If (l_InfVec_Bak .gt. 0) Then
         Call GetMem('InfVec_Bak','Free','Inte',ip_InfVec_Bak,
     &                                           l_InfVec_Bak)
         l_InfVec_Bak=0
      End If

      Return
      End
