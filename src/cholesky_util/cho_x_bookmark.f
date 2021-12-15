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
* Copyright (C) 2012, Thomas Bondo Pedersen                            *
************************************************************************
*  Cho_X_Bookmark
*
*> @brief
*>   Get number of Cholesky vectors needed to achieve a given integral target accuracy \p Thr
*> @author Thomas Bondo Pedersen, August 2012
*>
*> @details
*> Return number of Cholesky vectors in each symmetry block needed to
*> achieve an integral accuracy of \p Thr (&le; decomposition threshold).
*> The actual integral accuracy in each symmetry block is returned in
*> array \p delta (for 1C-CD, the accuracy is judged by 1-center
*> diagonal elements only!). Note that \p mSym may be smaller than \p nSym
*> (so that one may find a bookmark in irrep 1 only, for example).
*> Available for full CD and 1C-CD, but *not RI*.
*>
*> On exit, return codes signify:
*>
*> - \p irc = ``-1``: bookmarks not available
*> - \p irc =  ``0``: all OK
*> - \p irc =  ``1``: illegal input
*>
*> @param[in]  Thr   Target integral accuracy
*> @param[in]  mSym  Number of irreps
*> @param[out] nVec  Number of Cholesky vectors
*> @param[out] delta Actual integral accuracy
*> @param[out] irc   Return code
************************************************************************
      Subroutine Cho_X_Bookmark(Thr,mSym,nVec,delta,irc)
      Implicit None
      Real*8  Thr
      Integer mSym
      Integer nVec(mSym)
      Real*8  delta(mSym)
      Integer irc
#include "cho_para_info.fh"
#include "chobkm.fh"
#include "cholesky.fh"
#include "WrkSpc.fh"

      Character*14 SecNam
      Parameter (SecNam='Cho_X_Bookmark')

      Logical Found
      Logical DebugPrint
      Parameter (DebugPrint=.False.)

      Integer iSym
      Integer iRS
      Integer ip
      Integer l
      Integer n

      Integer i, j
      Integer nV
      Real*8  del
      nV(i,j)=iWork(ip_BkmVec-1+nRow_BkmVec*(j-1)+i)
      del(i,j)=Work(ip_BkmThr-1+nRow_BkmThr*(j-1)+i)

C     Set return code.
C     ----------------

      irc=0

C     Check that bookmarks are available.
C     -----------------------------------

      If (l_BkmVec.lt.1 .or. l_BkmThr.lt.1) Then
         irc=-1
         Return
      End If

C     Check input.
C     ------------

      If (sign(1.0d0,Thr).lt.0.0d0 .or. Thr.lt.ThrCom .or.
     &    mSym.lt.1 .or. mSym.gt.nSym) Then
         irc=1
         Return
      End If

C     Debug print.
C     ------------

      If (DebugPrint) Then
         Call Cho_Head(SecNam//': Bookmarks (nVec,delta)','-',80,6)
         Do iSym=1,nSym
            Write(6,'(A,I2,A)') 'Symmetry block',iSym,
     &                          ' Bookmarks (nVec,delta)'
            Write(6,'(5(1X,A,I6,A,D15.8,A))')
     &      ('(',nV(iRS,iSym),',',del(iRS,iSym),')',iRS=1,nRow_BkmThr)
         End Do
      End If

C     Find largest accuracy below Thr.
C     --------------------------------

      Do iSym=1,mSym
         iRS=0
         Found=.False.
         Do While (iRS.lt.nRow_BkmThr .and. .not.Found)
            iRS=iRS+1
            Found=del(iRS,iSym).le.Thr
         End Do
         If (.not.Found) Then
            Call Cho_Quit('Bug detected in '//SecNam,104)
         Else
            nVec(iSym)=nV(iRS,iSym)
            delta(iSym)=del(iRS,iSym)
         End If
      End Do

      ! Parallel run: take into account the distribution of vectors
      ! across nodes. Set nVec to the number of vectors needed on this
      ! node (process).
      If (Cho_Real_Par) Then
         l=nVec(1)
         Do iSym=2,mSym
            l=max(l,nVec(iSym))
         End Do
         Call GetMem('BkmScr','Allo','Inte',ip,l)
         Do iSym=1,mSym
            Call Cho_P_Distrib_Vec(1,nVec(iSym),iWork(ip),n)
            nVec(iSym)=n
         End Do
         Call GetMem('BkmScr','Free','Inte',ip,l)
      End If

#if defined (_DEBUGPRINT_)
      ! Self test
      Do iSym=1,mSym
         If (nVec(iSym).gt.NumCho(iSym)) Then
            Call Cho_Quit('Error in '//SecNam,104)
         End If
      End Do
#endif

      End
