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
      Subroutine Pos_QLast(Disc)
      use TList_Mod
      use SysDef, only: RtoI
      Implicit Real*8 (a-h,o-z)
      Integer iWR(2)
      Real*8 Dummy(1)
      Logical :: Copy=.True., NoCopy=.False.
*
      if(.NOT.Allocated(TskQ)) return

      Quad_ijkl  =TskQ(1,iTskCan)
      RST_triplet=TskQ(2,iTskCan)
      If (Quad_ijkl.eq.Not_Used) Return
*
*---- If already at the right position return
*
      If (Quad_ijkl.eq.QLast(1) .and.
     &    RST_triplet.eq.QLast(2)) Return
*
 1111 Continue
c     Call Diskat
      Call iRBuf(iWR,2,Copy)
      Call dRBuf(QLast,2,Copy)
      mInts=iWR(2)
      If (QLast(1).eq.Quad_ijkl .and.
     &    QLast(2).eq.RST_triplet) Then
         If (mInts.gt.0) Call dRBuf(Dummy,mInts,NoCopy)
         Disc = Disc + DBLE(2/RtoI + 2 + mInts)
         Return
      Else If (QLast(1).le.Quad_ijkl ) Then
         If (mInts.gt.0) Call dRBuf(Dummy,mInts,NoCopy)
         Disc = Disc + DBLE(2/RtoI + 2 + mInts)
         Go To 1111
      Else
         Write (6,*) 'Pos_QLast: batch is lost!'
         Write (6,'(A,2F10.1)') 'Index,1.0:  ',QLast(1),QLast(2)
         Write (6,'(A,2F10.1)') 'Looking for ',Quad_ijkl,RST_triplet
         Write (6,*) ' iTskCan,=',iTskCan
         Call RecPrt('TskQ',' ',TskQ,2,iTskCan)
         Write (6,*)
         Call XFlush(6)
         Call Abend()
      End If
*
      Write (6,*) 'Pos_QLast: Fatal problem!'
      Call XFlush(6)
      Call Abend
      End Subroutine Pos_QLast
