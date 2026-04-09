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
      use definitions, only: iwp, wp, u6
      use TList_Mod, only: iTskCan, Not_Used, TskQ, QLast
      use SysDef, only: RtoI
      Implicit None
      real(kind=wp), intent(inout):: Disc

      Integer(kind=iwp) iWR(2), mInts
      Real(kind=wp) Dummy(1), Quad_ijkl, RST_triplet
      Logical :: Copy=.True., NoCopy=.False.
*
      if(.NOT.Allocated(TskQ)) return

      Quad_ijkl  =TskQ(1,iTskCan)
      RST_triplet=TskQ(2,iTskCan)
      If (Quad_ijkl.eq.Not_Used) Return
*
*---- If already at the right position return
*
      If (Quad_ijkl==QLast(1) .and.
     &    RST_triplet==QLast(2)) Return
*
      Do

         Call iRBuf(iWR,2,Copy)
         Call dRBuf(QLast,2,Copy)
         mInts=iWR(2)
         If (QLast(1)==Quad_ijkl .and.
     &       QLast(2)==RST_triplet) Then
            If (mInts.gt.0) Call dRBuf(Dummy,mInts,NoCopy)
            Disc = Disc + DBLE(2/RtoI + 2 + mInts)
            Return
         Else If (QLast(1)<=Quad_ijkl ) Then
            If (mInts>0) Call dRBuf(Dummy,mInts,NoCopy)
            Disc = Disc + DBLE(2/RtoI + 2 + mInts)
            Cycle
         Else
            Write (u6,*) 'Pos_QLast: batch is lost!'
            Write (u6,'(A,2F10.1)') 'Index,1.0:  ',QLast(1),QLast(2)
            Write (u6,'(A,2F10.1)') 'Looking for ',Quad_ijkl,RST_triplet
            Write (u6,*) ' iTskCan,=',iTskCan
            Call RecPrt('TskQ',' ',TskQ,2,iTskCan)
            Write (u6,*)
            Call Abend()
       End If
*
      End Do

      Write (u6,*) 'Pos_QLast: Fatal problem!'
      Call XFlush(6)
      Call Abend()

      End Subroutine Pos_QLast
