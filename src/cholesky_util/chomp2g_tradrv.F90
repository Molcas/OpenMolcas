!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2010, Jonas Bostrom                                    *
!***********************************************************************
      SubRoutine ChoMP2g_TraDrv(irc,CMO,Diag,DoDiag)
!
!     Jonas Bostrom, Jan. 2010. (modified from ChoMP2_TraDrv)
!
!     Purpose: AO-to-MO (pq) transformation of Cholesky vectors
!              performed directly in reduced sets. This assumes
!              that the MP2 program has been appropriately initialized.
!
      use stdalloc
      use ChoMP2g
      Implicit Real*8 (a-h,o-z)
      Integer irc
      Real*8  CMO(*), Diag(*)
      Logical DoDiag, DoDiagbak
#include "cholesky.fh"
#include "chomp2.fh"

      Character(LEN=6), Parameter:: ThisNm = 'TraDrv'
      Character(LEN=14), Parameter:: SecNam = 'ChoMP2g_TraDrv'

      Real*8, Allocatable:: COrb1(:), COrb2(:)

      irc = 0

!     Reorder MO coefficients.
!     ------------------------

      DoDiagBak = DoDiag
      DoDiag = .false.
      nProdType = nMOType**2
      l_COrb = 0
      Do iSym = 1, nSym
         nAdrOff(iSym) = 0
      End Do
!
      Do iSym = 1, nSym
         Do i = 1, nProdType
            l_COrb = max(l_COrb,nMoAo(iSym,i))
         End Do
      End Do
      Call mma_allocate(COrb1,l_COrb,Label='COrb1')
      Call mma_allocate(COrb2,l_COrb,Label='COrb2')
!
      DoDiag = .True.
      Call ChoMP2g_MOReOrd(CMO,COrb1,COrb2,2,3)
      Call ChoMP2g_Tra(COrb1,COrb2,Diag,DoDiag,2,3)
      DoDiag = .False.
      Do iMoType = 1,3
         Do jMOType =  1, 3
            If((iMoType .eq. 2) .and. (jMoType .eq. 3)) Go To 50

            Call ChoMP2g_MOReOrd(CMO,COrb1,COrb2,iMOType,jMOType)
!           Transform vectors.
!           ------------------
            Call ChoMP2g_Tra(COrb1,COrb2,Diag,DoDiag,iMoType,jMoType)

 50         Continue
         End Do
      End Do

!     Deallocate reordered MO coefficients.
!     -------------------------------------

      DoDiag = DoDiagBak

      Call mma_deallocate(COrb2)
      Call mma_deallocate(COrb1)

      End
