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
      SUBROUTINE Cho_GetMQ(MQ,l_MQ,List_QShp,nQShp)
      use ChoSwp, only: iQuAB, nnBstRSh, iiBstRSh, IndRSh, IndRed
      Implicit Real*8 (a-h,o-z)

      Integer  l_MQ, nQShp
      Real*8   MQ(l_MQ)
      Integer  List_QShp(nQShp)

#include "cholesky.fh"
#include "choprint.fh"
#include "stdalloc.fh"

      Character(LEN=9), Parameter:: SecNam = 'Cho_GetMQ'
      Integer, Parameter:: iOpt = 2

      Integer, External:: Cho_P_LocalSP, Cho_F2SP

      Real*8, Allocatable:: Scr(:)
      Integer, Allocatable:: kOff_Shp(:)

      nTot = nQual(1)
      Do iSym = 2,nSym
         nTot = nTot + nQual(iSym)
      End Do
      If (nTot .lt. 1) Return ! this test makes sense for parallel run

      Call mma_allocate(kOff_Shp,nnShl,Label='kOff_Shp')

      iQoff=0
      Do jSym=1,nSym

         If (nQual(jSym) .lt. 1) goto 10 ! next symmetry

         Lint = 0
         Do iShp=1,nQShp  ! set only the needed offsets
            iL_ShpG = List_QShp(iShp) ! Shell pair
            iL_Shp = Cho_P_LocalSP(iL_ShpG) ! local shell pair
            kOff_Shp(iL_Shp) = Lint
            Lint = Lint + nnBstRSh(jSym,iL_Shp,2)
         End Do

         Call mma_allocate(Scr,Lint,Label='Scr')

C --- Read the integrals
C ----------------------
         Do jQ=1,nQual(jSym)

            iAdr = nnBstr(jSym,2)*(jQ-1)

            Do iShp=1,nQShp

               iL_ShpG = List_QShp(iShp)
               iL_Shp = Cho_P_LocalSP(iL_ShpG) ! local shell pair
               ipS = 1 + kOff_Shp(iL_Shp)
               iShpAdr = iAdr + iiBstRSh(jSym,iL_Shp,2)
               Lread = nnBstRSh(Jsym,iL_Shp,2)
               Call dDaFile(LuSel(JSym),iOpt,Scr(ipS),Lread,iShpAdr)

            End Do

C --- Extract the matrix of qualified integrals
C ---------------------------------------------
            Do iQ=1,nQual(jSym)

               iAB = iQuAB(iQ,jSym)  ! addr curr red set
               isAB = iAB - iiBstR(jSym,2)  ! symm. reduction
               iShABG = IndRsh(IndRed(iAB,2)) ! glob. SP it belongs to
               iShAB = Cho_P_LocalSP(Cho_F2SP(iShABG)) ! local SP
               iabSh = isAB - iiBstRSh(jSym,iShAB,2) ! addr within Shp
               ipfr = kOff_Shp(iShAB) + iabSh
               ipto = iQoff + nQual(jSym)*(jQ-1) + iQ
               MQ(ipto) = Scr(ipfr)

             End Do

         End Do

         iQoff = iQoff + nQual(jSym)**2

         Call mma_deallocate(Scr)

10       Continue

      End Do

      Call mma_deallocate(kOff_Shp)

      Return
      End
