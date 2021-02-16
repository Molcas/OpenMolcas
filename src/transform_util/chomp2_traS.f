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
* Copyright (C) 2004,2005, Giovanni Ghigo                              *
************************************************************************
      Subroutine ChoMP2_TraS(iSymL, iSym,jSym, NumV, CMO,NCMO,
     &                               lUCHFV, iStrtVec_AB, nFVec,nFBatch)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden                                    *
* Written:  October 2004                                               *
* Modified for Cholesky-MP2 May 2005                                   *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Integer NCMO
      Real*8 CMO(NCMO)

#include "rasdim.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"

      Real*8, Allocatable:: XAj(:), FAB(:,:)

* --- Memory to allocate & Nr. of Cholesky vectors transformable
*     A=Alpha(AO);  B=Beta(AO)

      IF ( .NOT. TCVXist(3,iSym,jSym)) RETURN

      Naj = 0
      Len_aj = 0
      Len_XAj = 0
      NFAB = nBas(iSym) * ( nBas(jSym) + 1 ) /2
      Len_XAj = nBas(iSym) * nIsh(jSym)
      Naj = nSsh(iSym) * nIsh(jSym)

*     Allocate memory for Transformed Cholesky Vectors - TCVx
      Len_aj = Naj * NumV
      iStrt_aj = 0
      Call GetMem('aj','ALLO','REAL',iStrt00_aj,Len_aj)
      iMemTCVX(3,iSym,jSym,1)=iStrt00_aj
      iMemTCVX(3,iSym,jSym,2)=Len_aj

      iStrt = 1
      Do i=1,iSym-1
         iStrt = iStrt + nBas(i) * nBas(i)
      EndDo
      jStrt = 1
      Do j=1,jSym-1
         jStrt = jStrt + nBas(j) * nBas(j)
      EndDo

* --- START LOOP iiVec   -----------------------------------------------
      DO iiVec = 1, NumV, nFVec
        NumFV=Max(nFVec,NumV-iiVec+1)
        iFBatch = (iiVec+nFVec-1)/nFVec

*       Allocate memory & Load Full Cholesky Vectors - CHFV

        iStrtVec_FAB = iStrtVec_AB + nFVec * (iFBatch-1)

        Call mma_allocate(FAB,NFAB,NumFV,Label='FAB')
        Call RdChoVec(FAB,NFAB,NumFV,iStrtVec_FAB,lUCHFV)

*  ---  Start Loop jVec  ---
        Do jVec=iiVec,iiVec+NumFV-1   ! Loop  jVec
          iVec = jVec - iiVec + 1

          iStrt_aj = iStrt00_aj + (jVec-1) * Naj

*     --- 1st Half-Transformation  iBeta(AO) -> q(MO) only occupied
          jStrt0MO = jStrt + nFro(jSym) * nBas(jSym)

C         From CHFV A(Alpha,Beta) to XAj(Alpha,jMO)
          Call mma_allocate(XAj,Len_XAj,Label='XAj')
          Call ProdsS_1(FAB(:,iVec), nBas(iSym),
     &                CMO(jStrt0MO),nIsh(jSym), XAj)

*     --- 2nd Half-Transformation  iAlpha(AO) -> p(MO)
          iStrt0MO = iStrt + (nFro(iSym)+nIsh(iSym)) * nBas(iSym)

C         From XAj(Alpha,jMO) to aj(a,j)
          Call ProdsS_2(XAj, nBas(iSym),nIsh(jSym),
     &              CMO(iStrt0MO),nSsh(iSym), Work(iStrt_aj))

*     --- End of Transformations

          Call mma_deallocate(XAj)

        EndDo
*  ---  End Loop  jVec  ---

        Call mma_deallocate(FAB)
      ENDDO
* --- END LOOP iiVec   -------------------------------------------------

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iSymL)
      End
