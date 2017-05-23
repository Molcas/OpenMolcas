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

#include "rasdim.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"
      Dimension CMO(NCMO)

* --- Memory to allocate & Nr. of Cholesky vectors transformable
*     A=Alpha(AO);  B=Beta(AO)

      IF ( .NOT. TCVXist(3,iSym,jSym)) RETURN

      Naj = 0
      Len_FAB = 0
      Len_aj = 0
      Len_XAj = 0
      NFAB = nBas(iSym) * ( nBas(jSym) + 1 ) /2
      Len_XAj = nBas(iSym) * nIsh(jSym)
      Naj = nSsh(iSym) * nIsh(jSym)

*     Allocate memory for Transformed Cholesky Vectors - TCVx
      Len_aj = Naj * NumV
      iStrt_aj = 0
      iStrt0_aj = 0
      Call GetMem('aj','ALLO','REAL',iStrt00_aj,Len_aj)
      iMemTCVX(3,iSym,jSym,1)=iStrt00_aj
      iMemTCVX(3,iSym,jSym,2)=Len_aj

* --- START LOOP iFBatch -----------------------------------------------
      DO iFBatch=1,nFBatch
        If (iFBatch.EQ.nFBatch) then
         NumFV = NumV - nFVec * (nFBatch-1)
        Else
         NumFV = nFVec
        EndIf
        iStrt0_aj = iStrt00_aj + (iFBatch-1) * nFVec * Naj

*       Allocate memory & Load Full Cholesky Vectors - CHFV
        Len_FAB = NFAB * NumFV
        iStrtVec_FAB = iStrtVec_AB + nFVec * (iFBatch-1)
        Call GetMem('FAB','Allo','Real',iStrt0_FAB,Len_FAB)
        Call RdChoVec(Work(iStrt0_FAB),NFAB,NumFV,iStrtVec_FAB,lUCHFV)

*  ---  Start Loop  iVec  ---
        Do iVec=1,NumFV   ! Loop  iVec
          iStrt_FAB = iStrt0_FAB + (iVec-1) * NFAB
          iStrt_aj = iStrt0_aj + (iVec-1) * Naj

*     --- 1st Half-Transformation  iBeta(AO) -> q(MO) only occupied
          jStrt0MO = 1
          Do j=1,jSym-1
            jStrt0MO = jStrt0MO + nBas(j) * nBas(j)
          EndDo
          jStrt0MO = jStrt0MO + nFro(jSym) * nBas(jSym)

C         From CHFV A(Alpha,Beta) to XAj(Alpha,jMO)
          Call GetMem('XAj','ALLO','REAL',iStrt0_XAj,Len_XAj)
          Call ProdsS_1(Work(iStrt_FAB), nBas(iSym),
     &                CMO(jStrt0MO),nIsh(jSym), Work(iStrt0_XAj))

*     --- 2nd Half-Transformation  iAlpha(AO) -> p(MO)
          iStrt0MO = 1
          Do i=1,iSym-1
            iStrt0MO = iStrt0MO + nBas(i) * nBas(i)
          EndDo
          iStrt0MO = iStrt0MO + (nFro(iSym)+nIsh(iSym)) * nBas(iSym)

C         From XAj(Alpha,jMO) to aj(a,j)
          Call ProdsS_2(Work(iStrt0_XAj), nBas(iSym),nIsh(jSym),
     &              CMO(iStrt0MO),nSsh(iSym), Work(iStrt_aj))

*     --- End of Transformations

          Call GetMem('XAj','FREE','REAL',iStrt0_XAj,Len_XAj)

        EndDo
*  ---  End Loop  iVec  ---

        Call GetMem('FAB','Free','Real',iStrt0_FAB,Len_FAB)
      ENDDO
* --- END LOOP iFBatch -------------------------------------------------

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iSymL)
      End
