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
      Subroutine ChoMP2_TraA(iSymL, iSym,jSym, NumV, CMO,NCMO,
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
      Logical TCVC,TCVCt

* --- Memory to allocate & Nr. of Cholesky vectors transformable
*     A=Alpha(AO);  B=Beta(AO)
      TCVC = .False.
      TCVCt= .False.

      NFAB = 0
      Naj = 0  ! C
      Nbi = 0  ! C"
      Len_FAB = 0
      Len_aj = 0  ! C
      Len_bi = 0  ! C"
      Len_XAj = 0   ! C
      Len_XBi = 0   ! C"
      iStrt_aj = 0  ! C
      iStrt_bi = 0  ! C"
      iStrt0_aj = 0  ! C
      iStrt0_bi = 0  ! C"
      NFAB = nBas(iSym) * nBas(jSym)

*     Allocate memory for Transformed Cholesky Vectors - TCVx
* TCV-C :
      If (TCVXist(3,iSym,jSym)) Then
        TCVC = .True.
        Len_XAj = nBas(iSym) * nIsh(jSym)
        Naj = nSsh(iSym) * nIsh(jSym)
        Len_aj = Naj * NumV
        Call GetMem('aj','ALLO','REAL',iStrt00_aj,Len_aj)
        iMemTCVX(3,iSym,jSym,1)=iStrt00_aj
        iMemTCVX(3,iSym,jSym,2)=Len_aj
      EndIf
* TCV-Ct:
      If (TCVXist(3,jSym,iSym)) Then
        TCVCt= .True.
        Len_XBi = nBas(jSym) * nIsh(iSym)
        Nbi = nSsh(jSym) * nIsh(iSym)
        Len_bi = Nbi * NumV
        Call GetMem('bi','ALLO','REAL',iStrt00_bi,Len_bi)
        iMemTCVX(3,jSym,iSym,1)=iStrt00_bi
        iMemTCVX(3,jSym,iSym,2)=Len_bi
      EndIf

* --- START LOOP iFBatch -----------------------------------------------
      DO iFBatch=1,nFBatch
        If (iFBatch.EQ.nFBatch) then
         NumFV = NumV - nFVec * (nFBatch-1)
        Else
         NumFV = nFVec
        EndIf
        If ( TCVC  ) iStrt0_aj = iStrt00_aj + (iFBatch-1) * nFVec * Naj
        If ( TCVCt ) iStrt0_bi = iStrt00_bi + (iFBatch-1) * nFVec * Nbi

*       Allocate memory & Load Full Cholesky Vectors - CHFV
        Len_FAB = NFAB * NumFV
        iStrtVec_FAB = iStrtVec_AB + nFVec * (iFBatch-1)
        Call GetMem('FAB','Allo','Real',iStrt0_FAB,Len_FAB)
        Call RdChoVec(Work(iStrt0_FAB),NFAB,NumFV,iStrtVec_FAB,lUCHFV)

*  ---  Start Loop  iVec  ---
        Do iVec=1,NumFV   ! Loop  iVec
          iStrt_FAB = iStrt0_FAB + (iVec-1) * NFAB
          If ( TCVC  ) iStrt_aj = iStrt0_aj + (iVec-1) * Naj
          If ( TCVCt ) iStrt_bi = iStrt0_bi + (iVec-1) * Nbi

*     --- 1st Half-Transformation  iBeta(AO) -> q(MO) only occupied
C         From CHFV A(Alpha,Beta) to XAj(Alpha,jMO)
          If ( TCVC ) then
            jStrt0MO = 1
            Do j=1,jSym-1
              jStrt0MO = jStrt0MO + nBas(j) * nBas(j)
            EndDo
            jStrt0MO = jStrt0MO + nFro(jSym) * nBas(jSym)
            Call GetMem('XAj','ALLO','REAL',iStrt0_XAj,Len_XAj)
            Call ProdsA_2(Work(iStrt_FAB), nBas(iSym),nBas(jSym),
     &                  CMO(jStrt0MO),nIsh(jSym), Work(iStrt0_XAj))
          EndIf
C         From CHFV A(Alpha,Beta) to XBi(Beta,iMO)
          If ( TCVCt ) then
            iStrt0MO = 1
            Do i=1,iSym-1
              iStrt0MO = iStrt0MO + nBas(i) * nBas(i)
            EndDo
            iStrt0MO = iStrt0MO + nFro(iSym) * nBas(iSym)
            Call GetMem('XBi','ALLO','REAL',iStrt0_XBi,Len_XBi)
            Call ProdsA_2t(Work(iStrt_FAB), nBas(iSym),nBas(jSym),
     &                  CMO(iStrt0MO),nIsh(iSym), Work(iStrt0_XBi))
          EndIf

*     --- 2nd Half-Transformation  iAlpha(AO) -> p(MO)
C         From XAj(Alpha,jMO) to aj(a,j)
          If ( TCVC ) then
            iStrt0MO = 1
            Do i=1,iSym-1
              iStrt0MO = iStrt0MO + nBas(i) * nBas(i)
            EndDo
            iStrt0MO = iStrt0MO + (nFro(iSym)+nIsh(iSym)) * nBas(iSym)
            Call ProdsA_1(Work(iStrt0_XAj), nBas(iSym),nIsh(jSym),
     &                  CMO(iStrt0MO),nSsh(iSym), Work(iStrt_aj))
          EndIf

C         From XBi(Beta,jMO) to bi(b,i)
          If ( TCVCt ) then
            jStrt0MO = 1
            Do j=1,jSym-1
              jStrt0MO = jStrt0MO + nBas(j) * nBas(j)
            EndDo
            jStrt0MO = jStrt0MO + (nFro(jSym)+nIsh(jSym)) * nBas(jSym)
            Call ProdsA_1(Work(iStrt0_XBi), nBas(jSym),nIsh(iSym),
     &                  CMO(jStrt0MO),nSsh(jSym), Work(iStrt_bi))
          EndIf

*     --- End of Transformations

          If ( TCVC ) then
            Call GetMem('XAj','FREE','REAL',iStrt0_XAj,Len_XAj)
          EndIf
          If ( TCVCt ) then
            Call GetMem('XBi','FREE','REAL',iStrt0_XBi,Len_XBi)
          EndIf

        EndDo
*  ---  End Loop  iVec  ---

        Call GetMem('FAB','Free','Real',iStrt0_FAB,Len_FAB)
      ENDDO
* --- END LOOP iFBatch -------------------------------------------------

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iSymL)
      End
