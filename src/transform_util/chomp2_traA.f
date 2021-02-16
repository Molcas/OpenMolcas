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
      Integer NCMO
      Real*8 CMO(NCMO)
#include "rasdim.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"
      Logical TCVC,TCVCt

      Real*8, Allocatable:: XAj(:), XBi(:), FAB(:,:)

* --- Memory to allocate & Nr. of Cholesky vectors transformable
*     A=Alpha(AO);  B=Beta(AO)
      TCVC = .False.
      TCVCt= .False.

      NFAB = 0
      Naj = 0  ! C
      Nbi = 0  ! C"
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

        Call mma_allocate(FAB,nFAB,NumFV,Label='FAB')
        Call RdChoVec(FAB,NFAB,NumFV,iStrtVec_FAB,lUCHFV)

*  ---  Start Loop  jVec  ---
        Do jVec=iiVec,iiVec+NumFV-1   ! Loop  jVec
          iVec = jVec - iiVec + 1

          If ( TCVC  ) iStrt_aj = iStrt00_aj + (jVec-1) * Naj
          If ( TCVCt ) iStrt_bi = iStrt00_bi + (jVec-1) * Nbi

*     --- 1st Half-Transformation  iBeta(AO) -> q(MO) only occupied
C         From CHFV A(Alpha,Beta) to XAj(Alpha,jMO)
          If ( TCVC ) then
            jStrt0MO = jStrt + nFro(jSym) * nBas(jSym)
            Call mma_allocate(XAj,Len_XAj,Label='XAj')
            Call ProdsA_2(FAB(:,iVec), nBas(iSym),nBas(jSym),
     &                  CMO(jStrt0MO),nIsh(jSym), XAj)
          EndIf
C         From CHFV A(Alpha,Beta) to XBi(Beta,iMO)
          If ( TCVCt ) then
            iStrt0MO = iStrt + nFro(iSym) * nBas(iSym)
            Call mma_allocate(XBi,Len_XBi,Label='XBi')
            Call ProdsA_2t(FAB(:,iVec), nBas(iSym),nBas(jSym),
     &                  CMO(iStrt0MO),nIsh(iSym), XBi)
          EndIf

*     --- 2nd Half-Transformation  iAlpha(AO) -> p(MO)
C         From XAj(Alpha,jMO) to aj(a,j)
          If ( TCVC ) then
            iStrt0MO = iStrt + (nFro(iSym)+nIsh(iSym)) * nBas(iSym)
            Call ProdsA_1(XAj, nBas(iSym),nIsh(jSym),
     &                  CMO(iStrt0MO),nSsh(iSym), Work(iStrt_aj))
          EndIf

C         From XBi(Beta,jMO) to bi(b,i)
          If ( TCVCt ) then
            jStrt0MO = jStrt + (nFro(jSym)+nIsh(jSym)) * nBas(jSym)
            Call ProdsA_1(XBi, nBas(jSym),nIsh(iSym),
     &                  CMO(jStrt0MO),nSsh(jSym), Work(iStrt_bi))
          EndIf

*     --- End of Transformations

          If (Allocated(XAj)) Call mma_deallocate(XAj)
          If (Allocated(XBi)) Call mma_deallocate(XBi)

        EndDo
*  ---  End Loop  jVec  ---

        Call mma_deallocate(FAB)
      ENDDO
* --- END LOOP iiVec   -------------------------------------------------

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iSymL)
      End
