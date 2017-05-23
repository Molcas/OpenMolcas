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
      Subroutine Cho_TraS(iSymL, iSym,jSym, NumV, CMO,NCMO,
     &                               lUCHFV, iStrtVec_AB, nFVec,nFBatch)
************************************************************************
*  This is the routine for the transformation from AO basis to MO      *
*  basis of the Cholesky Full Vectors when  iSym.EQ.jSym.              *
*  The new Transformed Cholesky Full Vectors L are :                   *
*                                                                      *
*   TCVA: L_ij  if DoTCVA=.True.                                       *
*   TCVB: L_tj  if DoTCVA=.True.                                       *
*   TCVC: L_aj                                                         *
*   TCVD: L_tu  if DoTCVA=.True.                                       *
*   TCVE: L_au  if DoTCVA=.True.                                       *
*   TCVF: L_ab  if DoTCVA=.True.                                       *
*   TCVG: L_jt  if DoTCVA=.True.                                       *
*  For generation of <pk|ql>  p,q: All MO, k,l: Occupied (i & t)       *
*  MO Indices  i,j: Inactive;   t,u: Active;   a,b: Secondary          *
*                                                                      *
*----------------------------------------------------------------------*
*  Author  :  Giovanni Ghigo                                           *
*             Lund University, Sweden & University di Torino, Italia   *
*  Written :  October 2004                                             *
*  Modified:  July 2005                                                *
************************************************************************
*
*   <DOC>
*     <Name>Cho\_TraS</Name>
*     <Syntax>Call Cho\_TraS(iSymL,iSym,jSym,NumV,CMO,NCMO,lUCHFV,
*      iStrtVec\_AB,  nFVec,nFBatch)
*     </Syntax>
*     <Arguments>
*      \Argument{iSymL}{Symmetry of the Cholesky vector}{Integer}{in}
*      \Argument{iSym}{Symmetry(i) of the Cholesky full vector}{Integer}
*      {in}
*      \Argument{jSym}{Symmetry(j) of the Cholesky full vector}{Integer}
*      {in}
*      \Argument{NumV}{Number of Cholesky vectors to transform in the
*      current batch}{Integer}{in}
*      \Argument{CMO}{MO coefficients}{Array Real*8}{in}
*      \Argument{NCMO}{Total number of MO coefficients}{Integer}{in}
*      \Argument{lUCHFV}{Unit number of the Cholesky full vector to
*      transform (CHFV}{Integer}{in}
*      \Argument{iStrtVec\_AB}{Current initial disk pointer of the
*      Cholesky full vector to transform (CHFV}{Integer}{in/out}
*      \Argument{nFVec}{Number of Cholesky vectors to transform in the
*      inner batch procedure (nFBatch)}{Integer}{in}
*      \Argument{nFBatch}{Number of cycles in the inner batch procedure}
*     {Integer}{in}
*     </Arguments>
*     <Purpose>
*      Routine for the transformation of the Cholesky vectors in
*      MO-based TCVx for case Sym(i).EQ.Sym(j).\\
*      Called by Cho\_TraCtl.
*     </Purpose>
*     <Dependencies>
*      The logical matrix TCVXist must be defined.
*     </Dependencies>
*     <Author>
*      G. Ghigo
*     </Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*      In the inner batch the Cholesky full vectors are transformed and
*      stored in memory. Adresses (1) and length (2) are stored in the
*      integer matrix iMemTCVX(iType, Sym(i),  Sym(j),1/2).
*      iType: 1=TCVA, 2=TCVB, 3=TCVC, 4=TCVD, 5=TCVE, 6=TCVF, 7=TCVG.
*      Types 1, 2 and 4-7 are generated only if DoTCVA=.True.
*      TCVC is always generated.\\
*      In the first half-transformation the vectors are contracted
*      only with the occupied (inactive and active) MO coefficients
*      for Sym(j). In the second half-transformation the vectors are
*      contracted with all MO coefficients.
*     </Description>
*    </DOC>
*
******************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)

#include "rasdim.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"
      Dimension CMO(NCMO)
      Logical TCVA,TCVB,TCVC,TCVD,TCVE,TCVF

* --- Memory to allocate & Nr. of Cholesky vectors transformable
*     A=Alpha(AO);  B=Beta(AO)

      TCVA = .False.
      TCVB = .False.
      TCVC = .False.
      TCVD = .False.
      TCVE = .False.
      TCVF = .False.

      Nij = 0
      Ntj = 0
      Njt = 0
      Naj = 0
      Ntu = 0
      Nau = 0
      Nab = 0

      Len_FAB = 0
      Len_ij = 0
      Len_tj = 0
      Len_jt = 0
      Len_aj = 0
      Len_tu = 0
      Len_au = 0
      Len_ab = 0

      Len_XAj = 0
      Len_XAu = 0
      Len_XAb = 0

      NFAB = nBas(iSym) * ( nBas(jSym) + 1 ) /2
      Len_ABSq = nBas(iSym) * nBas(jSym)

* TCV-A :
      If (TCVXist(1,iSym,jSym)) Then
        TCVA = .True.
        Len_XAj = nBas(iSym) * nIsh(jSym)
        Nij = nIsh(iSym) * nIsh(jSym)
      EndIf

* TCV-B :
      If (TCVXist(2,iSym,jSym)) Then
        TCVB = .True.
        Len_XAj = nBas(iSym) * nIsh(jSym)
        Ntj = nAsh(iSym) * nIsh(jSym)
        Njt = nIsh(jSym) * nAsh(iSym)
      EndIf

* TCV-C :
      If (TCVXist(3,iSym,jSym)) Then
        TCVC = .True.
        Len_XAj = nBas(iSym) * nIsh(jSym)
        Naj = nSsh(iSym) * nIsh(jSym)
      EndIf

* TCV-D :
      If (TCVXist(4,iSym,jSym)) Then
        TCVD = .True.
        Len_XAu = nBas(iSym) * nAsh(jSym)
        Ntu = nAsh(iSym) * nAsh(jSym)
      EndIf

* TCV-E :
      If (TCVXist(5,iSym,jSym)) Then
        TCVE = .True.
        Len_XAu = nBas(iSym) * nAsh(jSym)
        Nau = nSsh(iSym) * nAsh(jSym)
      EndIf

* TCV-F :
      If (TCVXist(6,iSym,jSym)) Then
        TCVF = .True.
        Len_XAb = nBas(iSym) * nSsh(jSym)
        Nab = nSsh(iSym) * nssh(jSym)
      EndIf

*     Allocate memory for Transformed Cholesky Vectors - TCVx
      Len_ij = Nij * NumV
      Len_tj = Ntj * NumV
      Len_jt = Njt * NumV
      Len_aj = Naj * NumV
      Len_tu = Ntu * NumV
      Len_au = Nau * NumV
      Len_ab = Nab * NumV
      iStrt_ij = 0
      iStrt_tj = 0
      iStrt_jt = 0
      iStrt_aj = 0
      iStrt_tu = 0
      iStrt_au = 0
      iStrt_ab = 0
      iStrt0_ij = 0
      iStrt0_tj = 0
      iStrt0_jt = 0
      iStrt0_aj = 0
      iStrt0_tu = 0
      iStrt0_au = 0
      iStrt0_ab = 0

      If (TCVA) then
        Call GetMem('ij','ALLO','REAL',iStrt00_ij,Len_ij)
        iMemTCVX(1,iSym,jSym,1)=iStrt00_ij
        iMemTCVX(1,iSym,jSym,2)=Len_ij
      EndIf
      If (TCVB) then
        Call GetMem('tj','ALLO','REAL',iStrt00_tj,Len_tj)
        iMemTCVX(2,iSym,jSym,1)=iStrt00_tj
        iMemTCVX(2,iSym,jSym,2)=Len_tj
        Call GetMem('jt','ALLO','REAL',iStrt00_jt,Len_jt)
        iMemTCVX(7,jSym,iSym,1)=iStrt00_jt
        iMemTCVX(7,jSym,iSym,2)=Len_jt
      EndIf
      If (TCVC) then
        Call GetMem('aj','ALLO','REAL',iStrt00_aj,Len_aj)
        iMemTCVX(3,iSym,jSym,1)=iStrt00_aj
        iMemTCVX(3,iSym,jSym,2)=Len_aj
      EndIf
      If (TCVD) then
        Call GetMem('tu','ALLO','REAL',iStrt00_tu,Len_tu)
        iMemTCVX(4,iSym,jSym,1)=iStrt00_tu
        iMemTCVX(4,iSym,jSym,2)=Len_tu
      EndIf
      If (TCVE) then
        Call GetMem('au','ALLO','REAL',iStrt00_au,Len_au)
        iMemTCVX(5,iSym,jSym,1)=iStrt00_au
        iMemTCVX(5,iSym,jSym,2)=Len_au
      EndIf
      If (TCVF) then
        Call GetMem('ab','ALLO','REAL',iStrt00_ab,Len_ab)
        iMemTCVX(6,iSym,jSym,1)=iStrt00_ab
        iMemTCVX(6,iSym,jSym,2)=Len_ab
      EndIf

* --- START LOOP iFBatch -----------------------------------------------
      DO iFBatch=1,nFBatch
        If (iFBatch.EQ.nFBatch) then
         NumFV = NumV - nFVec * (nFBatch-1)
        Else
         NumFV = nFVec
        EndIf
        If ( TCVA ) iStrt0_ij = iStrt00_ij + (iFBatch-1) * nFVec * Nij
        If ( TCVB ) iStrt0_tj = iStrt00_tj + (iFBatch-1) * nFVec * Ntj
        If ( TCVB ) iStrt0_jt = iStrt00_jt + (iFBatch-1) * nFVec * Njt
        If ( TCVC ) iStrt0_aj = iStrt00_aj + (iFBatch-1) * nFVec * Naj
        If ( TCVD ) iStrt0_tu = iStrt00_tu + (iFBatch-1) * nFVec * Ntu
        If ( TCVE ) iStrt0_au = iStrt00_au + (iFBatch-1) * nFVec * Nau
        If ( TCVF ) iStrt0_ab = iStrt00_ab + (iFBatch-1) * nFVec * Nab

*       Allocate memory & Load Full Cholesky Vectors - CHFV
        Len_FAB = NFAB * NumFV
        iStrtVec_FAB = iStrtVec_AB + nFVec * (iFBatch-1)
        Call GetMem('FAB','Allo','Real',iStrt0_FAB,Len_FAB)
        Call RdChoVec(Work(iStrt0_FAB),NFAB,NumFV,iStrtVec_FAB,lUCHFV)

*  ---  Start Loop  iVec  ---
        Do iVec=1,NumFV   ! Loop  iVec
          iStrt_FAB = iStrt0_FAB + (iVec-1) * NFAB
          If ( TCVA ) iStrt_ij = iStrt0_ij + (iVec-1) * Nij
          If ( TCVB ) iStrt_tj = iStrt0_tj + (iVec-1) * Ntj
          If ( TCVB ) iStrt_jt = iStrt0_jt + (iVec-1) * Njt
          If ( TCVC ) iStrt_aj = iStrt0_aj + (iVec-1) * Naj
          If ( TCVD ) iStrt_tu = iStrt0_tu + (iVec-1) * Ntu
          If ( TCVE ) iStrt_au = iStrt0_au + (iVec-1) * Nau
          If ( TCVF ) iStrt_ab = iStrt0_ab + (iVec-1) * Nab

*     --- 1st Half-Transformation  iBeta(AO) -> q(MO) only occupied
          jStrt0MO = 1
          Do j=1,jSym-1
            jStrt0MO = jStrt0MO + nBas(j) * nBas(j)
          EndDo
          jStrt0MO = jStrt0MO + nFro(jSym) * nBas(jSym)

C         From CHFV A(Alpha,Beta) to XAj(Alpha,jMO)
          If ( TCVA .or. TCVB .or. TCVC ) then
            Call GetMem('XAj','ALLO','REAL',iStrt0_XAj,Len_XAj)
            Call ProdsS_1(Work(iStrt_FAB), nBas(iSym),
     &                CMO(jStrt0MO),nIsh(jSym), Work(iStrt0_XAj))
          EndIf
          jStrt0MO = jStrt0MO + nIsh(jSym) * nBas(jSym)

C         From CHFV A(Alpha,Beta) to XAu(Alpha,uMO)
          If ( TCVD .or. TCVE ) then
            Call GetMem('XAu','ALLO','REAL',iStrt0_XAu,Len_XAu)
            Call ProdsS_1(Work(iStrt_FAB), nBas(iSym),
     &                CMO(jStrt0MO),nAsh(jSym), Work(iStrt0_XAu))
          EndIf
          jStrt0MO = jStrt0MO + nAsh(jSym) * nBas(jSym)

C         From CHFV A(Alpha,Beta) to XAb(Alpha,bMO)
          If ( TCVF ) then
            Call GetMem('XAb','ALLO','REAL',iStrt0_XAb,Len_XAb)
            Call ProdsS_1(Work(iStrt_FAB), nBas(iSym),
     &                CMO(jStrt0MO),nSsh(jSym), Work(iStrt0_XAb))
          EndIf

*     --- 2nd Half-Transformation  iAlpha(AO) -> p(MO)
          iStrt0MO = 1
          Do i=1,iSym-1
            iStrt0MO = iStrt0MO + nBas(i) * nBas(i)
          EndDo
          iStrt0MO = iStrt0MO + nFro(iSym) * nBas(iSym)

C         From XAj(Alpha,jMO) to ij(i,j)
          If ( TCVA ) then
            Call ProdsS_2(Work(iStrt0_XAj), nBas(iSym),nIsh(jSym),
     &                CMO(iStrt0MO),nIsh(iSym), Work(iStrt_ij))
          EndIf
          iStrt0MO = iStrt0MO + nIsh(iSym) * nBas(iSym)

C         From XAj(Alpha,jMO) to tj(t,j) & jt(j,t)
          If ( TCVB ) then
            Call ProdsS_2(Work(iStrt0_XAj), nBas(iSym),nIsh(jSym),
     &                CMO(iStrt0MO),nAsh(iSym), Work(iStrt_tj))
            Call Trnsps(nAsh(iSym),nIsh(jSym),
     &                                    Work(iStrt_tj),Work(iStrt_jt))
          EndIf

C         From XAu(Alpha,uMO) to tu(t,u)
          If ( TCVD ) then
            Call ProdsS_2(Work(iStrt0_XAu), nBas(iSym),nAsh(jSym),
     &                CMO(iStrt0MO),nAsh(iSym), Work(iStrt_tu))
          EndIf
          iStrt0MO = iStrt0MO + nAsh(iSym) * nBas(iSym)

C         From XAj(Alpha,jMO) to aj(a,j)
          If ( TCVC ) then
            Call ProdsS_2(Work(iStrt0_XAj), nBas(iSym),nIsh(jSym),
     &                CMO(iStrt0MO),nSsh(iSym), Work(iStrt_aj))
          EndIf

C         From XAu(Alpha,uMO) to au(a,u)
          If ( TCVE ) then
            Call ProdsS_2(Work(iStrt0_XAu), nBas(iSym),nAsh(jSym),
     &                CMO(iStrt0MO),nSsh(iSym), Work(iStrt_au))
          EndIf

C         From XAb(Alpha,bMO) to ab(a,b)
          If ( TCVF ) then
            Call ProdsS_2(Work(iStrt0_XAb), nBas(iSym),nSsh(jSym),
     &                CMO(iStrt0MO),nSsh(iSym), Work(iStrt_ab))
          EndIf

*     --- End of Transformations

          If ( TCVA .or. TCVB .or. TCVC ) then
            Call GetMem('XAj','FREE','REAL',iStrt0_XAj,Len_XAj)
          EndIf
          If ( TCVD .or. TCVE ) then
            Call GetMem('XAu','FREE','REAL',iStrt0_XAu,Len_XAu)
          EndIf
          If ( TCVF ) then
            Call GetMem('XAb','FREE','REAL',iStrt0_XAb,Len_XAb)
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

      Subroutine ProdsS_1(AB,iAB, CMO,nMO, Y)
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "WrkSpc.fh"
      Dimension AB(iAB*(iAB+1)/2), CMO(iAB,nMO), Y(iAB,nMO)
      Call GetMem('ABSq','ALLO','REAL',iAddABSq,iAB*iAB)
      Call SQUARE(AB,Work(iAddABSq),1,iAB,iAB)
      Call DGEMM_('N','N',iAB,nMO,iAB,1.0d0,Work(iAddABSq),iAB,CMO,iAB,
     &                                                    0.0d0,Y,iAB)
      Call GetMem('ABSq','FREE','REAL',iAddABSq,iAB*2)
      Return
      End

      Subroutine ProdsS_2(AB,iA,jMO, CMO,nMO, Y)
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Dimension AB(iA,jMO), CMO(iA,nMO), Y(nMO,jMO)
      Call DGEMM_('T','N',nMO,jMO,iA,1.0d0,CMO,iA,AB,iA,0.0d0,Y,nMO)
      Return
      End
