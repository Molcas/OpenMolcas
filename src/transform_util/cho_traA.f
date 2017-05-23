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
      Subroutine Cho_TraA(iSymL, iSym,jSym, NumV, CMO,NCMO,
     &                               lUCHFV, iStrtVec_AB, nFVec,nFBatch)
************************************************************************
*  This is the routine for the transformation from AO basis to MO      *
*  basis of the Cholesky Full Vectors when  iSym.NE.jSym.              *
*  The new Transformed Cholesky Full Vectors L are :                   *
*                                                                      *
*  IF DoFull.EQ.False. (CASPT2, MCLR,MBPT2)                            *
*   TCVA : L_ij, L_ji=T(L_ij) if DoTCVA=.True.                         *
*   TCVB : L_tj, L_ui         if DoTCVA=.True.                         *
*   TCVC : L_aj, L_bi                                                  *
*   TCVD : L_tu, L_ut=T(L_tu) if DoTCVA=.True.                         *
*   TCVE : L_au, L_bt         if DoTCVA=.True.                         *
*   TCVF : L_ab               if DoTCVA=.True.                         *
*   TCVG : L_jt, L_iu         if DoTCVA=.True.                         *
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
*     <Name>Cho\_TraA</Name>
*     <Syntax>Call Cho\_TraA(iSymL,iSym,jSym,NumV,CMO,NCMO,lUCHFV,
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
*      {Integer}{in}
*     </Arguments>
*     <Purpose>
*      Routine for the transformation of the Cholesky vectors in
*      MO-based TCVx for case Sym(i).NE.Sym(j).\\
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
      Logical TCVA,TCVB,TCVBt,TCVC,TCVCt,TCVD,TCVE,TCVEt,TCVF

* --- Memory to allocate & Nr. of Cholesky vectors transformable
*     A=Alpha(AO);  B=Beta(AO)
      TCVA = .False.
      TCVB = .False.
      TCVBt= .False.
      TCVC = .False.
      TCVCt= .False.
      TCVD = .False.
      TCVE = .False.
      TCVEt= .False.
      TCVF = .False.

      NFAB = 0
      Nij = 0  ! A
      Nji = 0  ! A"
      Ntj = 0  ! B
      Njt = 0  ! G
      Nui = 0  ! B"
      Niu = 0  ! G"
      Naj = 0  ! C
      Nbi = 0  ! C"
      Ntu = 0  ! D
      Nut = 0  ! D"
      Nau = 0  ! E
      Nbt = 0  ! E"
      Nab = 0  ! F

      Len_FAB = 0
      Len_ij = 0  ! A
      Len_ji = 0  ! A"
      Len_tj = 0  ! B
      Len_jt = 0  ! G
      Len_ui = 0  ! B"
      Len_iu = 0  ! G"
      Len_aj = 0  ! C
      Len_bi = 0  ! C"
      Len_tu = 0  ! D
      Len_ut = 0  ! D"
      Len_au = 0  ! E
      Len_bt = 0  ! E"
      Len_ab = 0  ! F

      Len_XAj = 0   ! A/A", B, C
      Len_XAu = 0   ! D/D", E
      Len_XAb = 0   ! F
      Len_XBi = 0   ! B",C"
      Len_XBt = 0   ! E"

      NFAB = nBas(iSym) * nBas(jSym)

* TCV-A :
      If (TCVXist(1,iSym,jSym)) Then
        TCVA = .True.
        Len_XAj = nBas(iSym) * nIsh(jSym)
        Nij = nIsh(iSym) * nIsh(jSym)
        Nji = nIsh(jSym) * nIsh(iSym)
      EndIf

* TCV-B :
      If (TCVXist(2,iSym,jSym)) Then
        TCVB = .True.
        Len_XAj = nBas(iSym) * nIsh(jSym)
        Ntj = nAsh(iSym) * nIsh(jSym)
        Njt = nIsh(jSym) * nAsh(iSym)
      EndIf
      If (TCVXist(2,jSym,iSym)) Then
        TCVBt= .True.
        Len_XBi = nBas(jSym) * nIsh(iSym)
        Nui = nAsh(jSym) * nIsh(iSym)
        Niu = nIsh(iSym) * nAsh(jSym)
      EndIf

* TCV-C :
      If (TCVXist(3,iSym,jSym)) Then
        TCVC = .True.
        Len_XAj = nBas(iSym) * nIsh(jSym)
        Naj = nSsh(iSym) * nIsh(jSym)
      EndIf
      If (TCVXist(3,jSym,iSym)) Then
        TCVCt= .True.
        Len_XBi = nBas(jSym) * nIsh(iSym)
        Nbi = nSsh(jSym) * nIsh(iSym)
      EndIf

* TCV-D :
      If (TCVXist(4,iSym,jSym)) Then
        TCVD = .True.
        Len_XAu = nBas(iSym) * nAsh(jSym)
        Ntu = nAsh(iSym) * nAsh(jSym)
        Nut = nAsh(jSym) * nAsh(iSym)
      EndIf

* TCV-E :
      If (TCVXist(5,iSym,jSym)) Then
        TCVE = .True.
        Len_XAu = nBas(iSym) * nAsh(jSym)
        Nau = nSsh(iSym) * nAsh(jSym)
      EndIf
      If (TCVXist(5,jSym,iSym)) Then
        TCVEt= .True.
        Len_XBt = nBas(jSym) * nAsh(iSym)
        Nbt = nSsh(jSym) * nAsh(iSym)
      EndIf

* TCV-F :
      If (TCVXist(6,iSym,jSym)) Then
        TCVF = .True.
        Len_XAb = nBas(iSym) * nSsh(jSym)
        Nab = nSsh(iSym) * nSsh(jSym)
      EndIf

*     Allocate memory for Transformed Cholesky Vectors - TCVx
      Len_ij = Nij * NumV
      Len_ji = Nji * NumV
      Len_tj = Ntj * NumV
      Len_jt = Njt * NumV
      Len_ui = Nui * NumV
      Len_iu = Niu * NumV
      Len_aj = Naj * NumV
      Len_bi = Nbi * NumV
      Len_tu = Ntu * NumV
      Len_ut = Nut * NumV
      Len_au = Nau * NumV
      Len_bt = Nbt * NumV
      Len_ab = Nab * NumV
      iStrt_ij = 0  ! A
      iStrt_ji = 0  ! A"
      iStrt_tj = 0  ! B
      iStrt_jt = 0  ! G
      iStrt_ui = 0  ! B"
      iStrt_iu = 0  ! G"
      iStrt_aj = 0  ! C
      iStrt_bi = 0  ! C"
      iStrt_tu = 0  ! D
      iStrt_ut = 0  ! D"
      iStrt_au = 0  ! E
      iStrt_bt = 0  ! E"
      iStrt_ab = 0  ! F
      iStrt0_ij = 0  ! A
      iStrt0_ji = 0  ! A"
      iStrt0_tj = 0  ! B
      iStrt0_jt = 0  ! G
      iStrt0_ui = 0  ! B"
      iStrt0_iu = 0  ! G"
      iStrt0_aj = 0  ! C
      iStrt0_bi = 0  ! C"
      iStrt0_tu = 0  ! D
      iStrt0_ut = 0  ! D"
      iStrt0_au = 0  ! E
      iStrt0_bt = 0  ! E"
      iStrt0_ab = 0  ! F
      If ( TCVA  ) then
        Call GetMem('ij','ALLO','REAL',iStrt00_ij,Len_ij)
        iMemTCVX(1,iSym,jSym,1)=iStrt00_ij
        iMemTCVX(1,iSym,jSym,2)=Len_ij
        Call GetMem('ji','ALLO','REAL',iStrt00_ji,Len_ji)
        iMemTCVX(1,jSym,iSym,1)=iStrt00_ji
        iMemTCVX(1,jSym,iSym,2)=Len_ji
      EndIf
      If ( TCVB  ) then
        Call GetMem('tj','ALLO','REAL',iStrt00_tj,Len_tj)
        iMemTCVX(2,iSym,jSym,1)=iStrt00_tj
        iMemTCVX(2,iSym,jSym,2)=Len_tj
        Call GetMem('jt','ALLO','REAL',iStrt00_jt,Len_jt)
        iMemTCVX(7,jSym,iSym,1)=iStrt00_jt
        iMemTCVX(7,jSym,iSym,2)=Len_jt
      EndIf
      If ( TCVBt ) then
        Call GetMem('ui','ALLO','REAL',iStrt00_ui,Len_ui)
        iMemTCVX(2,jSym,iSym,1)=iStrt00_ui
        iMemTCVX(2,jSym,iSym,2)=Len_ui
        Call GetMem('iu','ALLO','REAL',iStrt00_iu,Len_iu)
        iMemTCVX(7,iSym,jSym,1)=iStrt00_iu
        iMemTCVX(7,iSym,jSym,2)=Len_iu
      EndIf
      If ( TCVC  ) then
        Call GetMem('aj','ALLO','REAL',iStrt00_aj,Len_aj)
        iMemTCVX(3,iSym,jSym,1)=iStrt00_aj
        iMemTCVX(3,iSym,jSym,2)=Len_aj
      EndIf
      If ( TCVCt ) then
        Call GetMem('bi','ALLO','REAL',iStrt00_bi,Len_bi)
        iMemTCVX(3,jSym,iSym,1)=iStrt00_bi
        iMemTCVX(3,jSym,iSym,2)=Len_bi
      EndIf
      If ( TCVD  ) then
        Call GetMem('tu','ALLO','REAL',iStrt00_tu,Len_tu)
        iMemTCVX(4,iSym,jSym,1)=iStrt00_tu
        iMemTCVX(4,iSym,jSym,2)=Len_tu
        Call GetMem('ut','ALLO','REAL',iStrt00_ut,Len_ut)
        iMemTCVX(4,jSym,iSym,1)=iStrt00_ut
        iMemTCVX(4,jSym,iSym,2)=Len_ut
      EndIf
      If ( TCVE  ) then
        Call GetMem('au','ALLO','REAL',iStrt00_au,Len_au)
        iMemTCVX(5,iSym,jSym,1)=iStrt00_au
        iMemTCVX(5,iSym,jSym,2)=Len_au
      EndIf
      If ( TCVEt ) then
        Call GetMem('bt','ALLO','REAL',iStrt00_bt,Len_bt)
        iMemTCVX(5,jSym,iSym,1)=iStrt00_bt
        iMemTCVX(5,jSym,iSym,2)=Len_bt
      EndIf
      If ( TCVF  ) then
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
        If ( TCVA  ) iStrt0_ij = iStrt00_ij + (iFBatch-1) * nFVec * Nij
        If ( TCVA  ) iStrt0_ji = iStrt00_ji + (iFBatch-1) * nFVec * Nji
        If ( TCVB  ) iStrt0_tj = iStrt00_tj + (iFBatch-1) * nFVec * Ntj
        If ( TCVB  ) iStrt0_jt = iStrt00_jt + (iFBatch-1) * nFVec * Njt
        If ( TCVBt ) iStrt0_ui = iStrt00_ui + (iFBatch-1) * nFVec * Nui
        If ( TCVBt ) iStrt0_iu = iStrt00_iu + (iFBatch-1) * nFVec * Niu
        If ( TCVC  ) iStrt0_aj = iStrt00_aj + (iFBatch-1) * nFVec * Naj
        If ( TCVCt ) iStrt0_bi = iStrt00_bi + (iFBatch-1) * nFVec * Nbi
        If ( TCVD  ) iStrt0_tu = iStrt00_tu + (iFBatch-1) * nFVec * Ntu
        If ( TCVD  ) iStrt0_ut = iStrt00_ut + (iFBatch-1) * nFVec * Nut
        If ( TCVE  ) iStrt0_au = iStrt00_au + (iFBatch-1) * nFVec * Nau
        If ( TCVEt ) iStrt0_bt = iStrt00_bt + (iFBatch-1) * nFVec * Nbt
        If ( TCVF  ) iStrt0_ab = iStrt00_ab + (iFBatch-1) * nFVec * Nab

*       Allocate memory & Load Full Cholesky Vectors - CHFV
        Len_FAB = NFAB * NumFV
        iStrtVec_FAB = iStrtVec_AB + nFVec * (iFBatch-1)
        Call GetMem('FAB','Allo','Real',iStrt0_FAB,Len_FAB)
        Call RdChoVec(Work(iStrt0_FAB),NFAB,NumFV,iStrtVec_FAB,lUCHFV)

*  ---  Start Loop  iVec  ---
        Do iVec=1,NumFV   ! Loop  iVec
          iStrt_FAB = iStrt0_FAB + (iVec-1) * NFAB
          If ( TCVA  ) iStrt_ij = iStrt0_ij + (iVec-1) * Nij
          If ( TCVA  ) iStrt_ji = iStrt0_ji + (iVec-1) * Nji
          If ( TCVB  ) iStrt_tj = iStrt0_tj + (iVec-1) * Ntj
          If ( TCVB  ) iStrt_jt = iStrt0_jt + (iVec-1) * Njt
          If ( TCVBt ) iStrt_ui = iStrt0_ui + (iVec-1) * Nui
          If ( TCVBt ) iStrt_iu = iStrt0_iu + (iVec-1) * Niu
          If ( TCVC  ) iStrt_aj = iStrt0_aj + (iVec-1) * Naj
          If ( TCVCt ) iStrt_bi = iStrt0_bi + (iVec-1) * Nbi
          If ( TCVD  ) iStrt_tu = iStrt0_tu + (iVec-1) * Ntu
          If ( TCVD  ) iStrt_ut = iStrt0_ut + (iVec-1) * Nut
          If ( TCVE  ) iStrt_au = iStrt0_au + (iVec-1) * Nau
          If ( TCVEt ) iStrt_bt = iStrt0_bt + (iVec-1) * Nbt
          If ( TCVF  ) iStrt_ab = iStrt0_ab + (iVec-1) * Nab

*     --- 1st Half-Transformation  iBeta(AO) -> q(MO) only occupied
          jStrt0MO = 1
          Do j=1,jSym-1
            jStrt0MO = jStrt0MO + nBas(j) * nBas(j)
          EndDo
          jStrt0MO = jStrt0MO + nFro(jSym) * nBas(jSym)

C         From CHFV A(Alpha,Beta) to XAj(Alpha,jMO)
          If ( TCVA .or. TCVB .or. TCVC ) then
            Call GetMem('XAj','ALLO','REAL',iStrt0_XAj,Len_XAj)
            Call ProdsA_2(Work(iStrt_FAB), nBas(iSym),nBas(jSym),
     &                  CMO(jStrt0MO),nIsh(jSym), Work(iStrt0_XAj))
          EndIf
          jStrt0MO = jStrt0MO + nIsh(jSym) * nBas(jSym)

C         From CHFV A(Alpha,Beta) to XAu(Alpha,uMO)
          If ( TCVD .or. TCVE ) then
            Call GetMem('XAu','ALLO','REAL',iStrt0_XAu,Len_XAu)
            Call ProdsA_2(Work(iStrt_FAB), nBas(iSym),nBas(jSym),
     &                  CMO(jStrt0MO),nAsh(jSym), Work(iStrt0_XAu))
          EndIf
          jStrt0MO = jStrt0MO + nAsh(jSym) * nBas(jSym)

C         From CHFV A(Alpha,Beta) to XAb(Alpha,bMO)
          If ( TCVF ) then
            Call GetMem('XAb','ALLO','REAL',iStrt0_XAb,Len_XAb)
            Call ProdsA_2(Work(iStrt_FAB), nBas(iSym),nBas(jSym),
     &                  CMO(jStrt0MO),nSsh(jSym), Work(iStrt0_XAb))
          EndIf

          iStrt0MO = 1
          Do i=1,iSym-1
            iStrt0MO = iStrt0MO + nBas(i) * nBas(i)
          EndDo
          iStrt0MO = iStrt0MO + nFro(iSym) * nBas(iSym)

C         From CHFV A(Alpha,Beta) to XBi(Beta,iMO)
          If ( TCVBt .or. TCVCt ) then
            Call GetMem('XBi','ALLO','REAL',iStrt0_XBi,Len_XBi)
            Call ProdsA_2t(Work(iStrt_FAB), nBas(iSym),nBas(jSym),
     &                  CMO(iStrt0MO),nIsh(iSym), Work(iStrt0_XBi))
          EndIf
          iStrt0MO = iStrt0MO + nIsh(iSym) * nBas(iSym)

C         From CHFV A(Alpha,Beta) to XBt(Beta,tMO)
          If ( TCVEt ) then
            Call GetMem('XBt','ALLO','REAL',iStrt0_XBt,Len_XBt)
            Call ProdsA_2t(Work(iStrt_FAB), nBas(iSym),nBas(jSym),
     &                  CMO(iStrt0MO),nAsh(iSym), Work(iStrt0_XBt))
          EndIf

*     --- 2nd Half-Transformation  iAlpha(AO) -> p(MO)
          iStrt0MO = 1
          Do i=1,iSym-1
            iStrt0MO = iStrt0MO + nBas(i) * nBas(i)
          EndDo
          iStrt0MO = iStrt0MO + nFro(iSym) * nBas(iSym)

C         From XAj(Alpha,jMO) to ij(i,j)
          If ( TCVA ) then
            Call ProdsA_1(Work(iStrt0_XAj), nBas(iSym),nIsh(jSym),
     &                  CMO(iStrt0MO),nIsh(iSym), Work(iStrt_ij))
            Call Trnsps(nIsh(iSym),nIsh(jSym),
     &                                    Work(iStrt_ij),Work(iStrt_ji))
          EndIf
          iStrt0MO = iStrt0MO + nIsh(iSym) * nBas(iSym)

C         From XAj(Alpha,jMO) to tj(t,j)
          If ( TCVB ) then
            Call ProdsA_1(Work(iStrt0_XAj), nBas(iSym),nIsh(jSym),
     &                  CMO(iStrt0MO),nAsh(iSym), Work(iStrt_tj))
            Call Trnsps(nAsh(iSym),nIsh(jSym),
     &                                    Work(iStrt_tj),Work(iStrt_jt))
          EndIf

C         From XAu(Alpha,jMO) to tu(t,u)
          If ( TCVD ) then
            Call ProdsA_1(Work(iStrt0_XAu), nBas(iSym),nAsh(jSym),
     &                  CMO(iStrt0MO),nAsh(iSym), Work(iStrt_tu))
            Call Trnsps(nAsh(iSym),nAsh(jSym),
     &                                    Work(iStrt_tu),Work(iStrt_ut))
          EndIf
          iStrt0MO = iStrt0MO + nAsh(iSym) * nBas(iSym)

C         From XAj(Alpha,jMO) to aj(a,j)
          If ( TCVC ) then
            Call ProdsA_1(Work(iStrt0_XAj), nBas(iSym),nIsh(jSym),
     &                  CMO(iStrt0MO),nSsh(iSym), Work(iStrt_aj))
          EndIf

C         From XAu(Alpha,jMO) to au(a,u)
          If ( TCVE ) then
            Call ProdsA_1(Work(iStrt0_XAu), nBas(iSym),nAsh(jSym),
     &                  CMO(iStrt0MO),nSsh(iSym), Work(iStrt_au))
          EndIf

C         From XAb(Alpha,jMO) to ab(a,b)
          If ( TCVF ) then
            Call ProdsA_1(Work(iStrt0_XAb), nBas(iSym),nSsh(jSym),
     &                  CMO(iStrt0MO),nSsh(iSym), Work(iStrt_ab))
          EndIf

          jStrt0MO = 1
          Do j=1,jSym-1
            jStrt0MO = jStrt0MO + nBas(j) * nBas(j)
          EndDo
          jStrt0MO = jStrt0MO + ( nFro(jSym) + nIsh(jSym) ) * nBas(jSym)

C         From XBi(Beta,jMO) to ui(u,i)
          If ( TCVBt ) then
            Call ProdsA_1(Work(iStrt0_XBi), nBas(jSym),nIsh(iSym),
     &                  CMO(jStrt0MO),nAsh(jSym), Work(iStrt_ui))
            Call Trnsps(nAsh(jSym),nIsh(iSym),
     &                                    Work(iStrt_ui),Work(iStrt_iu))
          EndIf
          jStrt0MO = jStrt0MO + nAsh(jSym) * nBas(jSym)

C         From XBi(Beta,jMO) to bi(b,i)
          If ( TCVCt ) then
            Call ProdsA_1(Work(iStrt0_XBi), nBas(jSym),nIsh(iSym),
     &                  CMO(jStrt0MO),nSsh(jSym), Work(iStrt_bi))
          EndIf

C         From XBt(Beta,jMO) to bt(b,t)
          If ( TCVEt ) then
            Call ProdsA_1(Work(iStrt0_XBt), nBas(jSym),nAsh(iSym),
     &                  CMO(jStrt0MO),nSsh(jSym), Work(iStrt_bt))
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
          If ( TCVBt .or. TCVCt ) then
            Call GetMem('XBi','FREE','REAL',iStrt0_XBi,Len_XBi)
          EndIf
          If ( TCVEt ) then
            Call GetMem('XBt','FREE','REAL',iStrt0_XBt,Len_XBt)
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

      Subroutine ProdsA_1(AB,iA,iB, CMO,nMO, Y)
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Dimension AB(iA,iB), CMO(iA,nMO), Y(nMO,iB)
      Call DGEMM_('T','N',nMO,iB,iA,1.0d0,CMO,iA,AB,iA,0.0d0,Y,nMO)
      Return
      End

      Subroutine ProdsA_2(AB,iA,iB, CMO,nMO, Y)
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Dimension AB(iA,iB), CMO(iB,nMO), Y(iA,nMO)
      Call DGEMM_('N','N',iA,nMO,iB,1.0d0,AB,iA,CMO,iB,0.0d0,Y,iA)
      Return
      End

      Subroutine ProdsA_2t(AB,iA,iB, CMO,nMO, Y)
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Dimension AB(iA,iB), CMO(iA,nMO), Y(iB,nMO)
      Call DGEMM_('T','N',iB,nMO,iA,1.0d0,AB,iA,CMO,iA,0.0d0,Y,iB)
      Return
      End
