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
* Copyright (C) 1994, Per Ake Malmqvist                                *
*               2004,2005, Giovanni Ghigo                              *
************************************************************************
      Subroutine Mem_Est(iSymL,nVec,nFVec)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Torino University, Italy                                   *
*           January-February, July 2005                                *
*----------------------------------------------------------------------*
* Routine for the estimation of the memory usage and the number (nVec) *
* of Cholesky Vectors transformable.                                   *
* The routine also define (calling Def_TCVx) which TCVx create.        *
************************************************************************
*
*   <DOC>
*     <Name>Mem\_Est</Name>
*     <Syntax>Call Mem\_Est(iSymL,nVec,nFVec)
*     </Syntax>
*     <Arguments>
*      \Argument{iSymL}{Symmetry of the Cholesky vector}{Integer}{in}
*      \Argument{nVec}{Number of Cholesky vectors to transform in the
*      batch procedure}{Integer}{out}
*      \Argument{nFVec}{Number of Cholesky vectors to transform in the
*      inner batch procedure}{Integer}{out}
*     </Arguments>
*     <Purpose>
*      The routine makes an estimation of the optimal memory allocation
*      and defines the maximum number of Cholesky vectors to transfor
*      both for the main (nVec) and the inner (nFVec) batch
*      procedures.\\
*      Called by Cho\_TraCtl.
*     </Purpose>
*     <Dependencies>
*       Routine called: Def\_TCVx.
*     </Dependencies>
*     <Author>
*      G. Ghigo
*     </Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
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
      Character*16 Frmt1,Frmt3

      MulD2h(i,j) = iEor(i-1,j-1)+1

      nVec  = 0 ! Nr. of transformable vectors for batch.
      nFVec = 0 ! Nr. of full vectors for inner transformation batch

*   Memory for Reading (will change when reading reduced set)
      LenCHFV = 0     ! Mem for Full Vectors (p.v. = per vector)
*   Memory for Transformation
      LenTCVx   = 0   ! Mem for TCVx (p.v.)
      LenTmpTra = 0   ! Mem for First-half Transformation
*   Memory for Generation
      MaxInt   = 0    ! Mem for integrals
      MaxSlice = 0    ! Mem for Lx (p.v.)

* --- Memory for Transformation
      Nij = 0  ! A
      Nji = 0  ! A"
      Ntj = 0  ! B
      Nui = 0  ! B"
      Naj = 0  ! C
      Nbi = 0  ! C"
      Ntu = 0  ! D
      Nut = 0  ! D"
      Nau = 0  ! E
      Nbt = 0  ! E"
      Nab = 0  ! F
      Njt = 0  ! G
      Niu = 0  ! G"

      Do iSym=1,nSym
       If(nBas(iSym).GT.0) then
        Do jSym=1,iSym
         If(nBas(jSym).GT.0 .and. MulD2h(iSym,jSym).EQ.iSymL) then

          Len_YAj = 0  ! iSym=jSym: A, B, C
          Len_YAu = 0  ! iSym=jSym: D, E
          Len_YAb = 0  ! iSym=jSym: F
          Len_XAj = 0  ! iSym=/=jSym: A/A", B, C
          Len_XAu = 0  ! iSym=/=jSym: D/D", E
          Len_XAb = 0  ! iSym=/=jSym: F
          Len_XBi = 0  ! iSym=/=jSym: B",C"
          Len_XBt = 0  ! iSym=/=jSym: E"
          Len_ABSq = 0 ! Squared CHFV

          Call Def_TCVx(iSym,jSym)

          If (iSym.EQ.jSym) then

* TCV-A :
           If (TCVXist(1,iSym,jSym)) Then
             Len_YAj = nBas(iSym) * nIsh(jSym)
             Nij = Nij + nIsh(iSym) * nIsh(jSym)
           EndIf
* TCV-B :
           If (TCVXist(2,iSym,jSym)) Then
             Len_YAj = nBas(iSym) * nIsh(jSym)
             Ntj = Ntj + ( nAsh(iSym) * nIsh(jSym) )
* TCV-G :
             Njt = Njt + ( nIsh(jSym) * nAsh(iSym) )
           EndIf
* TCV-C :
           If (TCVXist(3,iSym,jSym)) Then
             Len_YAj = nBas(iSym) * nIsh(jSym)
             Naj = Naj + nSsh(iSym) * nIsh(jSym)
           EndIf
* TCV-D :
           If (TCVXist(4,iSym,jSym)) Then
             Len_YAu = nBas(iSym) * nAsh(jSym)
             Ntu = Ntu + nAsh(iSym) * nAsh(jSym)
           EndIf
* TCV-E :
           If (TCVXist(5,iSym,jSym)) Then
             Len_YAu = nBas(iSym) * nAsh(jSym)
             Nau = Nau + nSsh(iSym) * nAsh(jSym)
           EndIf
* TCV-F :
           If (TCVXist(6,iSym,jSym)) Then
             Len_YAb = nBas(iSym) * nSsh(jSym)
             Nab = Nab + nSsh(iSym) * nssh(jSym)
           EndIf

           LenCHFV = LenCHFV + ( nBas(iSym) * ( nBas(jSym) + 1 ) / 2 )
           Len_ABSq = nBas(iSym) * nBas(jSym)
           LenTmpTra = Max( LenTmpTra ,
     &                      ( Len_ABSq + Len_YAj + Len_YAu + Len_YAb ) )

          else

* TCV-A :
           If (TCVXist(1,iSym,jSym)) Then
             Len_XAj = nBas(iSym) * nIsh(jSym)
             Nij = Nij + ( nIsh(iSym) * nIsh(jSym) )
             Nji = Nji + ( nIsh(jSym) * nIsh(iSym) )
           EndIf
* TCV-B :
           If (TCVXist(2,iSym,jSym)) Then
             Len_XAj = nBas(iSym) * nIsh(jSym)
             Ntj = Ntj + ( nAsh(iSym) * nIsh(jSym) )
* TCV-G :
             Njt = Njt + ( nIsh(jSym) * nAsh(iSym) )
           EndIf
           If (TCVXist(2,jSym,iSym)) Then
             Len_XBi = nBas(jSym) * nIsh(iSym)
             Nui = Nui + ( nAsh(jSym) * nIsh(iSym) )
* TCV-G :
             Niu = Niu + ( nIsh(iSym) * nAsh(jSym) )
           EndIf
* TCV-C :
           If (TCVXist(3,iSym,jSym)) Then
             Len_XAj = nBas(iSym) * nIsh(jSym)
             Naj = Naj + nSsh(iSym) * nIsh(jSym)
           EndIf
           If (TCVXist(3,jSym,iSym)) Then
             Len_XBi = nBas(jSym) * nIsh(iSym)
             Nbi = Nbi + nSsh(jSym) * nIsh(iSym)
           EndIf
* TCV-D :
           If (TCVXist(4,iSym,jSym)) Then
             Len_XAu = nBas(iSym) * nAsh(jSym)
             Ntu = Ntu + ( nAsh(iSym) * nAsh(jSym) )
             Nut = Nut + ( nAsh(jSym) * nAsh(iSym) )
           EndIf
* TCV-E :
           If (TCVXist(5,iSym,jSym)) Then
             Len_XAu = nBas(iSym) * nAsh(jSym)
             Nau = Nau + nSsh(iSym) * nAsh(jSym)
           EndIf
           If (TCVXist(5,jSym,iSym)) Then
             Len_XBt = nBas(jSym) * nAsh(iSym)
             Nbt = Nbt + nSsh(jSym) * nAsh(iSym)
           EndIf
* TCV-F :
           If (TCVXist(6,iSym,jSym)) Then
             Len_XAb = nBas(iSym) * nSsh(jSym)
             Nab = Nab + nSsh(iSym) * nSsh(jSym)
           EndIf

           LenCHFV = LenCHFV + ( nBas(iSym) * nBas(jSym) )
           LenTmpTra = Max ( LenTmpTra ,
     &             ( Len_XAj + Len_XAu + Len_XAb + Len_XBi + Len_XBt ) )

          EndIf

         EndIf
        EndDo ! jSym
       EndIf
      EndDo ! iSym

      LenTCVx = Nij+Nji + Ntj+Nui + Naj+Nbi + Ntu+Nut + Nau+Nbt +
     &                                                     Nab + Niu+Njt

* --- Memory for Generation
      Do iSymI = 1, nSym
       Do iSymJ = 1, iSymI
        Do iSymA = 1, nSym
         Do iSymB = 1, iSymA
          iSymAI = MulD2h(iSymA,iSymI)
          iSymBJ = MulD2h(iSymB,iSymJ)
          Call LenInt(iSymI,iSymJ,iSymA,iSymB,nN_IJ,nN_AB,nN_Ex1,nN_Ex2)
          If (iSymAI.EQ.iSymL.and.iSymBJ.EQ.iSymL. and.
     &                                        nN_IJ*nN_AB.GT.0) then
            MaxInt    = Max(MaxInt, Max(nN_AB, Max(nN_Ex1, 2*nN_Ex2) ) )
            MaxSlice  = Max(MaxSlice, nOrb(iSymA)+nOrb(iSymB) )
          EndIf
         EndDo
        EndDo
       EndDo
      EndDo

      Call GetMem('MaxMem','MAX','REAL',KDUM,MEMX)
      MemFree0 = Max( MEMX - MEMX/10 , 0 )

      MemPerVec2 = LenTCVx + MaxSlice
      MemTmp2 = 2 * MaxInt

      MemMin = LenCHFV + LenTmpTra + LenTCVx + MaxSlice + 2 * MaxInt
      MemFree = MemFree0 - MemMin

      nVec = Min( ( MemFree / MemPerVec2 + 1 ) , NumCho(iSymL) )
      MemFree = MemFree - (nVec-1) * MemPerVec2
      nFVec = Min( ( MemFree / LenCHFV + 1 ) , NumCho(iSymL) )
      MemAlloc = MemMin + (nVec-1) * MemPerVec2 + (nFVec-1) * LenCHFV

      If (.NOT.IfTest.and.nVec.EQ.NumCho(iSymL).and.nFVec.GT.0) Return

      MaxNum = Max(MemFree0,MemAlloc)
      MaxSize = Max(9,Int(Log10(Dble(MaxNum)))+1)
      Write(Frmt1,'(I16)') MaxSize
      Frmt1='(A,1X,I'//Trim(AdjustL(Frmt1))//')'
      Write(Frmt3,'(I16)') MaxSize-5
      Frmt3='(A,1X,I'//Trim(AdjustL(Frmt3))//',A)'

      Write(6,*)
      Write(6,'(20A3)')('---',I=1,20)
      Write(6,05)' MEM for TRANSF/GENER of VECTORS in Sym:',iSymL
      Write(6,*)
      Write(6,Frmt1)' Mem (p.v.) for CHFV            :',LenCHFV
      Write(6,Frmt1)' Tmp Mem for transformation     :',LenTmpTra
      Write(6,Frmt1)' Mem (p.v.) for TCVx            :',LenTCVx
      Write(6,Frmt1)' Max Tmp Mem (p.v.) for generat.:',MaxSlice
      Write(6,Frmt1)' Max Mem for Integrals (twice)  :',2*MaxInt
      Write(6,Frmt1)
      Write(6,Frmt1)' TOTAL AVAILABLE MEMORY         :',MemFree0
      Write(6,Frmt1)' Minimal Memory required        :',MemMin
      Write(6,Frmt1)' Memory Required by generation  :',
     & (nVec-1)*MemPerVec2
      Write(6,Frmt1)' Memory Available for transform.:',MemFree
      Write(6,Frmt1)' Memory Required by transform.  :',
     & (nFVec-1)*LenCHFV
      Write(6,Frmt1)' TOTAL ALLOCATED MEMORY         :',MemAlloc
      Write(6,Frmt1)' Unemployed memory              :',
     & MemFree0-MemAlloc
      Write(6,*)
      Write(6,20)' Max nr. of Transformed vectors  :',nVec
      Write(6,20)' Max nr. of vectors in sub-batch :',nFVec
      Write(6,20)' Total Number of Cholesky vectors:',NumCho(iSymL)
      Write(6,*)
      Write(6,*) 'ESTIMATED MEMORY REQUIREMENTS'
      Write(6,Frmt3)'  Minimum:',2+MemMin / 119000,' MB'
      Write(6,Frmt3)'  Normal :',2+
     & (MemMin+MemPerVec2 * (NumCho(iSymL)-1) ) / 119000,' MB'
      Write(6,Frmt3)'  Maximum:',2+
     & (MemMin+(MemPerVec2+LenCHFV) * (NumCho(iSymL)-1)) / 119000,' MB'
      Write(6,'(20A3)')('---',I=1,20)
      Write(6,*)
      Call XFlush(6)
      If ( nVec.LT.NumCho(iSymL) ) then
        Write(6,*)' Batch procedure used. Increase memory if possible!'
        Write(6,*)
        Write(6,'(20A3)')('---',I=1,20)
        Call XFlush(6)
      EndIf
05    Format(A,1X,I1)
20    Format(A,1X,I6)

      Return
      End
      Subroutine Def_TCVx(iSym,jSym)
************************************************************************
* Author  :  Giovanni Ghigo                                            *
*            Lund University, Sweden                                   *
* Written :  October 2004                                              *
* Modified:  January 2005                                              *
*----------------------------------------------------------------------*
* Define which Transformed Cholesky Full Vectors (TCVx) to generate.   *
* TCVXist(iType,iSym,jSym) is .True. if the TCVx must be generated.    *
* iType(x):  1=A, 2=B, 3=C, 4=D, 5=E, 6=F, 7=G                         *
************************************************************************
*
*   <DOC>
*     <Name>Def\_TCVx</Name>
*     <Syntax>Call Def\_TCVx(iSym,jSym)
*     </Syntax>
*     <Arguments>
*      \Argument{iSym}{Symmetry(i) of the Cholesky full vector}{Integer}
*      {in}
*      \Argument{jSym}{Symmetry(j) of the Cholesky full vector}{Integer}
*      {in}
*     </Arguments>
*     <Purpose>
*      The routine defines which Transformed Cholesky (TCVx) to generate
*      setting .True. the logical
*      matrix TCVXist(iType,iSym,jSym) iType: 1=TCVA, 2=TCVB, 3=TCVC,
*      4=TCVD, 5=TCVE, 6=TCVF, 7=TCVG.\\
*      Called by Mem\_Est.
*     </Purpose>
*     <Dependencies>
*     </Dependencies>
*     <Author>
*      G. Ghigo
*     </Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*     </Description>
*    </DOC>
*
******************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "cho_tra.fh"

      If(nIsh(jSym).GT.0) then
        If(nIsh(iSym).GT.0 .and. DoTCVA) then
          TCVXist(1,iSym,jSym)=.True.
          TCVXist(1,jSym,iSym)=.True.  ! Aji = T(Aij) !
        EndIf
        If(nAsh(iSym).GT.0 .and. DoTCVA) then
          TCVXist(2,iSym,jSym)=.True.
          TCVXist(7,jSym,iSym)=.True.
        EndIf
        If(nSsh(iSym).GT.0) then
          TCVXist(3,iSym,jSym)=.True.
        EndIf
      EndIf

      If(nAsh(jSym).GT.0 .and. DoTCVA) then
        If(nIsh(iSym).GT.0 .and. iSym.NE.jSym) then
          TCVXist(2,jSym,iSym)=.True.
          TCVXist(7,iSym,jSym)=.True.
        EndIf
        If(nAsh(iSym).GT.0) then
          TCVXist(4,iSym,jSym)=.True.
          TCVXist(4,jSym,iSym)=.True.  ! Dji = T(Dij) !
        EndIf
        If(nSsh(iSym).GT.0) then
          TCVXist(5,iSym,jSym)=.True.
        EndIf
      EndIf

      If(nSsh(jSym).GT.0 .and. iSym.NE.jSym) then
        If(nIsh(iSym).GT.0) then
          TCVXist(3,jSym,iSym)=.True.
        EndIf
        If(nAsh(iSym).GT.0 .and. DoTCVA) then
          TCVXist(5,jSym,iSym)=.True.
        EndIf
      EndIf

      If(nSsh(jSym).GT.0 .and. nSsh(iSym).GT.0 .and. DoTCVA) then
        TCVXist(6,iSym,jSym)=.True.
      EndIf

      Return
      End

      Subroutine LenInt(iSymI,iSymJ,iSymA,iSymB,nProdIJ,
     &                                          nProdAB,nProdE1,nProdE2)
************************************************************************
* Author  :  Giovanni Ghigo                                            *
*            Lund University, Sweden                                   *
*----------------------------------------------------------------------*
* Return the Length of Coulomb (nProdAB), Exchanges (nProdE1 and       *
* nProdE2) matrices for each i,j and the length (nProdIJ) of the i,j   *
* matrix for each Symmetry Block (iSymI,iSymJ,iSymA,iSymB)             *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "cho_tra.fh"
      nProdIJ=0
      nProdAB=0
      nProdE1=0
      nProdE2=0
      nOccI = nOsh(iSymI)
      nOccJ = nOsh(iSymJ)
      nOrbA = nOrb(iSymA)
      nOrbB = nOrb(iSymB)
      nOrbA2= nOrb(iSymB)
      nOrbB2= nOrb(iSymA)
      nExtA = nSsh(iSymA)
      nExtB = nSsh(iSymB)
      nExtA2= nSsh(iSymB)
      nExtB2= nSsh(iSymA)
      If(iSymI.EQ.iSymJ) then
        nProdIJ=nOccI*(nOccJ+1)/2
      else
        nProdIJ=nOccI*nOccJ
      EndIf
      If(iSymA.EQ.iSymB) then
        nProdAB=nOrbA*(nOrbB+1)/2
      else
        If(iSymA.GT.iSymB) then
          nProdAB=nOrbA*nOrbB
        else
          nProdAB=0
        EndIf
      EndIf
      If(iSymA.GE.iSymB) then
        If(DoTCVA) then
          nProdE1=nOrbA*nOrbB
        else
          nProdE1=nExtA*nExtB
        EndIf
        nProdE2=0
      else
        nProdE1=0
        If(DoTCVA) then
          nProdE2=nOrbA*nOrbB
        else
          nProdE2=nExtA*nExtB
        EndIf
      EndIf
      Return
      End

      Subroutine Def_SubBlockE(iSymA,iSymB)
************************************************************************
* Author  :  Giovanni Ghigo                                            *
*            Lund University, Sweden                                   *
*----------------------------------------------------------------------*
* Define the SubBlocks to calculate in the Exchange matrix.            *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "cho_tra.fh"
      Do i=1,3
        Do j=1,3
          SubBlocks(i,j)=.False.
        EndDo
      EndDo
      If (DoTCVA .and. nIsh(iSymA).GT.0) then
        If (nIsh(iSymB).GT.0) SubBlocks(1,1)=.True.
        If (nAsh(iSymB).GT.0) SubBlocks(1,2)=.True.
        If (nSsh(iSymB).GT.0) SubBlocks(1,3)=.True.
      EndIf
      If (DoTCVA .and. nAsh(iSymA).GT.0) then
        If (nIsh(iSymB).GT.0) SubBlocks(2,1)=.True.
        If (nAsh(iSymB).GT.0) SubBlocks(2,2)=.True.
        If (nSsh(iSymB).GT.0) SubBlocks(2,3)=.True.
      EndIf
      If (DoTCVA .and. nSsh(iSymA).GT.0) then
        If (nIsh(iSymB).GT.0) SubBlocks(3,1)=.True.
        If (nAsh(iSymB).GT.0) SubBlocks(3,2)=.True.
      EndIf
      If ((nSsh(iSymA)*nSsh(iSymB)).GT.0) SubBlocks(3,3)=.True.
      Return
      End

      Subroutine Local_Triang(nRow,A)
C This routine is a modification of the Per-Ake's Triang routine
C found in src/caspt2/triang.f
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Dimension A(nRow**2)
c Convert a square matrix to triangular in-place.
      iFrom=1+nRow
      iTo=2
      Do i=2,nRow
        Call dCopy_(i,A(iFrom),1,A(iTo),1)
        iFrom=iFrom+nRow
        iTo=iTo+i
      EndDo
      Return
      End

      Subroutine PrintSquareMat(nRow,A)
* Prints a square matrix A(nRow,nRow)
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Dimension A(nRow**2)
      If (nRow.GT.8) Return
      iCount=0
      Do i=1,nRow
        If(nRow.EQ.1) Write(6,'(1F10.6)') (A(iCount+k),k=1,nRow)
        If(nRow.EQ.2) Write(6,'(2F10.6)') (A(iCount+k),k=1,nRow)
        If(nRow.EQ.3) Write(6,'(3F10.6)') (A(iCount+k),k=1,nRow)
        If(nRow.EQ.4) Write(6,'(4F10.6)') (A(iCount+k),k=1,nRow)
        If(nRow.EQ.5) Write(6,'(5F10.6)') (A(iCount+k),k=1,nRow)
        If(nRow.EQ.6) Write(6,'(6F10.6)') (A(iCount+k),k=1,nRow)
        If(nRow.EQ.7) Write(6,'(7F10.6)') (A(iCount+k),k=1,nRow)
        If(nRow.EQ.8) Write(6,'(8F10.6)') (A(iCount+k),k=1,nRow)
        iCount=iCount+nRow
      EndDo
      Return
      End

      Subroutine PrintDiagMat(nRow,A)
* Prints a diagonal matrix A(nRow,nRow)
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Dimension A(nRow*(nRow+1))
      If (nRow.GT.8) Return
      iCount=0
      Do i=1,nRow
        Write(6,'(8F10.6)') (A(iCount+k),k=1,i)
        iCount=iCount+i
      EndDo
      Return
      End
