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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine ChoMP2g_Unpack_Pairs(iSym,iVec,iSymP,nTq,Buf,lBuf,kBuf,Out,ldOut,nCol)
!***********************************************************************
!                                                                      *
! Purpose: unpack one MO Cholesky vector, stored as the nine orbital-  *
!          pair types of the MP2 gradient code, into a full block      *
!          L_pq^J with p in irrep iSymP (rows) and q in irrep          *
!          iSymQ = iSymP x iSym (columns).                             *
!                                                                      *
! Storage of the pair-type vectors (see ChoMP2g_Setup and              *
! ChoMP2g_TraVec): orbital types t = 1 (frozen), 2 (occupied),         *
! 3 (virtual). The vector of pair type 3*(t_q-1)+t_p in irrep iSym has *
! nMoMo(iSym,type) elements per Cholesky vector, ordered in blocks by  *
! the irrep of q; the block (iSymP,iSymQ) starts at                    *
! iMoMo(iSymP,iSymQ,type) and is nMO(iSymP,t_p) x nMO(iSymQ,t_q) with  *
! p running fastest.                                                   *
!                                                                      *
! In the output block the rows are the orbitals of irrep iSymP in the  *
! order frozen, occupied, virtual (nOrb(iSymP) rows), and the columns  *
! the orbitals of irrep iSymQ of the first nTq types in the same order *
! (nTq = 3: all orbitals; nTq = 2: frozen and occupied).               *
!                                                                      *
! Arguments:                                                           *
!   iSym  (in)  irrep of the Cholesky vectors                          *
!   iVec  (in)  index of the vector within the batch held in Buf       *
!   iSymP (in)  irrep of the row orbitals p                            *
!   nTq   (in)  number of column orbital types (2 or 3)                *
!   Buf   (in)  pair-type vectors of the batch; type iType starts at   *
!               Buf(1+kBuf(iType)) and holds nMoMo(iSym,iType) values  *
!               per vector                                             *
!   lBuf  (in)  length of Buf                                          *
!   kBuf  (in)  offsets of the nine pair types in Buf                  *
!   Out   (out) the block, ldOut x nCol (ldOut >= nOrb(iSymP))         *
!   ldOut (in)  leading dimension of Out                               *
!   nCol  (in)  number of columns of Out (>= orbitals of iSymQ of the  *
!               first nTq types)                                       *
!                                                                      *
!***********************************************************************

use Symmetry_Info, only: Mul
use ChoMP2, only: iMoMo, nFro, nMoMo, nOcc, nVir
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iSym, iVec, iSymP, nTq, lBuf, kBuf(9), ldOut, nCol
real(kind=wp), intent(in) :: Buf(lBuf)
real(kind=wp), intent(out) :: Out(ldOut,nCol)
integer(kind=iwp) :: i, iOff, iOffP(3), iOffQ(3), iSymQ, iTp, iTq, iType, n, nMOP(3), nMOQ(3)

iSymQ = Mul(iSym,iSymP)

! Number of orbitals of each type and their offset within the irrep
nMOP(:) = [nFro(iSymP),nOcc(iSymP),nVir(iSymP)]
nMOQ(:) = [nFro(iSymQ),nOcc(iSymQ),nVir(iSymQ)]
iOffP(:) = [0,nMOP(1),nMOP(1)+nMOP(2)]
iOffQ(:) = [0,nMOQ(1),nMOQ(1)+nMOQ(2)]

Out(:,:) = Zero
do iTq=1,nTq
  do i=1,nMOQ(iTq)
    do iTp=1,3
      n = nMOP(iTp)
      if (n == 0) cycle
      iType = 3*(iTq-1)+iTp
      iOff = kBuf(iType)+nMoMo(iSym,iType)*(iVec-1)+iMoMo(iSymP,iSymQ,iType)+(i-1)*n
      Out(iOffP(iTp)+1:iOffP(iTp)+n,iOffQ(iTq)+i) = Buf(iOff+1:iOff+n)
    end do
  end do
end do

end subroutine ChoMP2g_Unpack_Pairs
