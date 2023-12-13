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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_GetZ(irc,NVT,l_NVT,nBlock,l_nBlock,nV,l_nV1,l_nV2,iV1,l_iV11,l_iV12,ip_Z,l_Z1,l_Z2,Z,l_Z)
!
! Thomas Bondo Pedersen, April 2010.
!
! Purpose: get the Z vectors in core.
!
! Input:
!    NVT(i): total no. of vectors in symmetry i
!    nBlock(i): no. of vector blocks in symmetry i
!    nV(j,i): number of vectors in block j of symmetry i
!    iV1(j,i): first vector in block j of symmetry i
!    ip_Z(j,i): pointer to triangular block j of symmetry i
!               (here, j is a compound index j=iTri(k,l) for
!                block k,l)
!
! On exit, the Z vector blocks are stored in memory according to ip_Z.

use Index_Functions, only: iTri, nTri_Elem
use Cholesky, only: iiBstR, InfVec, LuPri, nnBstR, nSym, nSys_call, NumCho, TDECOM
#ifdef _DEBUGPRINT_
use Cholesky, only: iPrint
#endif
use Cholesky_procedures, only: Cho_X_GetIP_InfVec
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: l_NVT, NVT(l_NVT), l_nBlock, nBlock(l_nBlock), l_nV1, l_nV2, nV(l_NV1,l_NV2), l_iV11, l_iV12, &
                                 iV1(l_IV11,l_iV12), l_Z1, l_Z2, ip_Z(l_Z1,l_Z2), l_Z
real(kind=wp), intent(out) :: Z(l_Z)
integer(kind=iwp) :: idRS2RS, iJ, iLoc, iRed, iRedC, iSym, j, J_inBlock, jBlock, k, K_inBlock, kBlock, KK, KK1, KKK, kOffV, kOffZ, &
                     l_Wrk, mUsed, nVRead
real(kind=wp) :: C0, C1, W0, W1
integer(kind=iwp), pointer :: InfVct(:,:,:)
integer(kind=iwp), allocatable :: iRS2RS(:)
real(kind=wp), allocatable :: Wrk(:)
character(len=*), parameter :: SecNam = 'Cho_GetZ'
integer(kind=iwp), external :: Cho_iRange
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i, n, nBlock_Max, nnB
integer(kind=iwp), parameter :: myDebugInfo = 100
real(kind=wp), parameter :: Tol = 1.0e-14_wp
#endif

#ifndef _DEBUGPRINT_
#include "macros.fh"
unused_var(NVT)
#endif

! Set return code.
! ----------------

irc = 0

#ifdef _DEBUGPRINT_
! Check input variables
if ((l_NVT < nSym) .or. (l_nBlock < nSym) .or. (l_nV2 < nSym) .or. (l_iV12 < nSym) .or. (l_Z2 < nSym)) then
  irc = -1
  return
end if
nBlock_Max = nBlock(1)
do iSym=2,nSym
  nBlock_Max = max(nBlock_Max,nBlock(iSym))
end do
nnB = nTri_Elem(nBlock_Max)
if ((l_nV1 < nBlock_Max) .or. (l_iV11 < nBlock_Max) .or. (l_Z1 < nnB)) then
  irc = -2
  return
end if
if (iPrint >= myDebugInfo) then
  call Cho_Head('Entering '//SecNam//':','*',80,LuPri)
  write(LuPri,'(A,8I6)') 'NVT   :',(NVT(iSym),iSym=1,nSym)
  write(LuPri,'(A,8I6)') 'nBlock:',(nBlock(iSym),iSym=1,nSym)
end if
n = 0
do iSym=1,nSym
  if (iPrint >= myDebugInfo) then
    do i=1,nBlock(iSym)
      write(LuPri,'(2X,A,I6,A,I2,A,I6,A,I6)') 'Block',i,' of sym.',iSym,': nV=',nV(i,iSym),' iV1=',iV1(i,iSym)
    end do
  end if
  do j=1,nBlock(iSym)
    do i=j,nBlock(iSym)
      if (i == j) then
        k = nTri_Elem(nV(i,iSym))
      else
        k = nV(i,iSym)*nV(j,iSym)
      end if
      n = n+k
      if (iPrint >= myDebugInfo) then
        write(Lupri,'(5X,A,I6,A,I6,A,I2,A,I9)') 'iBlock',i,' jBlock',j,' of Sym.',iSym,': ip_Z=',ip_Z(iTri(i,j),iSym)
        write(LuPri,'(5X,A,I6,A,I9)') '--> Block dimension:',k,'  Block ends at:',ip_Z(iTri(i,j),iSym)+k-1
      end if
    end do
  end do
end do
k = nTri_Elem(NVT(1))
do iSym=2,nSym
  k = k+nTri_Elem(NVT(iSym))
end do
if (iPrint >= myDebugInfo) then
  write(Lupri,'(A,I8)') 'Total dimension of Z (from blocks):',n
  write(Lupri,'(A,I8)') 'Total dimension of Z (from NVT)   :',k
  write(Lupri,'(A,I8)') '                        Difference:',n-k
end if
if (n /= k) then
  irc = -3
  if (iPrint >= myDebugInfo) then
    write(LuPri,'(A)') 'Ooops! They disagree....'
    write(LuPri,'(A,I4)') 'Returning with code:',irc
  end if
  return
end if
#endif

! Zero result array.
! ------------------

do iSym=1,nSym
  do kBlock=1,nBlock(iSym)
    Z(ip_Z(iTri(kBlock,kBlock),iSym):ip_Z(iTri(kBlock,kBlock),iSym)+nTri_Elem(nV(kBlock,iSym))-1) = Zero
    do jBlock=kBlock+1,nBlock(iSym)
      Z(ip_Z(iTri(jBlock,kBlock),iSym):ip_Z(iTri(jBlock,kBlock),iSym)+nV(jBlock,iSym)*nV(kBlock,iSym)-1) = Zero
    end do
  end do
end do

! Scratch location in index arrays.
! ---------------------------------

iLoc = 3 ! do NOT change (used implicitly by reading routine)

! Get pointer to InfVec array for all vectors.
! Needed for parallel runs.
! --------------------------------------------

call Cho_X_GetIP_InfVec(InfVcT)

! Copy rs1 to location 2.
! -----------------------

! Note: location 2 must contain rs1 throughout this routine! It is
! an assumption which is never checked, so do NOT change this call
! or modify contents of location 2.
call Cho_X_RSCopy(irc,1,2)
if (irc /= 0) then
  write(LuPri,'(A,A,I5)') SecNam,': Cho_X_RSCopy returned code',irc
  irc = 1
  return ! exit
end if

! Get Z vectors.
! --------------

iRedC = -1
do iSym=1,nSym

  call mma_allocate(iRS2RS,nnBstR(iSym,1),Label='iRS2RS')
  call mma_maxDBLE(l_Wrk)
  call mma_allocate(Wrk,l_Wrk,Label='Wrk')
  iRS2RS(:) = 0
  idRS2RS = -2
  KK1 = 1
  do while (KK1 <= NumCho(iSym))
    nVRead = 0
    mUsed = 0
    call CWTime(C0,W0)
    call Cho_X_VecRd(Wrk,size(Wrk),KK1,NumCho(iSym),iSym,nVRead,iRedC,mUsed)
    if (nVRead < 1) then
      irc = 2
      return ! exit
    end if
    call CWTime(C1,W1)
    tDecom(1,2) = tDecom(1,2)+(C1-C0)
    tDecom(2,2) = tDecom(2,2)+(W1-W0)
    nSys_Call = nSys_Call+1
    kOffV = 0
    do KKK=0,nVRead-1
      KK = KK1+KKK
      iRed = InfVec(KK,2,iSym)
      if (iRedC /= iRed) then
        call Cho_X_SetRed(irc,iLoc,iRed)
        if (irc /= 0) then
          irc = 3
          return ! exit
        end if
        iRedC = iRed
      end if
      if (idRS2RS /= iRedC) then
        call Cho_RS2RS(iRS2RS,size(iRS2RS),2,iLoc,iRedC,iSym)
        idRS2RS = iRedC
      end if
      K = InfVec(KK,5,iSym)
      kBlock = Cho_iRange(K+1,iV1(1,iSym),nBlock(iSym),.true.)
#     ifdef _DEBUGPRINT_
      if ((kBlock < 1) .or. (kBlock > nBlock(iSym))) then
        call Cho_Quit('[BLOCK] Error detected in '//SecNam,104)
      end if
#     endif
      K_inBlock = K-iV1(kBlock,iSym)+1
      kOffZ = ip_Z(iTri(kBlock,kBlock),iSym)-1
      do J_inBlock=K_inBlock,nV(kBlock,iSym)
        J = iV1(kBlock,iSym)+J_inBlock-1
        iJ = iRS2RS(InfVcT(J,1,iSym)-iiBstR(iSym,1))
        Z(kOffZ+iTri(J_inBlock,K_inBlock)) = Wrk(kOffV+iJ)
      end do
      do jBlock=kBlock+1,nBlock(iSym)
        kOffZ = ip_Z(iTri(jBlock,kBlock),iSym)-1+nV(jBlock,iSym)*(K_inBlock-1)
        do J_inBlock=1,nV(jBlock,iSym)
          J = iV1(jBlock,iSym)+J_inBlock-1
          iJ = iRS2RS(InfVcT(J,1,iSym)-iiBstR(iSym,1))
          Z(kOffZ+J_inBlock) = Wrk(kOffV+iJ)
        end do
      end do
      kOffV = kOffV+nnBstR(iSym,iLoc)
    end do
    KK1 = KK1+nVRead
  end do
  call mma_deallocate(Wrk)
  call mma_deallocate(iRS2RS)

end do

! Synchronize result array across nodes.
! --------------------------------------

do iSym=1,nSym
  do kBlock=1,nBlock(iSym)
    call Cho_GAdGOp(Z(ip_Z(iTri(kBlock,kBlock),iSym)),nTri_Elem(nV(kBlock,iSym)),'+')
    do jBlock=kBlock+1,nBlock(iSym)
      call Cho_GAdGOp(Z(ip_Z(iTri(jBlock,kBlock),iSym)),nV(jBlock,iSym)*nV(kBlock,iSym),'+')
    end do
  end do
end do

#ifdef _DEBUGPRINT_
! Check that diagonal elements of Z are strictly positive.
! --------------------------------------------------------

if (iPrint >= myDebugInfo) then
  call Cho_Head(SecNam//': Diagonal of Z Vector Matrix','=',80,LuPri)
  write(Lupri,'(A)') ' '
end if
n = 0
do iSym=1,nSym
  do jBlock=1,nBlock(iSym)
    kOffZ = ip_Z(iTri(jBlock,jBlock),iSym)-1
    do J_inBlock=1,nV(jBlock,iSym)
      if (iPrint >= myDebugInfo) then
        write(LuPri,'(A,I2,A,I6,A,ES15.6,A,ES15.6)') 'Sym=',iSym,'  J=',iV1(jBlock,iSym)+J_inBlock-1,'  Z(J,J)=', &
                                                     Z(kOffZ+iTri(J_inBlock,J_inBlock)),'  Squared=', &
                                                     Z(kOffZ+iTri(J_inBlock,J_inBlock))**2
      end if
      if ((abs(Z(kOffZ+iTri(J_inBlock,J_inBlock))) < Tol) .or. (Z(kOffZ+iTri(J_inBlock,J_inBlock)) < -Tol)) then
        n = n+1
        if (iPrint >= myDebugInfo) write(LuPri,'(A)') '  --> Small or negative Z diagonal!'
      end if
    end do
  end do
end do
if (n /= 0) then
  irc = 20
  return ! return
end if
! Check diagonal elements
call Cho_CheckDiagFromZ(irc,NVT,l_NVT,nBlock,l_nBlock,nV,l_nV1,l_nV2,iV1,l_iV11,l_iV12,ip_Z,l_Z1,l_Z2,Z,l_Z,iPrint >= myDebugInfo)
if (irc /= 0) return ! return
#endif

end subroutine Cho_GetZ
