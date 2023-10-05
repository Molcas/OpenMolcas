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
! Copyright (C) 2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoLSOSMP2_Energy_Srt(N,w,t,EOcc,EVir,Delete,EMP2,irc)
!
! Thomas Bondo Pedersen, December 2012.
!
! Compute Laplace-SOS-MP2 energy correction from sorted Cholesky
! vectors (i.e., occupied orbitals processed in batches).

use Symmetry_Info, only: Mul
use Index_Functions, only: nTri_Elem
use Cholesky, only: nSym, NumCho
use ChoMP2, only: DecoMP2, iFirstS, iOcc, iVir, Laplace_BlockSize, Laplace_nGridPoints, LiT1am, LnOcc, LnT1am, lUnit, nBatch, &
                  nMP2Vec, nT1am, nVir
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: w(N), t(N), EOcc(*), EVir(*)
logical(kind=iwp), intent(in) :: Delete
real(kind=wp), intent(out) :: EMP2
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp) :: a, blast, bsize, i, iAddr, iBatch, iBlock, ii, iOpt, ip0, ip1, ipi, ipj, ipX, iSym, iSyma, iSymi, iVec, &
                     jBlock, l_Tot, l_X, lenX, Nai, nBlock, nEnrVec(8), nVeci, nVecj, q
real(kind=wp) :: Eq, lX, tq, wq, xb, xbp, xM, xn
#if !defined (_I8_) || defined (_DEBUGPRINT_)
real(kind=wp) :: Byte
character(len=2) :: Unt
#endif
real(kind=wp), allocatable :: V(:), X(:)
character(len=*), parameter :: SecNam = 'ChoLSOSMP2_Energy_Srt'
real(kind=wp), external :: dDot_

! init return code
irc = 0

! init energy
EMP2 = Zero

! check input (incl. common block variables)
if (nBatch < 2) then
  irc = -1
  return
end if
if (N /= Laplace_nGridPoints) then
  irc = -2
  return
end if
if (Laplace_BlockSize < 1) then
  irc = -3
  return
end if

! set number of vectors
if (DecoMP2) then
  nEnrVec(1:nSym) = nMP2Vec(1:nSym)
else
  nEnrVec(1:nSym) = NumCho(1:nSym)
end if

! allocate X
lX = Zero
do iSym=1,nSym
  if ((nT1am(iSym) > 0) .and. (nEnrVec(iSym) > 0)) then
    bsize = min(Laplace_BlockSize,nEnrVec(iSym))
    nBlock = (nEnrVec(iSym)-1)/bsize+1
    blast = nEnrVec(iSym)-bsize*(nBlock-1)
    xM = real(nEnrVec(iSym),kind=wp)
    xn = real(nBlock,kind=wp)
    xb = real(bsize,kind=wp)
    xbp = real(blast,kind=wp)
    lX = max(lX,Half*(xM*(xM+One)+(xn-One)*xb*(xb-One)+xbp*(xbp-One)))
  end if
end do
l_X = int(lX)
#if !defined (_I8_) || defined (_DEBUGPRINT_)
if (l_X < 0) then
  write(Lupri,'(A,A)') SecNam,': dimension of X matrix is negative!'
  write(Lupri,'(A,I15)') 'l_X=',l_X
  if (lX > Zero) then
    write(LuPri,'(A)') 'This seems to be an integer overflow!'
    call Cho_RWord2Byte(lX,Byte,Unt)
    write(LuPri,'(A,1P,D15.6,A,D15.6,1X,A,A)') 'In double precision, lX=',lX,' words (',Byte,Unt,')'
  end if
  irc = 1
  return
end if
#endif
call mma_allocate(X,l_X,Label='X')

! allocate vector array
Nai = 0
do iBatch=1,nBatch
  do iSym=1,nSym
    Nai = max(Nai,LnT1am(iSym,iBatch)*nEnrVec(iSym))
  end do
end do
call mma_allocate(V,Nai,Label='V')

! compute energy correction
do q=1,N
  ! init energy for this q
  Eq = Zero
  ! scale grid point by 1/2
  tq = Half*t(q)
  ! scale weight by 2 (only lower blocks of X computed)
  wq = Two*w(q)
  do iSym=1,nSym
    if (nEnrVec(iSym) > 0) then
      ! init X for this symmetry
      bsize = min(Laplace_BlockSize,nEnrVec(iSym))
      nBlock = (nEnrVec(iSym)-1)/bsize+1
      blast = nEnrVec(iSym)-bsize*(nBlock-1)
      lenX = nTri_Elem(nEnrVec(iSym))+(nBlock-1)*nTri_Elem(bsize-1)+nTri_Elem(blast-1)
#     ifdef _DEBUGPRINT_
      if (lenX > l_X) then
        call WarningMessage(2,SecNam//': insufficient X allocation')
        write(u6,'(A,2(1X,I10))') 'lenX,l_X=',lenX,l_X
        call Abend()
      end if
#     endif
      X(1:lenX) = Zero
      do iBatch=1,nBatch
        Nai = LnT1am(iSym,iBatch)
        if (Nai > 0) then
          ! open file, read vectors, close file
          ! do not delete file here - it may be needed later
          ! (because of the loop over q)
          call ChoMP2_OpenB(1,iSym,iBatch)
          iOpt = 2
          l_Tot = Nai*nEnrVec(iSym)
          iAddr = 1
#         ifdef _DEBUGPRINT_
          if (l_Tot > size(V)) then
            call WarningMessage(2,SecNam//': insufficient V allocation')
            write(u6,'(A,2(1X,I10))') 'l_Tot,SIZE(V)=',l_Tot,size(V)
            call Abend()
          end if
#         endif
          call dDAFile(lUnit(iSym,iBatch),iOpt,V,l_Tot,iAddr)
          call ChoMP2_OpenB(2,iSym,iBatch)
          ! scale vectors
          do iVec=1,nEnrVec(iSym)
            ip0 = Nai*(iVec-1)
            do iSymi=1,nSym
              if (LnOcc(iSymi,iBatch) > 0) then
                iSyma = Mul(iSym,iSymi)
                ip1 = ip0+LiT1am(iSyma,iSymi,iBatch)
                do i=0,LnOcc(iSymi,iBatch)-1
                  ii = iFirstS(iSymi,iBatch)+i
                  V(ip1+nVir(iSyma)*i+1:ip1+nVir(iSyma)*(i+1)) = exp(EOcc(iOcc(iSymi)+ii)*tq)* &
                                                                 V(ip1+nVir(iSyma)*i+1:ip1+nVir(iSyma)*(i+1))
                end do
                do a=1,nVir(iSyma)
                  call dScal_(LnOcc(iSymi,iBatch),exp(-EVir(iVir(iSyma)+a)*tq),V(ip1+a),nVir(iSyma))
                end do
              end if
            end do
          end do
          ! loop over vector blocks to compute
          ! X(J,K) += sum_ai L(ai,J)*L(ai,K)*exp(-(e(a)-e(i))*t(q)/2)
          ipX = 1
          do jBlock=1,nBlock
            ipj = 1+Nai*Laplace_BlockSize*(jBlock-1)
            if (jBlock == nBlock) then
              nVecj = nEnrVec(iSym)-Laplace_BlockSize*(nBlock-1)
            else
              nVecj = Laplace_BlockSize
            end if
            do iBlock=jBlock,nBlock
              ipi = 1+Nai*Laplace_BlockSize*(iBlock-1)
              if (iBlock == nBlock) then
                nVeci = nEnrVec(iSym)-Laplace_BlockSize*(nBlock-1)
              else
                nVeci = Laplace_BlockSize
              end if
              call dGEMM_('T','N',nVeci,nVecj,Nai,One,V(ipi),Nai,V(ipj),Nai,One,X(ipX),nVeci)
              ipX = ipX+nVeci*nVecj
            end do
          end do
#         ifdef _DEBUGPRINT_
          if (lenX /= (ipX-1)) then
            call WarningMessage(2,SecNam//': dimension problem [1]')
            write(u6,'(A,I10,A,I10)') 'lenX=',lenX,' ipX-1=',ipX-1
            call Abend()
          end if
#         endif
        end if
      end do
      ! compute energy contribution
      ! Eq += sum_JK [X(J,K)]**2
      ipX = 1
      do jBlock=1,nBlock
        if (jBlock == nBlock) then
          nVecj = nEnrVec(iSym)-Laplace_BlockSize*(nBlock-1)
        else
          nVecj = Laplace_BlockSize
        end if
        do iBlock=jBlock,nBlock
          if (iBlock == nBlock) then
            nVeci = nEnrVec(iSym)-Laplace_BlockSize*(nBlock-1)
          else
            nVeci = Laplace_BlockSize
          end if
          if (iBlock == jBlock) then
            Eq = Eq+Half*dDot_(nVeci*nVecj,X(ipX),1,X(ipX),1)
          else
            Eq = Eq+dDot_(nVeci*nVecj,X(ipX),1,X(ipX),1)
          end if
          ipX = ipX+nVeci*nVecj
        end do
      end do
#     ifdef _DEBUGPRINT_
      if (lenX /= (ipX-1)) then
        call WarningMessage(2,SecNam//': dimension problem [2]')
        write(u6,'(A,I10,A,I10)') 'lenX=',lenX,' ipX-1=',ipX-1
        call Abend()
      end if
#     endif
    end if
  end do
  ! scale Eq
  Eq = -wq*Eq
  ! Accumulate in EMP2
  EMP2 = EMP2+Eq
end do

! deallocations
call mma_deallocate(V)
call mma_deallocate(X)

! delete files if requested
if (Delete) then
  do iBatch=1,nBatch
    do iSym=1,nSym
      call ChoMP2_OpenB(1,iSym,iBatch)
      call ChoMP2_OpenB(3,iSym,iBatch)
    end do
  end do
end if

end subroutine ChoLSOSMP2_Energy_Srt
