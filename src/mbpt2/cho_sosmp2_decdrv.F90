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
! Copyright (C) 2007, Francesco Aquilante                              *
!***********************************************************************

subroutine Cho_SOSmp2_DecDrv(irc,DelOrig,Diag)
! Francesco Aquilante, May 2007.
!
! Purpose: decompose M(ai,bj) = (ai|bj)^2 for use in
!          SOS-MP2 approach.
!
! DelOrig: flag for deleting files with original vectors after
!          decomposition completes.

use Cholesky, only: lBuf, nSym, NumCho, Span
use ChoMP2, only: ChkDecoMP2, Incore, MxQual_Def, MxQualMP2, nMP2Vec, NowSym, nT1am, lUnit_F, OldVec, SpanMP2, ThrMP2, Verbose
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
logical(kind=iwp), intent(in) :: DelOrig
real(kind=wp), intent(in) :: Diag(*)
integer(kind=iwp) :: iAdr, iBin, iClos(2), iOpt, iSym, iTyp, kOffD, lB, Left, lTot, MxQual, nBin, nDim, nInC
real(kind=wp) :: ErrStat(3), RMS, Thr, XMn, XMx
logical(kind=iwp) :: Failed
integer(kind=iwp), allocatable :: iPivot(:), iQual(:)
real(kind=wp), allocatable :: Bin(:), Buf(:), Qual(:)
logical(kind=iwp), parameter :: Restart = .false.
character(len=*), parameter :: SecNam = 'Cho_SOSmp2_DecDrv'
external :: Cho_SOSmp2_Col, ChoMP2_Vec

! Initializations.
! ----------------

irc = 0

if (Verbose) then
  nBin = 18
  call mma_allocate(Bin,nBin,label='Bin')
else
  nBin = 0
end if

do iSym=1,nSym
  nMP2Vec(iSym) = 0
  InCore(iSym) = .false.
end do
NowSym = -999999

if (DelOrig) then
  iClos(1) = 3  ! signals close and delete original vectors
else
  iClos(1) = 2  ! signals close and keep original vectors
end if
iClos(2) = 2    ! signals close and keep result vectors

! Print.
! ------

if (Verbose) then
  write(u6,*)
  call Cho_Head('Cholesky Decomposition of  M(ai,bj) = (ai|bj)^2 for SOS-MP2','=',80,u6)
  write(u6,'(/,1X,A)') 'Configuration of decomposition:'
  write(u6,'(1X,A,1P,D15.6)') 'Threshold: ',ThrMP2
  write(u6,'(1X,A,1P,D15.6)') 'Span     : ',SpanMP2
  if (ChkDecoMP2) then
    write(u6,'(1X,A)') 'Full decomposition check activated.'
  end if
end if

! Start symmetry loop.
! --------------------

kOffD = 1
do iSym=1,nSym

  nDim = nT1am(iSym)
  if ((nDim > 0) .and. (NumCho(iSym) > 0)) then

    if (Verbose .and. (nBin > 0)) then
      Bin(1) = 1.0e2_wp
      do iBin=2,nBin
        Bin(iBin) = Bin(iBin-1)*0.1_wp
      end do
      write(u6,'(//,1X,A,I2,A,I9)') '>>> Cholesky decomposing symmetry block ',iSym,', dimension: ',nDim
      write(u6,'(/,1X,A)') 'Analysis of initial diagonal:'
      call Cho_AnaSize(Diag(kOffD),nDim,Bin,nBin,u6)
    end if

    ! Open files.
    ! -----------

    do iTyp=1,2
      call ChoMP2_OpenF(1,iTyp,iSym)
    end do

    ! Setup decomposition.
    ! --------------------

    NowSym = iSym

    if (MxQualMP2 /= MxQual_Def) then ! user-defined
      MxQual = min(max(MxQualMP2,1),nDim)
    else ! default
      if (nDim > 10) then
        MxQual = max(min(nDim/10,MxQualMP2),1)
      else
        MxQual = max(min(nDim,MxQualMP2),1)
      end if
    end if
#   if !defined (_I8_)
    lTstBuf = (nDim+MxQual)*MxQual
    lTstQua = nDim*(MxQual+1)
    do while (((lTstBuf < 0) .or. (lTstQua < 0)) .and. (MxQual > 0))
      MxQual = MxQual-1
      lTstBuf = (nDim+MxQual)*MxQual
      lTstQua = nDim*(MxQual+1)
    end do
    if (MxQual < 1) then
      write(u6,*) SecNam,': MxQual causes integer overflow!'
      write(u6,*) SecNam,': parameters:'
      write(u6,*) 'Symmetry block: ',iSym
      write(u6,*) 'Dimension     : ',nDim
      write(u6,*) 'MxQual        : ',MxQual
      irc = -99
      call finalize()
      return
    end if
#   endif

    call mma_allocate(Qual,nDim*(MxQual+1),label='Qual')
    call mma_allocate(iQual,MxQual,label='iQual')
    call mma_allocate(iPivot,nDim,label='iPivot')

    call mma_maxDBLE(lB)
    lBuf = min((nDim+MxQual)*MxQual,lB)
    Left = lB-lBuf
    nInC = Left/nDim
    if (nInC >= NumCho(iSym)) then
      InCore(iSym) = .true.
      lTot = nDim*NumCho(iSym)
      call mma_allocate(OldVec,lTot,Label='OldVec')
      iOpt = 2
      iAdr = 1
      call ddaFile(lUnit_F(iSym,1),iOpt,OldVec,lTot,iAdr)
    end if
    call mma_maxDBLE(lBuf)
    call mma_allocate(Buf,lBuf,label='DecBuf')

    ! Decompose this symmetry block.
    ! ------------------------------

    Thr = ThrMP2
    Span = SpanMP2
    call ChoDec(Cho_SOSmp2_Col,ChoMP2_Vec,Restart,Thr,Span,MxQual,Diag(kOffD),Qual,Buf,iPivot,iQual,nDim,lBuf,ErrStat, &
                nMP2Vec(iSym),irc)
    if (irc /= 0) then
      write(u6,*) SecNam,': ChoDec returned ',irc,'   Symmetry block: ',iSym
      call finalize()
      return
    end if
    XMn = ErrStat(1)
    XMx = ErrStat(2)
    RMS = ErrStat(3)
    if (Verbose) then
      write(u6,'(/,1X,A)') '- decomposition completed!'
      write(u6,'(1X,A,I9,A,I9,A)') 'Number of vectors needed: ',nMP2Vec(iSym),' (number of AO vectors: ',NumCho(iSym),')'
      write(u6,'(1X,A)') 'Error statistics for (ai|ai)^2 [min,max,rms]:'
      write(u6,'(1X,1P,3(D15.6,1X))') XMn,XMx,RMS
    end if
    Failed = (abs(Xmn) > Thr) .or. (abs(XMx) > thr) .or. (RMS > Thr)
    if (Failed) then
      if (.not. Verbose) then
        write(u6,'(1X,A)') 'Error statistics for (ai|ai)^2 [min,max,rms]:'
        write(u6,'(1X,1P,3(D15.6,1X))') XMn,XMx,RMS
      end if
      write(u6,*) SecNam,': (ai|bj)^2 decomposition failed!'
      irc = -9999
      call finalize()
      return
    end if

    ! If requested, check decomposition.
    ! ----------------------------------

    if (ChkDecoMP2) then
      write(u6,*)
      write(u6,*) SecNam,': Checking M(ai,bj)=(ai|bj)^2 CD.'
      write(u6,*) 'Symmetry block: ',iSym
      write(u6,*) 'Threshold, Span, MxQual: ',Thr,Span,MxQual
      write(u6,*) 'Error statistics for (ai|ai)^2 [min,max,rms]:'
      write(u6,*) ErrStat(:)
      call Cho_SOSmp2_DecChk(irc,iSym,Qual,nDim,MxQual,Buf,lBuf,ErrStat)
      if (irc /= 0) then
        write(u6,*) SecNam,': ChoMP2_DecChk returned ',irc,'   Symmetry block: ',iSym
        call SysAbendMsg(SecNam,'SOS-MP2 decomposition failed!',' ')
      else
        XMn = ErrStat(1)
        XMx = ErrStat(2)
        RMS = ErrStat(3)
        Failed = Failed .or. (abs(Xmn) > Thr) .or. (abs(XMx) > Thr) .or. (RMS > Thr)
        write(u6,*) 'Error statistics for (ai|bj)^2 [min,max,rms]:'
        write(u6,*) XMn,XMx,RMS
        if (Failed) then
          write(u6,*) '==> DECOMPOSITION FAILURE <=='
          irc = -9999
          call finalize()
          return
        else
          write(u6,*) '==> DECOMPOSITION SUCCESS <=='
        end if
        call xFlush(u6)
      end if
    end if

    ! Free memory.
    ! ------------

    call mma_deallocate(Buf)
    if (InCore(iSym)) call mma_deallocate(OldVec)
    call mma_deallocate(Qual)
    call mma_deallocate(iQual)
    call mma_deallocate(iPivot)

    ! Close (possibly deleting original) files.
    ! -----------------------------------------

    do iTyp=1,2
      call ChoMP2_OpenF(iClos(iTyp),iTyp,iSym)
    end do

    ! Update pointer to diagonal block.
    ! ---------------------------------

    kOffD = kOffD+nT1am(iSym)

  else

    if (Verbose) then
      write(u6,'(//,1X,A,I2,A)') '>>> Symmetry block',iSym,' is empty!'
    end if

  end if

end do

call finalize()

return

contains

subroutine finalize()
  integer(kind=iwp) :: iSym, iTyp
  if (irc /= 0) then ! make sure files are closed before exit
    do iSym=1,nSym
      do iTyp=1,2
        call ChoMP2_OpenF(2,iTyp,iSym)
      end do
    end do
  end if
  if (allocated(Bin)) call mma_deallocate(Bin)
end subroutine finalize

end subroutine Cho_SOSmp2_DecDrv
