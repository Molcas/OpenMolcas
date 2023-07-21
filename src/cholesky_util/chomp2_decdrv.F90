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
! Copyright (C) 2004,2008, Thomas Bondo Pedersen                       *
!***********************************************************************

subroutine ChoMP2_DecDrv(irc,DelOrig,Diag,CD_Type)
!
! Thomas Bondo Pedersen, Dec. 2004.
! * Amplitude extension, Thomas Bondo Pedersen, Dec. 2007/Jan. 2008.
!
! Purpose: decompose (ai|bj) integrals or
!          amplitudes [(ai|bj)/e(a)-e(i)+e(b)-e(j)]
!          for use in Cholesky MP2 program.
!
! Arguments:
! irc...... OUT: return code - if non-zero, decomposition failed.
!           Caller MUST check this!
! DelOrig.. INP: flag for deleting files with original vectors after
!           decomposition completes.
! Diag..... INP: integral diagonal (ai|ai) or
!           amplitude diagonal (ai|ai)/2[e(a)-e(i)]
! CD_Type.. INP: string
!           'Integrals'  - integral decomposition
!           'Amplitudes' - amplitude decomposition
!
! Other input such as orbital energies are read from mbpt2 include
! files.

use ChoMP2, only: OldVec
use ChoMP2_dec, only: Incore, NowSym, iOption_MP2CD
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer irc
logical DelOrig
real*8 Diag(*)
character(len=*) CD_Type
external ChoMP2_Col, ChoMP2_Vec
integer :: IOPTION, ISYM, LERRSTAT, nBin, kOffD, nDim, iBin, iTyp, MxQual, LEFT, lB, nInc, lTot, iOpt, iAdr
real*8 :: THR, XMN, XMX, RMS
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
character(len=6), parameter :: ThisNm = 'DecDrv'
character(len=13), parameter :: SecNam = 'ChoMP2_DecDrv'
logical, parameter :: Restart = .false.
logical Failed, ConventionalCD
integer, parameter :: nOption = 2
character(len=18) Option
integer iClos(2)
integer MxCDVec(8)
real*8, allocatable :: ErrStat(:), Bin(:), Qual(:), Buf(:)
integer, allocatable :: iQual(:), iPivot(:)

! Initializations.
! ----------------

irc = 0

#ifdef _DEBUGPRINT_
ChkDecoMP2 = .true.
#endif

if (len(CD_Type) < 1) then ! input error
  iOption = 0
  Option = 'Unknown !?!?!     '
else
  if ((CD_Type(1:1) == 'i') .or. (CD_Type(1:1) == 'I')) then
    iOption = 1
    Option = '(ai|bj) integrals '
  else if ((CD_Type(1:1) == 'a') .or. (CD_Type(1:1) == 'A')) then
    iOption = 2
    Option = 'MP2 amplitudes    '
  else
    iOption = nOption+1
    Option = 'Unknown !?!?!     '
  end if
end if
if ((iOption < 1) .or. (iOption > nOption)) then
  irc = -98
  write(6,*) SecNam,': illegal input option (argument CD_Type)'
  return
end if
iOption_MP2CD = iOption  ! copy to include file chomp2_dec.fh
!-TBP:
! Frankie,
! I use the array MxCDVec(iSym) to decide whether the decomposition
! of a given symmetry block is conventional or "MaxVec":
!ConventionalCD = (MxCDVec(iSym) < 1) .or. (MxCDVec(iSym) >= nT1Am(iSym))
! You may want to calculate the MxCDVec values somewhere else and store
! them in an include-file (say, chomp2_dec.fh).
! For now, I simply define the array here:
do iSym=1,nSym
  MxCDVec(iSym) = nT1Am(iSym)
end do

lErrStat = 3
call mma_allocate(ErrStat,lErrStat,Label='ErrStat')
if (Verbose) then
  nBin = 18
  call mma_allocate(Bin,nBin,Label='Bin')
else
  nBin = 0
  call mma_allocate(Bin,1,Label='Bin')
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
iClos(2) = 2     ! signals close and keep result vectors

! Print.
! ------

if (Verbose) then
  write(6,*)
  call Cho_Head('Cholesky decomposition of '//Option,'=',80,6)
  write(6,'(/,1X,A)') 'Configuration of decomposition:'
  write(6,'(1X,A,1P,D15.6)') 'Threshold: ',ThrMP2
  write(6,'(1X,A,1P,D15.6)') 'Span     : ',SpanMP2
  if (ChkDecoMP2) write(6,'(1X,A)') 'Full decomposition check activated.'
end if

! Start symmetry loop.
! --------------------

kOffD = 1
do iSym=1,nSym

  nDim = nT1am(iSym)
  if ((nDim > 0) .and. (NumCho(iSym) > 0)) then

    ConventionalCD = (MxCDVec(iSym) < 1) .or. (MxCDVec(iSym) >= nDim)
    if (Verbose .and. (nBin > 0)) then
      Bin(1) = 1.0d2
      do iBin=2,nBin
        Bin(iBin) = Bin(iBin-1)*1.0D-1
      end do
      if (ConventionalCD) then
        write(6,'(//,1X,A,I2,A,I9)') '>>> Conventional Cholesky decomposition of symmetry block ',iSym,', dimension: ',nDim
      else
        write(6,'(//,1X,A,I2,A,I9)') '>>> MaxVec Cholesky decomposition of symmetry block ',iSym,', dimension: ',nDim
      end if
      write(6,'(/,1X,A)') 'Analysis of initial diagonal:'
      call Cho_AnaSize(Diag(kOffD),nDim,Bin,nBin,6)
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
#   ifndef _I8_
    lTstBuf = (nDim+MxQual)*MxQual
    lTstQua = nDim*(MxQual+1)
    do while (((lTstBuf < 0) .or. (lTstQua < 0)) .and. (MxQual > 0))
      MxQual = MxQual-1
      lTstBuf = (nDim+MxQual)*MxQual
      lTstQua = nDim*(MxQual+1)
    end do
    if (MxQual < 1) then
      write(6,*) SecNam,': MxQual causes integer overflow!'
      write(6,*) SecNam,': parameters:'
      write(6,*) 'Symmetry block: ',iSym
      write(6,*) 'Dimension     : ',nDim
      write(6,*) 'MxQual        : ',MxQual
      irc = -99
      Go To 1 ! exit
    end if
#   endif

    call mma_allocate(Qual,nDim*(MxQual+1),Label='Qual')
    call mma_allocate(iQual,MxQual,Label='iQual')
    call mma_allocate(iPivot,nDim,Label='iPivot')

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
    call mma_allocate(Buf,lBuf,Label='Buf')

    ! Decompose this symmetry block.
    ! ------------------------------

    Thr = ThrMP2
    Span = SpanMP2
    if (ConventionalCD) then
      call ChoDec(ChoMP2_Col,ChoMP2_Vec,Restart,Thr,Span,MxQual,Diag(kOffD),Qual,Buf,iPivot,iQual,nDim,lBuf,ErrStat,nMP2Vec(iSym), &
                  irc)
      if (irc /= 0) then
        write(6,*) SecNam,': ChoDec returned ',irc,'   Symmetry block: ',iSym
        Go To 1 ! exit...
      end if
    else
      call ChoDec_MxVec(ChoMP2_Col,ChoMP2_Vec,MxCDVec(iSym),Restart,Thr,Span,MxQual,Diag(kOffD),Qual,Buf,iPivot,iQual,nDim,lBuf, &
                        ErrStat,nMP2Vec(iSym),irc)
      if (irc /= 0) then
        write(6,*) SecNam,': ChoDec_MxVec returned ',irc,'   Symmetry block: ',iSym
        Go To 1 ! exit...
      end if
    end if
    XMn = ErrStat(1)
    XMx = ErrStat(2)
    RMS = ErrStat(3)
    if (Verbose) then
      write(6,'(/,1X,A)') '- decomposition completed!'
      write(6,'(1X,A,I9,A,I9,A)') 'Number of vectors needed: ',nMP2Vec(iSym),' (number of AO vectors: ',NumCho(iSym),')'
      if (.not. ConventionalCD) write(6,'(1X,A,I9)') 'Max. number of vectors allowed: ',MxCDVec(iSym)
      write(6,'(1X,A)') 'Error statistics for diagonal [min,max,rms]:'
      write(6,'(1X,1P,3(D15.6,1X))') XMn,XMx,RMS
    end if
    if (ConventionalCD) then
      Failed = (abs(Xmn) > Thr) .or. (abs(XMx) > thr) .or. (RMS > Thr)
    else
      Failed = abs(Xmn) > Thr
    end if
    if (Failed) then
      if (.not. Verbose) then
        write(6,'(1X,A)') 'Error statistics for diagonal [min,max,rms]:'
        write(6,'(1X,1P,3(D15.6,1X))') XMn,XMx,RMS
      end if
      write(6,'(A,A,A,A)') SecNam,': decomposition of ',Option,' failed!'
      irc = -9999
      Go To 1 ! exit
    end if

    ! If requested, check decomposition.
    ! ----------------------------------

    if (ChkDecoMP2) then
      write(6,*)
      write(6,'(A,A,A)') SecNam,': Checking decomposition of ',Option
      write(6,*) 'Symmetry block: ',iSym
      write(6,*) 'Threshold, Span, MxQual: ',Thr,Span,MxQual
      write(6,*) 'Error statistics for diagonal [min,max,rms]:'
      write(6,*) ErrStat(:)
      call ChoMP2_DecChk(irc,iSym,Qual,nDim,MxQual,Buf,lBuf,ErrStat)
      if (irc /= 0) then
        if (irc == -123456) then
          write(6,*) ' -- Sorry, full decomposition check not yet implemented --'
          irc = 0
        else
          write(6,*) SecNam,': ChoMP2_DecChk returned ',irc,'   Symmetry block: ',iSym
          call ChoMP2_Quit(SecNam,'decomposition failed!',' ')
        end if
      else
        XMn = ErrStat(1)
        XMx = ErrStat(2)
        RMS = ErrStat(3)
        write(6,'(A,A,A)') 'Error statistics for ',Option,' [min,max,rms]:'
        write(6,*) XMn,XMx,RMS
        Failed = Failed .or. (abs(Xmn) > Thr) .or. (abs(XMx) > Thr) .or. (RMS > Thr)
        if (ConventionalCD) then
          if (Failed) then
            write(6,*) '==> DECOMPOSITION FAILURE <=='
            irc = -9999
            Go To 1 ! exit
          else
            write(6,*) '==> DECOMPOSITION SUCCESS <=='
          end if
        else
          if (Failed) then
            write(6,*) '==> DECOMPOSITION SUCCESS <== (by definition)'
          else
            write(6,*) '==> DECOMPOSITION SUCCESS <=='
          end if
        end if
        call xFlush(6)
      end if
    end if

    ! Free memory.
    ! ------------

    call mma_deallocate(Buf)
    if (InCore(iSym)) call mma_deallocate(OldVec)
    call mma_deallocate(iPivot)
    call mma_deallocate(iQual)
    call mma_deallocate(Qual)

    ! Close (possibly deleting original) files.
    ! -----------------------------------------

    do iTyp=1,2
      call ChoMP2_OpenF(iClos(iTyp),iTyp,iSym)
    end do

    ! Update pointer to diagonal block.
    ! ---------------------------------

    kOffD = kOffD+nT1am(iSym)

  else

    if (Verbose) write(6,'(//,1X,A,I2,A)') '>>> Symmetry block',iSym,' is empty!'

  end if

end do

1 continue
if (irc /= 0) then ! make sure files are closed before exit
  do iSym=1,nSym
    do iTyp=1,2
      call ChoMP2_OpenF(2,iTyp,iSym)
    end do
  end do
end if

call mma_deallocate(Bin)
call mma_deallocate(ErrStat)

end subroutine ChoMP2_DecDrv
