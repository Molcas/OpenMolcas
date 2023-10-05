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

use Cholesky, only: lBuf, nSym, NumCho, Span
use ChoMP2, only: ChkDecoMP2, Incore, iOption_MP2CD, lUnit_F, MxQual_Def, MxQualMP2, nMP2Vec, NowSym, nT1am, OldVec, SpanMP2, &
                  ThrMP2, Verbose
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
logical(kind=iwp), intent(in) :: DelOrig
real(kind=wp), intent(in) :: Diag(*)
character(len=*) :: CD_Type
integer(kind=iwp) :: iAdr, iBin, iClos(2), iOpt, IOPTION, ISYM, iTyp, kOffD, lB, LEFT, LERRSTAT, lTot, MxCDVec(8), MxQual, nBin, &
                     nDim, nInc
real(kind=wp) :: THR, XMN, XMX, RMS
logical(kind=iwp) :: Failed, ConventionalCD
character(len=18) :: Option
integer(kind=iwp), allocatable :: iPivot(:), iQual(:)
real(kind=wp), allocatable :: Bin(:), Buf(:), ErrStat(:), Qual(:)
integer(kind=iwp), parameter :: nOption = 2
logical(kind=iwp), parameter :: Restart = .false.
character(len=*), parameter :: SecNam = 'ChoMP2_DecDrv'
external :: ChoMP2_Col, ChoMP2_Vec

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
  write(u6,*) SecNam,': illegal input option (argument CD_Type)'
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
MxCDVec(1:nSym) = nT1Am(1:nSym)

lErrStat = 3
call mma_allocate(ErrStat,lErrStat,Label='ErrStat')
if (Verbose) then
  nBin = 18
  call mma_allocate(Bin,nBin,Label='Bin')
else
  nBin = 0
  call mma_allocate(Bin,1,Label='Bin')
end if

nMP2Vec(1:nSym) = 0
InCore(1:nSym) = .false.
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
  write(u6,*)
  call Cho_Head('Cholesky decomposition of '//Option,'=',80,u6)
  write(u6,'(/,1X,A)') 'Configuration of decomposition:'
  write(u6,'(1X,A,1P,D15.6)') 'Threshold: ',ThrMP2
  write(u6,'(1X,A,1P,D15.6)') 'Span     : ',SpanMP2
  if (ChkDecoMP2) write(u6,'(1X,A)') 'Full decomposition check activated.'
end if

! Start symmetry loop.
! --------------------

kOffD = 1
do iSym=1,nSym

  nDim = nT1am(iSym)
  if ((nDim > 0) .and. (NumCho(iSym) > 0)) then

    ConventionalCD = (MxCDVec(iSym) < 1) .or. (MxCDVec(iSym) >= nDim)
    if (Verbose .and. (nBin > 0)) then
      Bin(1) = 1.0e2_wp
      do iBin=2,nBin
        Bin(iBin) = Bin(iBin-1)*1.0e-1_wp
      end do
      if (ConventionalCD) then
        write(u6,'(//,1X,A,I2,A,I9)') '>>> Conventional Cholesky decomposition of symmetry block ',iSym,', dimension: ',nDim
      else
        write(u6,'(//,1X,A,I2,A,I9)') '>>> MaxVec Cholesky decomposition of symmetry block ',iSym,', dimension: ',nDim
      end if
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
#   ifndef _I8_
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
      call Finish_this()
      return
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
        write(u6,*) SecNam,': ChoDec returned ',irc,'   Symmetry block: ',iSym
        call Finish_this()
        return
      end if
    else
      call ChoDec_MxVec(ChoMP2_Col,ChoMP2_Vec,MxCDVec(iSym),Restart,Thr,Span,MxQual,Diag(kOffD),Qual,Buf,iPivot,iQual,nDim,lBuf, &
                        ErrStat,nMP2Vec(iSym),irc)
      if (irc /= 0) then
        write(u6,*) SecNam,': ChoDec_MxVec returned ',irc,'   Symmetry block: ',iSym
        call Finish_this()
        return
      end if
    end if
    XMn = ErrStat(1)
    XMx = ErrStat(2)
    RMS = ErrStat(3)
    if (Verbose) then
      write(u6,'(/,1X,A)') '- decomposition completed!'
      write(u6,'(1X,A,I9,A,I9,A)') 'Number of vectors needed: ',nMP2Vec(iSym),' (number of AO vectors: ',NumCho(iSym),')'
      if (.not. ConventionalCD) write(u6,'(1X,A,I9)') 'Max. number of vectors allowed: ',MxCDVec(iSym)
      write(u6,'(1X,A)') 'Error statistics for diagonal [min,max,rms]:'
      write(u6,'(1X,1P,3(D15.6,1X))') XMn,XMx,RMS
    end if
    if (ConventionalCD) then
      Failed = (abs(Xmn) > Thr) .or. (abs(XMx) > thr) .or. (RMS > Thr)
    else
      Failed = abs(Xmn) > Thr
    end if
    if (Failed) then
      if (.not. Verbose) then
        write(u6,'(1X,A)') 'Error statistics for diagonal [min,max,rms]:'
        write(u6,'(1X,1P,3(D15.6,1X))') XMn,XMx,RMS
      end if
      write(u6,'(A,A,A,A)') SecNam,': decomposition of ',Option,' failed!'
      irc = -9999
      call Finish_this()
      return
    end if

    ! If requested, check decomposition.
    ! ----------------------------------

    if (ChkDecoMP2) then
      write(u6,*)
      write(u6,'(A,A,A)') SecNam,': Checking decomposition of ',Option
      write(u6,*) 'Symmetry block: ',iSym
      write(u6,*) 'Threshold, Span, MxQual: ',Thr,Span,MxQual
      write(u6,*) 'Error statistics for diagonal [min,max,rms]:'
      write(u6,*) ErrStat(:)
      call ChoMP2_DecChk(irc,iSym,Qual,nDim,MxQual,Buf,lBuf,ErrStat)
      if (irc /= 0) then
        if (irc == -123456) then
          write(u6,*) ' -- Sorry, full decomposition check not yet implemented --'
          irc = 0
        else
          write(u6,*) SecNam,': ChoMP2_DecChk returned ',irc,'   Symmetry block: ',iSym
          call SysAbendMsg(SecNam,'decomposition failed!',' ')
        end if
      else
        XMn = ErrStat(1)
        XMx = ErrStat(2)
        RMS = ErrStat(3)
        write(u6,'(A,A,A)') 'Error statistics for ',Option,' [min,max,rms]:'
        write(u6,*) XMn,XMx,RMS
        Failed = Failed .or. (abs(Xmn) > Thr) .or. (abs(XMx) > Thr) .or. (RMS > Thr)
        if (ConventionalCD) then
          if (Failed) then
            write(u6,*) '==> DECOMPOSITION FAILURE <=='
            irc = -9999
            call Finish_this()
            return
          else
            write(u6,*) '==> DECOMPOSITION SUCCESS <=='
          end if
        else
          if (Failed) then
            write(u6,*) '==> DECOMPOSITION SUCCESS <== (by definition)'
          else
            write(u6,*) '==> DECOMPOSITION SUCCESS <=='
          end if
        end if
        call xFlush(u6)
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

    if (Verbose) write(u6,'(//,1X,A,I2,A)') '>>> Symmetry block',iSym,' is empty!'

  end if

end do

call Finish_this()

contains

subroutine Finish_this()

  integer(kind=iwp) :: iSym, iTyp

  if (irc /= 0) then ! make sure files are closed before exit
    do iSym=1,nSym
      do iTyp=1,2
        call ChoMP2_OpenF(2,iTyp,iSym)
      end do
    end do
  end if

  call mma_deallocate(Bin)
  call mma_deallocate(ErrStat)

end subroutine Finish_this

end subroutine ChoMP2_DecDrv
