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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************
!  Cho_X_Init
!
!> @brief
!>   Initialize Cholesky vector information for external use.
!> @author Thomas Bondo Pedersen
!>
!> @details
!> This routine reads and processes the information
!> stored on the runfile/restart files by the Cholesky
!> decomposition utility.
!> This routine is also used for setting up the environment in
!> density fitting (DF or RI) runs.
!> Index arrays are allocated and initialized.
!> All information is stored in the module Cholesky
!>
!> \p BufFrac is the fraction of total available memory that will be
!> allocated as Cholesky vector buffer. For example, \p BufFrac = ``0.35``
!> implies that 35% of the total available memory (after the
!> allocations of this routine) will be allocated for Cholesky
!> vectors. The vectors will be read into the buffer as part of
!> the initialization. The reading routines can then be used as
!> usual; the buffer is automatically taken care of by the reading
!> routines. If \p BufFrac is less than or equal to zero, no buffer
!> will be used.
!>
!> Return codes:
!>
!> - \p irc = ``-2``: Local DF not implemented here
!> - \p irc = ``-1``: Cholesky flag not found on runfile
!> - \p irc =  ``0``: initialization success
!> - \p irc =  ``1``: runfile info corrupted
!> - \p irc =  ``2``: restart file info corrupted
!> - \p irc =  ``3``: inconsistent include file(s) detected (typically an internal error/bug)
!> - \p irc =  ``4``: error in parallel setup
!>
!> @note
!> The two-electron repulsion integrals must have been decomposed by Seward.
!>
!> @param[out] irc     Return code
!> @param[in]  BufFrac Fraction of memory to be used as buffer
!***********************************************************************

subroutine Cho_X_Init(irc,BufFrac)

use Index_Functions, only: nTri_Elem
use Para_Info, only: Is_Real_Par
use Cholesky, only: BkmThr, BkmVec, Cho_AdrVec, Cho_Fake_Par, Cho_IOVec, ChoIniCheck, Cho_Real_Par, iBas, iBasSh, iiBstRSh, &
                    iiBstRSh_Hidden, IndRed, IndRed_Hidden, IndRSh, IndRSh_Hidden, IPRINT, iRS2F, iShlSO, iSOShl, iSP2F, LuCho, &
                    LuMap, LuPri, LuRed, LuRst, MaxRed, MaxVec, mmBstRT, Mx2Sh, MxOrSh, MySP, N1_VecRd, N2_VecRd, n_MySP, nBas, &
                    nBasSh, nBasT, nBstSh, nCol_BkmThr, nCol_BkmVec, nDGM_call, nDimRS, nnBstRSh, nnBstRSh_Hidden, nnBstR, &
                    nnBstRT, nnShl, nnShl_SP, nnShl_Tot, nRow_BkmThr, nRow_BkmVec, nShell, nSym, nSys_call, NumCho, NumChT, &
                    RUN_EXTERNAL, RUN_MODE
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(in) :: BufFrac
integer(kind=iwp) :: ChoIsIni, iErr, ijShl, iLoc, iRed, iShl, iSym, jShl, l, Numij
real(kind=wp) :: Frac
logical(kind=iwp) :: DidCholesky, DoDummy, FirstCall = .true., isDF
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: is1CCD, l_Max
real(kind=wp) :: Byte
character(len=2) :: Unt
#endif
integer(kind=iwp), allocatable :: BkmDim(:)
character(len=*), parameter :: SecNam = 'Cho_X_Init'

#ifdef _DEBUGPRINT_
call mma_maxDBLE(l_Max)
call Cho_Word2Byte(l_Max,8,Byte,Unt)
write(u6,*) '>>>>> Available memory on entry to ',SecNam,': ',l_Max,' = ',Byte,Unt
#endif

! Check that this is a Cholesky run.
! ----------------------------------

call DecideOnCholesky(DidCholesky)
if (.not. DidCholesky) then
  call Finish_this(-1)
  return
end if

! Check if already initialized.
! -----------------------------

if (FirstCall) then ! it cannot be already done
  FirstCall = .false.
else ! might be already done
  call Get_iScalar('ChoIni',ChoIsIni)
  if (ChoIsIni == ChoIniCheck) then ! already done
    irc = 0
    call Finish_this(0)
    return
  end if
end if

! Check if this is density fitting (DF).
! --------------------------------------

call DecideOnDF(isDF)

! Define entries in Cholesky module.
! ----------------------------------

call Cho_X_SetInc(irc)
if (irc /= 0) then  ! include file inconsistency detected
  call Finish_this(3)
  return
end if

! Set parallel info (picked up from para_info).
! ---------------------------------------------

CHO_FAKE_PAR = .false.
Cho_Real_Par = Is_Real_Par() .and. (.not. CHO_FAKE_PAR)

! Define n_MySP.
! --------------

n_MySP = 0

! Set run mode to "external".
! ---------------------------

Run_Mode = Run_External

! Set output unit used by the decomposition core routines.
! --------------------------------------------------------

LuPri = u6

! Set print level to -5 (ensuring that Cho_X_Checkdiag will
! print information, if called). All other routines will be
! silent.
! ---------------------------------------------------------

iPrint = -5

! Get nSym: the number of irreps.
! -------------------------------

call Get_iScalar('nSym',nSym)
if ((nSym < 1) .or. (nSym > 8)) then
  write(u6,*) SecNam,': nSym out of bounds: ',nSym
  call Finish_this(1)
  return
end if

! Get Cho_AdrVec: addressing of vector files (1=WA, 2=DA).
! --------------------------------------------------------

call Get_iScalar('ChoVec Address',Cho_AdrVec)

! Open files with red. set and vector info.
! -----------------------------------------

LURED = 0
LUCHO(1:NSYM) = 0
LURST = 0
LUMAP = 0
call Cho_OpenVR(1,2)

! Set vector I/O model etc.
! -------------------------

Cho_IOVec = 3
N1_VecRd = 2
N2_VecRd = 3
nSys_Call = 0
nDGM_Call = 0

! Get info (derived) from the runfile.
! nBas  : #basis functions in each irrep
! iSOShl: shell index for each basis function (SO)
! NumCho: #Cholesky vectors in each irrep
! MaxVec: max. element in NumCho (used to allocate InfVec)
! --------------------------------------------------------

call Get_iArray('nBas',nBas,nSym)
iBas(1) = 0
nBasT = nBas(1)
do iSym=2,nSym
  iBas(iSym) = nBasT
  nBasT = nBasT+nBas(iSym)
end do
if (nBasT < 1) then
  write(u6,*) SecNam,': nBasT out of bounds: ',nBasT
  call Finish_this(1)
  return
end if
call mma_allocate(iSOShl,nBasT,Label='iSOShl')
call Get_iArray('ISOSHL',iSOShl,nBasT)

call Get_iArray('NumCho',NumCho,nSym)
NumChT = sum(NumCho(1:nSym))
MaxVec = NumCho(1)
do iSym=2,nSym
  MaxVec = max(MaxVec,NumCho(iSym))
end do

! Read Cholesky restart file. Allocate InfRed and InfVec.
! nShell: #shells
! nnShl_Tot : total #shell pairs
! nnShl : #shell pairs contributing in diagonal
! MaxRed: #reduced sets (used to allocate InfRed)
! InfRed: InfRed(i) is the disk address of reduced set i
! InfVec: InfVec(i,1,iSym) is the parent index of vector i in irrep
!                          iSym in first reduced set
!         InfVec(i,2,iSym) is the reduced set of this vector
!         InfVec(i,3,iSym) is the disk address for reading this
!                          vector
!         InfVec(i,4,iSym) is the WA disk address of this vector.
!         InfVec(i,5,iSym) in a parallel run is the global index
!                          of the i-th vector in iSym for this node
!                          (in serial: InfVec(i,5,iSym) = i)
! Note: InfVec(i,5,iSym) is treated in a quite dirty way here:
!   for DF:
!       not defined in serial, thus it will be defined by the call
!       to Cho_X_DefineInfVec_5.
!       defined in parallel (simply read from disk and not modified
!       by Cho_X_DefineInfVec_5).
!   for Cholesky:
!       not defined in serial, thus it will be defined by the call
!       to Cho_X_DefineInfVec_5.
!       not defined in parallel, thus it will be defined by the call
!       to Cho_X_DefineInfVec_5. It is later modified by
!       Cho_X_Init_Par.
! ------------------------------------------------------------------

ierr = 0
call Cho_X_RdRst(ierr)
if (ierr /= 0) then
  call Finish_this(2)
  return
end if
nnShl_Tot = nTri_Elem(nShell)
call Cho_X_DefineInfVec_5(isDF)

! nnShl_SP makes it possible to use function Cho_F2SP.
! (Stored in module Cholesky)
! ----------------------------------------------------

nnShl_SP = nnShl

! Allocate and initialize index arrays.
! -------------------------------------

call mma_allocate(iiBstRSh_Hidden,nSym,nnShl,3,Label='iiBstRSh_Hidden')
iiBstRSh => iiBstRSh_Hidden
call mma_allocate(nnBstRSh_Hidden,nSym,nnShl,3,Label='nnBstRSh_Hidden')
nnBstRSh => nnBstRSh_Hidden
call Cho_RstD_GetInd1()
mmBstRT = nnBstRT(1)

call mma_allocate(IndRed_Hidden,nnBstRT(1),3,Label='IndRed_Hidden')
IndRed => IndRed_Hidden
call mma_allocate(IndRSh_Hidden,nnBstRT(1),Label='IndRSh_Hidden')
IndRSh => IndRSh_Hidden
call Cho_RstD_GetInd2()

call mma_allocate(iSP2F,nnShl,Label='iSP2F')
call Cho_RstD_GetInd3(iSP2F,size(iSP2F))

! Allocate and read bookmarks (if available on runfile).
! ------------------------------------------------------

if (isDF) then
  nRow_BkmVec = 0
  nCol_BkmVec = 0
  nRow_BkmThr = 0
  nCol_BkmThr = 0
else
  l = 4
  call mma_allocate(BkmDim,l,Label='BkmDim')
  call Get_iArray('Cholesky BkmDim',BkmDim,l)
  nRow_BkmVec = BkmDim(1)
  nCol_BkmVec = BkmDim(2)
  nRow_BkmThr = BkmDim(3)
  nCol_BkmThr = BkmDim(4)
  call mma_deallocate(BkmDim)
  if ((nRow_BkmVec > 0) .and. (nCol_BkmVec > 0) .and. (nRow_BkmThr > 0) .and. (nCol_BkmThr > 0)) then
    call mma_allocate(BkmVec,nRow_BkmVec,nCol_BkmVec,Label='BkmVec')
    call Get_iArray('Cholesky BkmVec',BkmVec,size(BkmVec))
    call mma_allocate(BkmThr,nRow_BkmThr,nCol_BkmThr,Label='BkmThr')
    call Get_dArray('Cholesky BkmThr',BkmThr,size(BkmThr))
  else
    nRow_BkmVec = 0
    nCol_BkmVec = 0
    nRow_BkmThr = 0
    nCol_BkmThr = 0
  end if
end if

! mySP is set because a few core routines may use it.
! After the decomposition is done, it must be a trivial mapping and
! the user (programmer) should not worry about it at all.
! -----------------------------------------------------------------

call mma_allocate(MySP,nnShl,Label='MySP')
do ijShl=1,nnShl
  MySP(ijShl) = ijShl
end do

! Copy reduced set 1 to location 2.
! ---------------------------------

call Cho_RSCopy(1,2)

! Get dimensions of reduced sets.
! -------------------------------

call mma_allocate(nDimRS,nSym,MaxRed,Label='nDimRS')
nDimRS(:,1) = nnBstR(1:nSym,1)
iLoc = 3
do iRed=2,MaxRed
  call Cho_GetRed(iRed,iLoc,.false.)
  call Cho_SetRedInd(iLoc)
  nDimRS(:,iRed) = nnBstR(1:nSym,iLoc)
end do

! Copy reduced set 1 to location 3.
! ---------------------------------

call Cho_RSCopy(1,3)

! Derive:
! nBasSh: #basis functions in sym. block of shell.
! nBstSh: #basis functions in each shell.
! Mx2Sh : max. shell pair dimension.
! iShlSO: index of SO within its shell.
! ------------------------------------------------

call mma_allocate(iBasSh,nSym,nShell,Label='iBasSh')
call mma_allocate(nBasSh,nSym,nShell,Label='nBasSh')
call mma_allocate(nBstSh,nShell,Label='nBstSh')
call Cho_SetSh(iBasSh,nBasSh,nBstSh,iBas,nBas,iSOShl,nSym,nShell,nBasT)

MxOrSh = nBstSh(1)
do iShl=2,nShell
  MxOrSh = max(MxOrSh,nBstSh(iShl))
end do

Mx2Sh = 0
do ijShl=1,nnShl
  call Cho_InvPck(iSP2F(ijShl),iShl,jShl,.true.)
  if (iShl == jShl) then
    Numij = nTri_Elem(nBstSh(iShl))
  else
    Numij = nBstSh(iShl)*nBstSh(jShl)
  end if
  Mx2Sh = max(Mx2Sh,Numij)
end do

call mma_allocate(iShlSO,nBasT,Label='iShlSO')
call Cho_SetSh2(iShlSO,iSOShl,nBstSh,nBasT,nShell)

! Allocate and compute mapping RS1->Full.
! ---------------------------------------

call mma_allocate(iRS2F,2,nnBstRT(1),Label='iRS2F')
call Cho_RSTOF(iRS2F,2,nnBstRT(1),1)

! Allocate integer scratch array used for Cholesky vector I/O.
! ------------------------------------------------------------

DoDummy = .not. ((Cho_IOVec == 1) .or. (Cho_IOVec == 2) .or. (Cho_IOVec == 3) .or. (Cho_IOVec == 4))
call Cho_Allo_iScr(DoDummy)

! Setup for parallel runs.
! ------------------------

call Cho_X_Init_Par(irc,isDF)
if (irc /= 0) then
  call Finish_this(4)
  return
end if

#ifdef _DEBUGPRINT_
! Debug: test bookmarks.
! Note that 1C-CD flag must be available on runfile
! (make sure _DEBUGPRINT_ is defined also in Cho_Final().
! -------------------------------------------------------

if (allocated(BkmVec) .and. allocated(BkmThr)) then
  call Get_iScalar('1C-CD',is1CCD)
  call Cho_TestBookmark(irc,.true.,is1CCD == 1)
  if (irc /= 0) call Cho_Quit('Bookmark test failed!',104)
end if
#endif

! Allocate and initialize (i.e. read vectors) vector buffer.
! ----------------------------------------------------------

Frac = min(max(BufFrac,Zero),One)
call Cho_VecBuf_Init(Frac,nnBstR(1,1)) ! allocate
call Cho_VecBuf_Ini2() ! read

! Normal exit point.
! ------------------

ChoIsIni = ChoIniCheck
call Put_iScalar('ChoIni',ChoIsIni)
irc = 0
call Finish_this(0)

contains

subroutine Finish_this(num)

  integer(kind=iwp), intent(in) :: num

  ! Error branches.
  ! ===============
  select case (num)
    case (-1)
      ! Cholesky flag not found on runfile
      irc = -1
      write(u6,'(//,A,A,//)') SecNam,': two-electron integrals not Cholesky decomposed!'
    case (1)
      ! Bad info obtained from runfile
      irc = 1
      write(u6,'(//,A,A,//)') SecNam,': WARNING: error reading runfile!'
    case (2)
      ! restart info corrupted
      irc = 2
      write(u6,'(//,A,A)') SecNam,': WARNING: error reading restart info!'
      write(u6,'(A,A,I6,//)') SecNam,': return code from read:',ierr
    case (3)
      ! include file inconsistency detected
      irc = 3
      write(u6,'(//,A,A,//)') SecNam,': WARNING: include file inconsistency detected!'
    case (4)
      ! error in parallel setup
      irc = 4
      write(u6,'(//,A,A,//)') SecNam,': WARNING: error in parallel setup!'
    case default
  end select

  ! Return.
  ! =======

# ifdef _DEBUGPRINT_
  call mma_maxDBLE(l_Max)
  call Cho_Word2Byte(l_Max,8,Byte,Unt)
  write(u6,*) '>>>>> Available memory on exit from ',SecNam,': ',l_Max,' = ',Byte,Unt
  call xFlush(u6)
# endif

end subroutine Finish_this

end subroutine Cho_X_Init
