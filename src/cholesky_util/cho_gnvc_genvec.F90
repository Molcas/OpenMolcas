!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Cho_GnVc_GenVec(Diag,xInt,lInt,nVecRS,iVecRS,RS2RS,mSym,mPass,iPass1,NumPass)
!
! Purpose: generate Cholesky vectors from raw integral columns.

use Data_Structures, only: Alloc1DiArray_Type
use Cholesky, only: iiBstR, IndRed, INF_PASS, INF_PROGRESS, InfVec, iOff_Col, IPRINT, LuPri, MaxRed, nDimRS, nnBstR, nnZTot, &
                    nQual, nSym, NumCho, NumChT, tDecom, XnPass
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lInt, mSym, mPass, nVecRS(mSym,mPass), iVecRS(mSym,mPass), iPass1, NumPass
real(kind=wp), intent(inout) :: Diag(*), xInt(lInt)
type(Alloc1DiArray_Type) :: RS2RS(8)
integer(kind=iwp) :: i, iAB, ii, iOff1(8), iOff2(8), iP, ip_Scr, iPass, iPass2, irc, iSym, iV, iVec, iVec1, iVecT, jAB, jj, jPass, &
                     jV, jVec, jVec0, kAB, kOff, kOff0, kOff1, kOff2, l_VecTmp, l_Wrk, lAB, LenLin, lOff0, lTot, MxSubtr, nAB, &
                     nBin, nConv, Need, nNeg, nNegT, nPass, NumCho_OLD(8), NumVec
real(kind=wp) :: Bin1, C1, C2, Fac, olDiag, Step, W1, W2, XC, xM, xMax, xMin
real(kind=wp), allocatable :: VecTmp(:), Wrk(:)
character(len=*), parameter :: SecNam = 'Cho_GnVc_GenVec'

! Check input.
! ------------

if (NumPass < 1) return

if (mSym /= nSym) call Cho_Quit('Input error [1] in '//SecNam,103)

if (iPass1 < 1) call Cho_Quit('Input error [2] in '//SecNam,103)

iPass2 = iPass1+NumPass-1
if (iPass2 > mPass) call Cho_Quit('Input error [3] in '//SecNam,103)

nPass = XnPass
if (mPass /= nPass) call Cho_Quit('Input error [4] in '//SecNam,103)

NumVec = 0
do iPass=iPass1,iPass2
  NumVec = NumVec+sum(nVecRS(1:nSym,iPass))
end do
if (NumVec < 1) return ! exit

! Subtract previous vectors.
! --------------------------

call mma_maxDBLE(l_Wrk)
call mma_allocate(Wrk,l_Wrk,Label='Wrk')
do iSym=1,nSym
  kOff = iOff_Col(iSym)+1
  call Cho_Subtr(xInt(kOff),Wrk,l_Wrk,iSym)
end do
call mma_deallocate(Wrk)

! Initialize vector generation.
! -----------------------------

iOff1(1:nSym) = iOff_Col(1:nSym)+1
iOff2(1:nSym) = iOff_Col(1:nSym)+1

l_VecTmp = 0
do iPass=iPass1,iPass2
  do iSym=1,nSym
    Need = nnBstR(iSym,2)*nVecRS(iSym,iPass)
    l_VecTmp = max(l_VecTmp,Need)
  end do
end do
MxSubtr = 0
do iPass=iPass1,iPass2
  do iSym=1,nSym
    nAB = sum(nVecRS(iSym,iPass+1:iPass2))
    Need = nAB*nVecRS(iSym,iPass)
    MxSubtr = max(MxSubtr,Need)
  end do
end do
l_VecTmp = max(l_VecTmp,MxSubtr)
call mma_allocate(VecTmp,l_VecTmp,Label='VecTmp')

! Copy reduced set iPass1 to location 3.
! --------------------------------------

irc = 0
call Cho_X_RSCopy(irc,2,3)
if (irc /= 0) then
  write(Lupri,*) SecNam,': Cho_X_RSCopy returned ',irc
  call Cho_Quit('Error termination in '//SecNam,104)
end if

! Decomposition pass loop.
! ------------------------

do iPass=iPass1,iPass2

  ! Print header.
  ! -------------

  LenLin = 0 ! to avoid compiler warnings
  if (iPrint >= INF_PROGRESS) then
    call Cho_Head(SecNam//': Generation of Vectors from Map','=',80,Lupri)
    write(Lupri,'(/,A,I5)') 'Integral pass number',iPass
    write(Lupri,'(A,8I8)') '#Cholesky vec.: ',(NumCho(iSym),iSym=1,nSym)
    write(Lupri,'(A,8I8)') '#qualified    : ',(nVecRS(iSym,iPass),iSym=1,nSym)
    write(Lupri,'(A,8I8)') 'Current  dim. : ',(nnBstR(iSym,3),iSym=1,nSym)
    write(Lupri,'(A,8I8)') 'Original dim. : ',(nnBstR(iSym,1),iSym=1,nSym)
    write(Lupri,'(/,A,/,A)') '           #Vectors             Treated Diagonal', &
                             'Sym.     Sym.     Total     Index     Before      After   Conv. Neg.   New Max'
    LenLin = 79
    write(Lupri,'(80A)') ('-',i=1,LenLin)
    call XFlush(Lupri)
    NumCho_OLD(1:nSym) = NumCho(1:nSym)
  else if (iPrint >= INF_PASS) then
    write(Lupri,'(/,A,I5)') 'Integral pass number',iPass
    write(LUPRI,'(A,8I8)') '#Cholesky vec.: ',(NumCho(iSym),iSym=1,nSym)
    write(LUPRI,'(A,8I8)') '#qualified    : ',(nVecRS(iSym,iPass),iSym=1,nSym)
    call XFlush(Lupri)
    NumCho_OLD(1:nSym) = NumCho(1:nSym)
  end if

  ! Zero entries in integral matrix that are not part of this
  ! reduced set.
  ! ---------------------------------------------------------

  if (iPass > iPass1) then
    do iSym=1,nSym
      do iV=1,nVecRS(iSym,iPass)
        lTot = nnBstR(iSym,2)
        VecTmp(1:lTot) = Zero
        lOff0 = iOff1(iSym)+nnBstR(iSym,2)*(iV-1)-1
        do lAB=1,nnBstR(iSym,3)
          jAB = IndRed(iiBstR(iSym,3)+lAB,3)-iiBstR(iSym,1)
          kAB = RS2RS(jAB)%A(iSym)
          VecTmp(kAB) = xInt(lOff0+kAB)
        end do
        XInt(lOff0+1:lOff0+lTot) = VecTmp(1:lTot)
      end do
    end do
  end if

  ! Write reduced set info for this pass (index arrays are stored
  ! at location 3).
  ! -------------------------------------------------------------

  call Cho_P_PutRed(iPass,3)

  ! Start symmetry loop.
  ! --------------------

  do iSym=1,nSym

    if (nVecRS(iSym,iPass) >= 1) then

      ! Generate vectors for this pass and symmetry.
      ! --------------------------------------------

      do iV=1,nVecRS(iSym,iPass)

        kOff0 = iOff1(iSym)+nnBstR(iSym,2)*(iV-1)-1

        iVec = iVecRS(iSym,iPass)+iV-1
        iAB = InfVec(iVec,1,iSym) ! addr in 1st red. set

        XC = Diag(iAB)
        if (abs(Diag(iAB)) > 1.0e-14_wp) then ! TODO/FIXME
          Fac = One/sqrt(abs(Diag(iAB)))
        else
          Fac = 1.0e7_wp
        end if
        xInt(kOff0:kOff0+nnBstR(iSym,2)) = Fac*xInt(kOff0:kOff0+nnBstR(iSym,2))

        do i=1,nnBstR(iSym,2)
          ii = iiBstR(iSym,2)+i
          jj = IndRed(ii,2)
          if (Diag(jj) == Zero) xInt(kOff0+i) = Zero
        end do

        do i=1,nnBstR(iSym,2)
          ii = iiBstR(iSym,2)+i
          jj = IndRed(ii,2)
          Diag(jj) = Diag(jj)-xInt(kOff0+i)*xInt(kOff0+i)
        end do

        olDiag = Diag(iAB)
        Diag(iAB) = Zero
        call Cho_ChkDia(Diag,iSym,xMin,xMax,xM,nNegT,nNeg,nConv)
        nNZTot = nNZTot+nNeg

        do jV=iV+1,nVecRS(iSym,iPass)
          jVec = iVecRS(iSym,iPass)+jV-1
          jAB = InfVec(jVec,1,iSym)
          kOff2 = iOff1(iSym)+nnBstR(iSym,2)*(jV-1)-1
          Fac = -xInt(kOff0+RS2RS(jAB-iiBstR(iSym,1))%A(iSym))
          xInt(kOff2+1:kOff2+nnBstR(iSym,2)) = xInt(kOff2+1:kOff2+nnBstR(iSym,2))+Fac*xInt(kOff0+1:kOff0+nnBstR(iSym,2))
        end do

        call Cho_SetVecInf(iVec,iSym,iAB,iPass,3)

        if (iPrint >= INF_PROGRESS) then
          iVecT = NumChT+iV
          write(Lupri,'(I3,3(1X,I9),2(1X,D11.3),2(1X,I4),1X,D11.3)') iSym,iVec,iVecT,iAB,XC,olDiag,nConv,nNeg,xM
        end if

      end do

      ! Subtract contributions to later vectors.
      ! ----------------------------------------

      nAB = nQual(iSym)+sum(nVecRS(iSym,iPass1:iPass))
      if (nAB > 0) then
        ip_Scr = 1
        iP = iPass
        jVec0 = -1
        do while ((iP < iPass2) .and. (jVec0 < 0))
          iP = iP+1
          jVec0 = iVecRS(iSym,iP)-1
        end do
        if (jVec0 < 0) call Cho_Quit('jVec0 < 0 in '//SecNam,103) ! should never happen
        do iV=1,nVecRS(iSym,iPass)
          kOff1 = ip_Scr+nAB*(iV-1)-1
          kOff2 = iOff1(iSym)+nnBstR(iSym,2)*(iV-1)-1
          do iAB=1,nAB
            jVec = jVec0+iAB
            jAB = InfVec(jVec,1,iSym)
            kAB = RS2RS(jAB-iiBstR(iSym,1))%A(iSym)
            VecTmp(kOff1+iAB) = xInt(kOff2+kAB)
          end do
        end do
        kOff1 = iOff1(iSym)
        kOff2 = iOff1(iSym)+nnBstR(iSym,2)*nVecRS(iSym,iPass)
        call DGEMM_('N','T',nnBstR(iSym,2),nAB,nVecRS(iSym,iPass),-One,xInt(kOff1),nnBstR(iSym,2),VecTmp,nAB,One,xInt(kOff2), &
                    nnBstR(iSym,2))
      end if

      ! Reorder vectors to appropriate reduced set.
      ! Skipped for iPass1, as they are already in correct storage.
      ! -----------------------------------------------------------

      if (iPass > iPass1) then
        lTot = nnBstR(iSym,2)*nVecRS(iSym,iPass)
        VecTmp(1:lTot) = xInt(iOff1(iSym):iOff1(iSym)+lTot-1)
        do iV=1,nVecRS(iSym,iPass)
          kOff0 = iOff2(iSym)+nnBstR(iSym,3)*(iV-1)-1
          lOff0 = nnBstR(iSym,2)*(iV-1)
          do kAB=1,nnBstR(iSym,3)
            jAB = IndRed(iiBstR(iSym,3)+kAB,3)
            lAB = RS2RS(jAB-iiBstR(iSym,1))%A(iSym)
            xInt(kOff0+kAB) = VecTmp(lOff0+lAB)
          end do
        end do
      end if

      ! Update vector counters.
      ! -----------------------

      NumCho(iSym) = NumCho(iSym)+nVecRS(iSym,iPass)
      NumChT = NumChT+nVecRS(iSym,iPass)

      ! Update pointer arrays.
      ! iOff1: pointer to integral columns (in xInt).
      ! iOff2: pointer to vectors (also in xInt).
      ! ---------------------------------------------

      iOff1(iSym) = iOff1(iSym)+nnBstR(iSym,2)*nVecRS(iSym,iPass)
      iOff2(iSym) = iOff2(iSym)+nnBstR(iSym,3)*nVecRS(iSym,iPass)
    end if

    ! Cycle point for empty symmetry.
    ! -------------------------------

    if (iPrint >= INF_PROGRESS) call XFlush(Lupri)

  end do ! symmetry

  ! Print.
  ! ------

  if (iPrint >= INF_PROGRESS) then
    NumCho_OLD(1:nSym) = NumCho(1:nSym)-NumCho_OLD(1:nSym)
    write(Lupri,'(80A)') ('-',I=1,LenLin)
    write(Lupri,'(A,8I8)') '#vec. gener.  : ',(NumCho_OLD(iSym),iSym=1,nSym)
    call XFlush(Lupri)
  else if (iPrint >= INF_PASS) then
    NumCho_OLD(1:nSym) = NumCho(1:nSym)-NumCho_OLD(1:nSym)
    write(Lupri,'(A,8I8)') '#vec. gener.  : ',(NumCho_OLD(iSym),iSym=1,nSym)
    call XFlush(Lupri)
  end if

  ! Analyze diagonal.
  ! -----------------

  if (iPrint >= INF_PASS) then
    Bin1 = 1.0e2_wp
    Step = 1.0e-1_wp
    nBin = 18
    call Cho_AnaDia(Diag,Bin1,Step,nBin,.false.)
  end if

  ! Set next (iPass+1) reduced set at location 2.
  ! Reduced set iPass1 is now stored at location 3.
  ! -----------------------------------------------

  call Cho_SetRed(Diag)
  jPass = iPass+1
  call Cho_SetRSDim(nDimRS,nSym,MaxRed,jPass,2)
  if (iPrint >= INF_PASS) then
    call Cho_PrtRed(2)
    call XFlush(Lupri)
  end if

  ! Swap locations so that:
  ! location 2 contains reduced set iPass1 and
  ! location 3 contains next (iPass+1) reduced set.
  ! -----------------------------------------------

  irc = 0
  call Cho_X_RSSwap(irc,2,3)
  if (irc /= 0) then
    write(Lupri,*) SecNam,': Cho_X_RSSwap returned ',irc
    call Cho_Quit('Error termination in '//SecNam,104)
  end if

end do ! integral pass

! Deallocate temporary vector array.
! ----------------------------------

call mma_deallocate(VecTmp)

! Write vectors to disk.
! ----------------------

call CWTime(C1,W1)
do iSym=1,nSym
  NumVec = sum(nVecRS(iSym,iPass1:iPass2))
  if (NumVec > 0) then
    iPass = iPass1
    iVec1 = iVecRS(iSym,iPass)
    do while ((iVec1 < 1) .and. (iPass < iPass2))
      iPass = iPass+1
      iVec1 = iVecRS(iSym,iPass)
    end do
    if (iVec1 < 1) then
      call Cho_Quit('Logical error in '//SecNam,103)
    else
      call Cho_PutVec2(xInt(iOff_Col(iSym)+1),NumVec,iVec1,iSym)
    end if
  end if
end do
call CWTime(C2,W2)
tDecom(1,2) = tDecom(1,2)+C2-C1
tDecom(2,2) = tDecom(2,2)+W2-W1

! Write restart files.
! --------------------

call Cho_P_WrRstC(iPass2)

! Store next (iPass2+1) reduced set at location 2.
! ------------------------------------------------

irc = 0
call Cho_X_RSCopy(irc,3,2)
if (irc /= 0) then
  write(Lupri,*) SecNam,': Cho_X_RSCopy returned ',irc
  call Cho_Quit('Error termination in '//SecNam,104)
end if

end subroutine Cho_GnVc_GenVec
