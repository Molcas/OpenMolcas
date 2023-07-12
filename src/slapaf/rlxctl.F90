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

subroutine RlxCtl(iStop)
!***********************************************************************
!     Program for determination of the new molecular geometry          *
!***********************************************************************

use Slapaf_Info, only: BMx, BSet, Coor, Cubic, Cx, dqInt, E_Delta, Fallback, Free_Slapaf, GNrm, HSet, HUpMet, iNeg, iRef, &
                       isFalcon, iter, Lbl, lCtoF, lNmHss, lTherm, mRowH, mTROld, mTtAtm, MxItr, nDimBC, NmIter, Numerical, nWNdw, &
                       PrQ, qInt, Request_Alaska, Request_RASSI, Shift, UpMeth, User_Def
use Chkpnt, only: Chkpnt_close, Chkpnt_open, Chkpnt_update
use kriging_mod, only: Kriging, nspAI
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iStop
#include "print.fh"
integer(kind=iwp) :: i, iErr, iOff, j, kIter, Lu, LuSpool, mInt, mIntEff, nGB, nHQ, nHX, nHX2, nIntCoor, nKtB, nQQ
real(kind=wp) :: rDum(1), ThrGrd
logical(kind=iwp) :: Do_ESPF, Error, Found, GoOn, Just_Frequencies, NewCarDone
character :: Step_trunc
real(kind=wp), allocatable :: GB(:), HQ(:), HX(:), KtB(:)
integer(kind=iwp), external :: AixRm

!                                                                      *
!***********************************************************************
!                                                                      *
Lu = u6
Just_Frequencies = .false.
NewCarDone = .false.
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Process the input

LuSpool = 21
call SpoolInp(LuSpool)

call RdCtl_Slapaf(LuSpool,.false.)
mInt = nDimBC-mTROld

call Close_LuSpool(LuSpool)

call Chkpnt_open()
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
if (Request_Alaska .or. Request_RASSI) then

  ! Alaska/RASSI only
  iStop = 3
  call Free_Slapaf()
  return
else if (isFalcon) then
  iStop = 1
  call Free_Slapaf()
  return
end if
!                                                                      *
!***********************************************************************
!                                                                      *
PrQ = .not. Request_Alaska
!                                                                      *
!***********************************************************************
!                                                                      *
if (lCtoF .and. PrQ) call Def_CtoF(.false.)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the Wilson B-matrices, which describe the transformations
! between internal and Cartesian coordinates. Values of the
! Internal coordinates are computed too.

HSet = .true.
BSet = .true.
kIter = iter

! Compute number of steps for numerical differentiation

NmIter = 1
! Numerical only for some rows
if (allocated(mRowH)) NmIter = size(mRowH)+1
if (lNmHss) NmIter = 2*mInt+1   ! Full numerical
if (Cubic) NmIter = 2*mInt**2+1 ! Full cubic

if (lTherm .and. (iter == 1)) call Put_dArray('Initial Coordinates',Coor,size(Coor))

! Fix the definition of internal during numerical differentiation
if (lNmHss .and. (iter < NmIter) .and. (iter /= 1)) nPrint(122) = 5

! Do not overwrite numerical Hessian
if ((lNmHss .or. allocated(mRowH)) .and. ((iter > NmIter) .or. (iter < NmIter))) HSet = .false.

! Set logical to indicate status during numerical differentiation
Numerical = lNmHss .and. (iter <= NmIter) .and. (iter /= 1)

if (Numerical) nWndw = NmIter
iRef = 0
call BMtrx(size(Coor,2),Coor,iter,mTtAtm,nWndw)
nQQ = size(qInt,1)

nPrint(30) = nPrint(30)-1

call Put_dArray('BMtrx',BMx,size(Coor)*nQQ)
call Put_iScalar('No of Internal coordinates',nQQ)
!                                                                      *
!***********************************************************************
!                                                                      *
call Reset_ThrGrd(Iter,mTtAtm,ThrGrd)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Compute the norm of the Cartesian gradient.

call G_Nrm(nQQ,GNrm,iter,dqInt,mIntEff)
if (nPrint(116) >= 6) call ListU(Lu,Lbl,dqInt,nQQ,iter)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Accumulate gradient for complete or partial numerical
! differentiation of the Hessian.

if ((allocated(mRowH) .and. (iter < NmIter)) .or. (lNmHss .and. (iter < NmIter))) then

  if (allocated(mRowH) .and. (iter < NmIter)) then

    !------------------------------------------------------------------*
    !    I) Update geometry for selected numerical differentiation.    *
    !------------------------------------------------------------------*

    call Freq1()
    UpMeth = 'RowH  '
  else

    !------------------------------------------------------------------*
    !    II) Update geometry for full numerical differentiation.       *
    !------------------------------------------------------------------*

    call NwShft()
    UpMeth = 'NumHss'
  end if

  call MxLbls(nQQ,dqInt(:,iter),Shift(:,iter),Lbl)
  iNeg(:) = -99
  HUpMet = ' None '
  nPrint(116) = nPrint(116)-3
  nPrint(52) = nPrint(52)-1  ! Status
  nPrint(53) = nPrint(53)-1
  nPrint(54) = nPrint(54)-1
  write(u6,*) ' Accumulate the gradient for selected numerical differentiation.'
  write(u6,'(1x,i5,1x,a,1x,i5)') iter,'of',NmIter
  E_Delta = Zero
  Step_Trunc = ' '
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! Compute updated geometry in Internal coordinates

  Step_Trunc = ' '
  E_Delta = zero
  if (allocated(mRowH) .or. lNmHss) kIter = iter-(NmIter-1)
# ifdef UNIT_MM
  call Init_UpdMask(nInter)
# endif

  ! Update geometry

  if (Kriging .and. (Iter >= nspAI)) then
    call Update_Kriging(Step_Trunc,nWndw)
    if ((Step_Trunc == '#') .and. Fallback) then
      call Update_sl(Step_Trunc,nWndw/2,kIter)
    else
      NewCarDone = .true.
    end if
  else
    call Update_sl(Step_Trunc,nWndw,kIter)
  end if

# ifdef UNIT_MM
  call Free_UpdMask()
# endif
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Transform the new internal coordinates to Cartesians
! (if not already done by Kriging)

if (NewCarDone) then
  Coor(:,:) = Cx(:,:,Iter+1)
else
  PrQ = .false.
  Error = .false.
  iRef = 0
  call NewCar(Iter,size(Coor,2),Coor,mTtAtm,Error)
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! If this is a ESPF QM/MM job, the link atom coordinates are updated

Do_ESPF = .false.
call DecideOnESPF(Do_ESPF)
if (Do_ESPF) then
  call LA_Morok(size(Coor,2),Coor,2)
  Cx(:,:,Iter+1) = Coor(:,:)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Adjust some print levels

if ((lNmHss .or. allocated(mRowH)) .and. (iter == NmIter)) then

  ! If only frequencies no more output

  nPrint(21) = 5 ! Hessian already printed calling Update_sl
  if (kIter > MxItr) then
    Just_Frequencies = .true.
    nPrint(116) = nPrint(116)-3
    nPrint(52) = nPrint(52)-1
    nPrint(53) = nPrint(53)-1
    nPrint(54) = nPrint(54)-1
  end if
end if

! Fix correct reference structure in case of Numerical Hessian optimization.

if ((lNmHss .or. allocated(mRowH)) .and. (kIter == 1)) Cx(:,:,iter) = Cx(:,:,1)

! Print statistics and check on convergence

GoOn = (lNmHss .and. (iter < NmIter)) .or. (allocated(mRowH) .and. (iter < NmIter))
Numerical = (lNmHss .or. allocated(mRowH)) .and. (iter <= NmIter)
call Convrg(iter,kIter,nQQ,iStop,MxItr,mIntEff,mTtAtm,GoOn,Step_Trunc,Just_Frequencies)

!***********************************************************************
!                                                                      *
!                           EPILOGUE                                   *
!                                                                      *
!***********************************************************************

! Write information to files

Numerical = (lNmHss .or. allocated(mRowH)) .and. (iter <= NmIter)

call DstInf(iStop,Just_Frequencies)

if (lCtoF) call Def_CtoF(.true.)
if ((.not. User_Def) .and. ((lNmHss .and. (iter >= NmIter)) .or. (.not. lNmHss))) call cp_SpcInt()

! After a numerical frequencies calculation, restore the original
! runfile, but save the useful data (gradient and Hessian)

if (lNmHss .and. (iter >= NmIter)) then
  call f_Inquire('RUNBACK',Found)
  if (Found) then
    ! Read data
    nGB = size(Coor)
    call mma_allocate(GB,nGB,Label='GB')
    call Get_dArray_chk('GRAD',GB,nGB)

    call Qpg_dArray('Hss_X',Found,nHX)
    call mma_allocate(HX,nHX,Label='HX')
    call Get_dArray('Hss_X',HX,nHX)

    call Qpg_dArray('Hss_Q',Found,nHQ)
    call mma_allocate(HQ,nHQ,Label='HQ')
    call Get_dArray('Hss_Q',HQ,nHQ)

    call Qpg_dArray('KtB',Found,nKtB)
    call mma_allocate(KtB,nKtB,Label='KtB')
    call Get_dArray('KtB',KtB,nKtB)

    call Get_iScalar('No of Internal coordinates',nIntCoor)
    ! Write data in backup file
    call NameRun('RUNBACK')
    call Put_dArray('GRAD',GB,nGB)
    call Put_dArray('Hss_X',HX,nHX)
    call Put_dArray('Hss_Q',HQ,nHQ)
    call Put_dArray('Hss_upd',rdum,0)
    call Put_dArray('Hess',HQ,nHQ)
    call Put_dArray('KtB',KtB,nKtB)
    call Put_iScalar('No of Internal coordinates',nIntCoor)
    ! Pretend the Hessian is analytical
    nHX2 = int(sqrt(real(nHX,kind=wp)))
    iOff = 0
    do i=1,nHX2
      do j=1,i
        iOff = iOff+1
        HX(iOff) = HX((i-1)*nHX2+j)
      end do
    end do
#   ifdef _DEBUGPRINT_
    call TriPrt('AnalHess',' ',HX,nHX2)
#   endif

    call Put_AnalHess(HX,iOff)
    call NameRun('#Pop')

    call mma_deallocate(GB)
    call mma_deallocate(HX)
    call mma_deallocate(HQ)
    call mma_deallocate(KtB)

    ! Restore and remove the backup runfile

    call fCopy('RUNBACK','RUNFILE',iErr)
    if (iErr /= 0) call Abend()
    if (AixRm('RUNBACK') /= 0) call Abend()
  end if
end if

! Remove the GRADS file

call f_Inquire('GRADS',Found)
if (Found) then
  if (AixRm('GRADS') /= 0) call Abend()
end if

call Chkpnt_update()
call Chkpnt_close()
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocate memory

call Free_Slapaf()

! Terminate the calculations.

return

end subroutine RlxCtl
