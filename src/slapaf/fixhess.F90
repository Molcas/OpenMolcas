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

subroutine FixHess(H,nH,iOptC,MF,GNrm,nsAtom,AnalHess,AllowFindTS)

use Index_Functions, only: nTri_Elem
use Slapaf_Info, only: iNeg, GNrm_Threshold, Mode
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Ten, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nH, nsAtom
real(kind=wp), intent(inout) :: H(nH,nH), MF(3,nsAtom)
integer(kind=iwp), intent(inout) :: iOptC
real(kind=wp), intent(in) :: GNrm
logical(kind=iwp), intent(in) :: AnalHess, AllowFindTS
integer(kind=iwp) :: i, ij, iLow, iStatus, iTest, j, jNeg, nRP, NumVal, nVStep
real(kind=wp) :: dRx, Fact, Fix_Val, rLow, rq, SumHii, temp, Test
logical(kind=iwp) :: Corrected, Found
real(kind=wp), allocatable :: EVal(:), FixVal(:), LowVal(:), LowVec(:,:), Rx(:,:), Tmp(:,:), Vect(:)
real(kind=wp), parameter :: HHigh = One, HTh = 1.0e-3_wp, ZTh = 1.0e-12_wp
real(kind=wp), external :: DDot_
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
logical(kind=iwp) :: Too_Small
#endif

#ifdef _DEBUGPRINT_
write(u6,*) 'AnalHess=',AnalHess
call RecPrt('FixHess: H(Start)',' ',H,nH,nH)
Lu = u6
Too_Small = .false.
#endif

Corrected = .false.

call mma_allocate(EVal,nTri_Elem(nH),Label='EVal')

! Copy elements for H

SumHii = Zero
ij = 0
do i=1,nH
  do j=1,i
    ij = ij+1
    EVal(ij) = H(i,j)
  end do
  SumHii = SumHii+H(i,i)
end do

#ifdef _DEBUGPRINT_
write(Lu,*) 'FixHess: SumHii=',SumHii
call RecPrt('FixHess: Hessian',' ',H,nH,nH)
call TriPrt('FixHess: H',' ',EVal,nH)
#endif

! Compute eigenvalues and eigenvectors

! Davidson procedure to compute only the lowest eigenvalues
! For small matrices we can afford to solve it directly
! (tested with NumVal=2 in all cases)
if (nH <= 30) then
  NumVal = nH
else
  NumVal = 2
end if
nVStep = 2
Found = .false.
call mma_allocate(LowVal,NumVal,Label='LowVal')
call mma_allocate(LowVec,nH,NumVal,Label='LowVec')
LowVec(:,:) = Zero
! Stop when the highest eigenvalue found is larger than HTh * 10
do while (.not. Found)
  call Davidson(EVal,nH,NumVal,LowVal,LowVec,iStatus)
# ifdef _DEBUGPRINT_
  call RecPrt(' Eigenvalues',' ',LowVal,1,NumVal)
  call RecPrt(' Eigenvectors',' ',LowVec,nH,NumVal)
# endif
  if (iStatus > 0) call SysWarnMsg('FixHess','Davidson procedure did not converge','')
  if ((LowVal(NumVal) > Ten*HTh) .or. (NumVal >= nH)) then
    Found = .true.
  else
    ! Increase the number of eigenpairs to compute
    call mma_allocate(Tmp,nH,NumVal,Label='Tmp')
    Tmp(:,:) = LowVec(:,:)
    call mma_deallocate(LowVal)
    call mma_deallocate(LowVec)
    ! At some point, start doubling the number
    if (NumVal >= 16) NVStep = NumVal
    nVStep = min(nVStep,nH-NumVal)
    NumVal = NumVal+nVStep
    call mma_allocate(LowVal,NumVal,Label='LowVal')
    LowVal(:) = Zero
    call mma_allocate(LowVec,nH,NumVal,Label='LowVec')
    LowVec(:,:) = Zero
    LowVec(:,1:NumVal-nVStep) = Tmp(:,:)
    call mma_deallocate(Tmp)
  end if
end do
call mma_deallocate(EVal)

! Apply corrections if any ...

call mma_allocate(FixVal,NumVal,Label='FixVal')
#ifdef _DEBUGPRINT_
call RecPrt(' Eigenvalues',' ',LowVal,1,NumVal)
call RecPrt(' Eigenvectors',' ',LowVec,nH,NumVal)
#endif
iNeg(1) = 0
jNeg = 0
rLow = Ten
iLow = 0
! with sorted eigenvalues, jNeg=iNeg, iLow=1
do i=1,NumVal
  temp = LowVal(i)
  FixVal(i) = temp
  if (temp < rlow) then
    rlow = temp
    iLow = i
  end if
  ! No fixes if the Hessian is analytical
  if ((.not. AnalHess) .and. (abs(temp) < HTh)) then
#   ifdef _DEBUGPRINT_
    Too_Small = .true.
#   endif
    Corrected = .true.

    ! For redundant coordinates we will have some
    ! eigenvalues which are zero due to the redundancy.

    if (abs(temp) < ZTh) then
      temp = Zero
      FixVal(i) = Zero
    else
      FixVal(i) = sign(HTh,temp)
    end if
  end if
  if (temp < Zero) then
    iNeg(1) = iNeg(1)+1
    jNeg = i
    if (((.not. AnalHess) .or. btest(iOptC,8)) .and. btest(iOptC,7) .and. (.not. btest(iOptC,12))) then

      ! Change the sign and if just too large reduce the value
      ! to the default of HHigh.

      FixVal(i) = min(HHigh,abs(FixVal(i)))
      Corrected = .true.
    end if
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! If FindTS and we have got a negative eigenvalue signal that
! we are in the TS regime.

if (AllowFindTS) then
  call qpg_darray('TanVec',Found,nRP)
  if (btest(iOptC,12) .and. ((iNeg(1) >= 1) .or. (Mode >= 0)) .and. ((GNrm <= GNrm_Threshold) .or. Found)) then
    if (.not. btest(iOptC,13)) then
      write(u6,*) '**************************'
      write(u6,*) '* Enable TS optimization *'
      write(u6,*) '**************************'
    end if
    iOptC = ibset(iOptC,13)
  end if
  if (btest(iOptC,13)) iOptC = ibclr(iOptC,7)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
if (Too_Small) then
  write(Lu,*)
  write(Lu,*) ' Some too small eigenvalues has been corrected'
  write(Lu,*)
end if
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! MINIMA
!                                                                      *
!***********************************************************************
!                                                                      *
if (btest(iOptC,7)) then

  if ((iNeg(1) /= 0) .and. (.not. btest(iOptC,8))) then
    Corrected = .true.
#   ifdef _DEBUGPRINT_
    write(Lu,*) ' Some negative eigenvalues has been corrected'
    write(Lu,*) 'iNeg=',iNeg(1)
    write(Lu,*)
#   endif
  end if
!                                                                      *
!***********************************************************************
!                                                                      *
! TRANSITION STATE SEARCH
!                                                                      *
!***********************************************************************
!                                                                      *
else if (btest(iOptC,3) .or. btest(iOptC,13)) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! 1 negative eigenvalue

  if (iNeg(1) == 1) then

    if (Mode <= 0) then

      ! Store the eigenvector which we are following

      Mode = jNeg
      call ReacX(LowVec(:,Mode),nH,MF,3*nsAtom)
#     ifdef _DEBUGPRINT_
      write(Lu,'(A,I3)') ' Store Original mode:',Mode
      call RecPrt(' Reaction mode',' ',MF,3,nsAtom)
#     endif

    else

      ! Check that it is the correct eigenvector!

#     ifdef _DEBUGPRINT_
      call RecPrt(' Old Reaction mode',' ',MF,3,nsAtom)
#     endif
      iTest = 0
      Test = Zero
      call mma_allocate(Rx,3,nsAtom,Label='Rx')
      do i=1,NumVal
        call ReacX(LowVec(:,i),nH,Rx,3*nsAtom)
        dRx = sqrt(DDot_(3*nsAtom,Rx,1,Rx,1))
        rq = abs(DDot_(3*nsAtom,MF,1,Rx,1))/dRx
        if (rq > Test) then
          iTest = i
          Test = rq
        end if
        Temp = FixVal(i)
#       ifdef _DEBUGPRINT_
        write(u6,*) '<old|new>,H_new=',rq,Temp
#       endif
      end do
      call mma_deallocate(Rx)

      ! Only iTest and jNeg may be touched
      if (iTest == jNeg) then
        Mode = jNeg
      else
#       ifdef _DEBUGPRINT_
        write(Lu,*) ' Warning: wrong eigenvector has negative eigenvalue.'
#       endif
        ! Keep the old vector if there is significant overlap
        ! Note: there could be a better vector not in the computed set
        if ((.not. AnalHess) .and. (Test > Half)) then
          Mode = iTest
#         ifdef _DEBUGPRINT_
          write(Lu,*) 'Keep old eigenvector!',Mode
#         endif
          FixVal(jNeg) = abs(FixVal(jNeg))
          Corrected = .true.
          ! Prefer the new eigenvector if the Hessian is analytical
          ! or if the best overlap is poor
        else
          Mode = jNeg
#         ifdef _DEBUGPRINT_
          write(Lu,*) 'Take new eigenvector!',Mode
#         endif
        end if
      end if

      FixVal(Mode) = -abs(FixVal(Mode))
      call ReacX(LowVec(:,Mode),nH,MF,3*nsAtom)
#     ifdef _DEBUGPRINT_
      write(Lu,'(A,1X,I3)') ' Store mode:',Mode
      call RecPrt(' New Reaction mode',' ',MF,3,nsAtom)
#     endif

    end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! 0 negative eigenvalues

  else if (iNeg(1) == 0) then

    if (Mode < 1) then

      ! Store the eigenvector which we are following

      Mode = iLow
      call ReacX(LowVec(:,Mode),nH,MF,3*nsAtom)
#     ifdef _DEBUGPRINT_
      write(Lu,'(A,I3)') ' Store Original mode:',Mode
      call RecPrt(' Reaction mode',' ',MF,3,nsAtom)
#     endif

    else

      ! Find the eigenvector with the best overlap

#     ifdef _DEBUGPRINT_
      call RecPrt(' Old Reaction mode',' ',MF,3,nsAtom)
#     endif
      iTest = 0
      Test = Zero
      call mma_allocate(Rx,3,nsAtom,Label='Rx')
      do i=1,NumVal
        call ReacX(LowVec(:,i),nH,Rx,3*nsAtom)
        dRx = sqrt(DDot_(3*nsAtom,Rx,1,Rx,1))
        rq = abs(DDot_(3*nsAtom,MF,1,Rx,1))/dRx
        if (rq > Test) then
          iTest = i
          Test = rq
        end if
        Temp = FixVal(i)
#       ifdef _DEBUGPRINT_
        write(u6,*) '<old|new>,H_new=',rq,Temp
#       endif
      end do
      call mma_deallocate(Rx)

      ! Keep the old vector if there is significant overlap
      ! Note: there could be a better vector not in the computed set
      if (Test > Half) then
        Mode = iTest
        ! Prefer the lowest eigenvector if the best overlap is poor
      else
        Mode = iLow
#       ifdef _DEBUGPRINT_
        write(Lu,*) ' Warning: no good overlap among the computed set of eigenvectors.'
#       endif
      end if

      call ReacX(LowVec(:,Mode),nH,MF,3*nsAtom)
#     ifdef _DEBUGPRINT_
      write(Lu,'(A,1X,I3)') ' Store mode:',Mode
      call RecPrt(' New Reaction mode',' ',MF,3,nsAtom)
#     endif

    end if

    FixVal(Mode) = -Half*abs(FixVal(Mode))
    Corrected = .true.
#   ifdef _DEBUGPRINT_
    write(Lu,'(A,I2,A)') ' No negative eigenvalue, correction: mode ',Mode,' was changed to negative'
    write(Lu,*)
#   endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! 2 or more negative eigenvalues

  else if (iNeg(1) >= 2) then

    if (Mode < 1) then

      ! Store the eigenvector which we are following

      Mode = iLow
      call ReacX(LowVec(:,Mode),nH,MF,3*nsAtom)
#     ifdef _DEBUGPRINT_
      write(Lu,'(A,I3)') ' Store Original mode:',Mode
      call RecPrt(' Reaction mode',' ',MF,3,nsAtom)
#     endif

    else

      ! Find the eigenvector with the best overlap

#     ifdef _DEBUGPRINT_
      call RecPrt(' Old Reaction mode',' ',MF,3,nsAtom)
#     endif
      iTest = 0
      Test = Zero
      call mma_allocate(Rx,3,nsAtom,Label='Rx')
      do i=1,NumVal
        call ReacX(LowVec(:,i),nH,Rx,3*nsAtom)
        dRx = sqrt(DDot_(3*nsAtom,Rx,1,Rx,1))
        rq = abs(DDot_(3*nsAtom,MF,1,Rx,1))/dRx
        if (rq > Test) then
          iTest = i
          Test = rq
        end if
        Temp = FixVal(i)
#       ifdef _DEBUGPRINT_
        write(u6,*) '<old|new>,H_new=',rq,Temp
#       endif
      end do
      call mma_deallocate(Rx)

      ! Keep the old vector if there is significant overlap
      ! Note: there could be a better vector not in the computed set
      if (Test > Half) then
        Mode = iTest
        ! Prefer the lowest eigenvector if the best overlap is poor
      else
        Mode = iLow
#       ifdef _DEBUGPRINT_
        write(Lu,*) ' Warning: no good overlap among the computed set of eigenvectors.'
#       endif
      end if

      call ReacX(LowVec(:,Mode),nH,MF,3*nsAtom)
#     ifdef _DEBUGPRINT_
      write(Lu,'(A,1X,I3)') ' Store mode:',Mode
      call RecPrt(' New Reaction mode',' ',MF,3,nsAtom)
#     endif

    end if

    ! Caution! Negative eigenvalues which are not assigned
    ! to the reaction mode will be increased by two orders of magnitude!

    Fact = Ten**2
    do i=1,NumVal
      Temp = FixVal(i)
      if (i == Mode) then
        FixVal(i) = -abs(Temp)
      else if (Temp < 0) then
        FixVal(i) = abs(Temp)*Fact
      else
        FixVal(i) = abs(Temp)
      end if
    end do
    Corrected = .true.
#   ifdef _DEBUGPRINT_
    write(Lu,'(A,I2,A)') ' Too many negative eigenvalue, correction: mode ',Mode,' was kept'
    write(Lu,*)
#   endif

  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
#ifdef _DEBUGPRINT_
else
  write(u6,*) 'No Hessian massage!'
#endif
end if
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(Lu,*)
write(Lu,*) ' Analysis of the Hessian'
write(Lu,*)
call RecPrt(' Eigenvalues',' ',FixVal,1,NumVal)
call RecPrt(' Eigenvectors',' ',LowVec,nH,NumVal)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Recompute the Hessian if needed

iNeg(2) = iNeg(1)
if (Corrected) then

# ifdef _DEBUGPRINT_
  call RecPrt(' Corrected eigenvalues',' ',FixVal,1,NumVal)
  call RecPrt(' Hessian',' ',H,nH,nH)
# endif

  call mma_allocate(Vect,nH,Label='Vect')
  iNeg(1) = 0
  do i=1,NumVal
    if (FixVal(i) < Zero) iNeg(1) = iNeg(1)+1
    Fix_Val = FixVal(i)-LowVal(i)
    if (abs(Fix_Val) > ZTh) then
      Fix_Val = FixVal(i)+LowVal(i)

      ! H |i>

      call dGeMV_('N',nH,nH,One,H,nH,LowVec(:,i),1,Zero,Vect,1)

      ! H' = (I-|i><i|) H (I-|i><i|) + val_new |i><i|
      !    = H - H|i> <i| - |i> <i|H + (val_old + val_new) |i> <i|
      !
      ! (since |i> is an eigenvector this could be used instead):
      ! H' = H + (val_new - val_old) |i> <i|

      do j=1,nH
        H(:,j) = H(:,j)-Vect(:)*LowVec(j,i)-Vect(j)*LowVec(:,i)+Fix_Val*LowVec(:,i)*LowVec(j,i)
      end do
    end if
  end do
  call mma_deallocate(Vect)

end if

#ifdef _DEBUGPRINT_
call RecPrt('FixHess: Hessian',' ',H,nH,nH)
#endif

call mma_deallocate(FixVal)
call mma_deallocate(LowVal)
call mma_deallocate(LowVec)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end
