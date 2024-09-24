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
! Copyright (C) 1994,2017, Roland Lindh                                *
!               1995, Per-Olof Widmark                                 *
!               1995, Markus P. Fuelscher                              *
!               1995, Piotr Borowski                                   *
!               1995, Martin Schuetz                                   *
!***********************************************************************

!#define _DEBUGPRINT_
!#define _NEW_
subroutine DIIS_x(nD,CInter,nCI,QNRStp,Ind)
!***********************************************************************
!                                                                      *
!     purpose: Accelerate convergence using DIIS method                *
!                                                                      *
!     modified by:                                                     *
!     P.O. Widmark, M.P. Fuelscher, P. Borowski & M.Schuetz            *
!     University of Lund, Sweden, 1995                                 *
!                                                                      *
!     Derived from code for c1- and c2-DIIS as implemented by          *
!     R. Lindh in Slapaf and SCF in 1994.                              *
!                                                                      *
!***********************************************************************

use InfSO, only: IterSO, Energy
use InfSCF, only: TimFld, mOV, kOptim, Iter, C1DIIS, AccCon, Iter_Start, kOV
use Constants, only: Zero, One, Ten
#ifdef _NEW_
use Constants, only: Half
#endif
use MxDM, only: MxOptm
use stdalloc, only: mma_allocate, mma_deallocate
use LnkLst, only: LLx

implicit none
integer nCI, nD
real*8 CInter(nCI,nD)
integer Ind(MxOptm)
logical QNRstp
! Define local variables
real*8, dimension(:,:), allocatable :: EVector, Bij
real*8, dimension(:), allocatable :: EValue, Err1, Err2, Scratch
!real*8, dimension(:), allocatable :: Err3, Err4
real*8 GDiis(MxOptm+1), BijTri(MxOptm*(MxOptm+1)/2)
real*8 EMax, Fact, ee2, ee1, E_Min_G, Dummy, Alpha, B11
real*8 :: E_Min = Zero
integer iVec, jVec, kVec, nBij, nFound
integer :: i, j
!integer :: iPos
integer :: ipBst, ij, iErr, iDiag, iDum
real*8 :: cpu1, cpu2
real*8 :: tim1, tim2, tim3
logical :: Case1 = .false., Case2 = .false., Case3 = .false.
! threshold to determine numerical imbalance.
real*8 :: delta = 1.0D-4
real*8 :: delta_E = 1.0D-4
! Factor for checking that the diagonal B elements decline.
real*8, parameter :: Fact_Decline = 15.0d0
real*8 :: ThrCff = Ten
#ifdef _NOT_USED_
real*8 :: f1 = Half, f2 = One/Half
#endif
real*8 :: c2, Bii_Min, DD, DD1
real*8, external :: DDot_
character(len=80) Text, Fmt
#ifdef _DEBUGPRINT_
real*8 cDotV
#endif
interface
  subroutine OptClc_X(CInter,nCI,nD,Array,mOV,Ind,MxOptm,kOptim,kOV,LL,DD)
    implicit none
    integer nCI, nD, mOV, MxOptm, kOptim, kOV(2), LL
    real*8 CInter(nCI,nD), Array(mOV)
    integer Ind(MxOptm)
    real*8, optional :: DD
  end subroutine OptClc_X
end interface

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

call Timing(Cpu1,Tim1,Tim2,Tim3)
100 continue
Case1 = .false.
Case2 = .false.
Case3 = .false.

! Select from the kOptim last iterations

do i=1,kOptim
  Ind(i) = iter-kOptim+i
end do
#ifdef _DEBUGPRINT_
!write(6,*) 'Iter, Iter_Start=',Iter,Iter_Start
write(6,*) 'kOptim=',kOptim
write(6,*) 'Ind(i):',(Ind(i),i=1,kOptim)
#endif

! The following piece of code computes the DIIS coeffs
! with error vectors chosen as the grds (DIIS only)
! or as delta=-Hinv*grd (QNR/DIIS)

! Allocate memory for error vectors (gradient or delta)

call mma_allocate(Err1,mOV,Label='Err1')
call mma_allocate(Err2,mOV,Label='Err2')
!call mma_allocate(Err3,mOV,Label='Err3')
!call mma_allocate(Err4,mOV,Label='Err4')
nBij = kOptim+1
call mma_allocate(Bij,nBij,nBij)
call FZero(Bij,nBij**2)

! Compute norms, <e_i|e_j>

#ifdef _DEBUGPRINT_
write(6,*) 'kOV(:)=',kOV
write(6,*) 'mOV   =',mOV
call RecPrt('Energy',' ',Energy,1,iter)
#endif
E_Min_G = Zero
E_Min = Zero
Bii_min = 1.0D+99
do i=1,kOptim
  call ErrV(mOV,Ind(i),QNRStp,Err1)
  !if (QNRStp) call ErrV(mOV,Ind(i),.false.,Err3)
# ifdef _DEBUGPRINT_
  call NrmClc(Err1,mOV,'Diis  ','Err(i) ')
# endif
  do j=1,i-1

    call ErrV(mOV,Ind(j),QNRStp,Err2)
    !if (QNRStp) call ErrV(mOV,Ind(j),.false.,Err4)
#   ifdef _DEBUGPRINT_
    call NrmClc(Err2,mOV,'Diis  ','Err(j)  ')
#   endif
    if (QNRStp) then
#     ifdef _NEW_
      !Bij(i,j) = Half*dble(nD)*(DDot_(mOV,Err1,1,Err4,1)+DDot_(mOV,Err3,1,Err2,1))
#     else
      Bij(i,j) = dble(nD)*DDot_(mOV,Err1,1,Err2,1)
#     endif
    else
      Bij(i,j) = dble(nD)*DDot_(mOV,Err1,1,Err2,1)
    end if
    Bij(j,i) = Bij(i,j)
  end do
  if (QNRStp) then
#   ifdef _NEW_
    !Bij(i,i) = dble(nD)*DDot_(mOV,Err1,1,Err3,1)
#   else
    Bij(i,i) = dble(nD)*DDot_(mOV,Err1,1,Err1,1)
#   endif
  else
    Bij(i,i) = dble(nD)*DDot_(mOV,Err1,1,Err1,1)
  end if
  E_min = min(E_min,Energy(Ind(i)))
  if (Bij(i,i) < Bii_Min) then
    E_min_G = Energy(Ind(i))
    Bii_Min = Bij(i,i)
  end if
end do

i = kOptim
! Monitor if the sequence of norms of the error vectors and their
! corresponding energies are consistent with a single convex
! potential energy minimum.

! Case 1
! Matrix elements are just too large. This is probably due to
! that the BFGS update is ill-conditioned.
! Case1 = ((Bii_Min > One) .and. (kOptim > 1))
Case1 = (((Bij(i,i) > One) .or. (Bii_Min > One)) .and. (kOptim > 1) .and. (IterSO > 1))

! Case 2
! Check if we are sliding off a shoulder, that is, we have a
! lowering of the energy while the norm of the error vector
! increase.
Case2 = ((Bij(i,i) > Bii_Min) .and. (Energy(Ind(i))+delta_E < E_Min_G) .and. (kOptim > 1))

! Case 3
! Check if elements are in decending order
Case3 = .false.
do i=1,kOptim-1
  if (Bij(i,i) < 1.0D-6) cycle
  if (Fact_Decline*Bij(i,i) < Bij(i+1,i+1)) Case3 = .true.
end do
if (Energy(Ind(i)) >= E_min) Case3 = .false.

if (qNRStp .and. (Case1 .or. Case2 .or. Case3)) then
# ifdef _DEBUGPRINT_
  write(6,*) 'Case1=',Case1
  write(6,*) 'Case2=',Case2
  write(6,*) 'Case3=',Case3
  write(6,*) '   RESETTING kOptim!!!!'
  write(6,*) '   Calculation of the norms in Diis :'
  Fmt = '(6(G0.12,2x))'
  Text = 'B-matrix squared in Diis :'
  call RecPrt(Text,Fmt,Bij,nBij,nBij)
  write(6,'(A,2(G0.9,2x))') 'Bij(i,i),      Bii_Min=',Bij(i,i),Bii_Min
  write(6,'(A,2F16.6)') 'Energy(Ind(i)),E_Min_G=',Energy(Ind(i)),E_Min_G
  write(6,*) Energy(Ind(i)),E_Min_g
  write(6,*)

# endif
  ! Rest the depth of the DIIS and the BFGS update.
  if (Case3) then
    write(6,*) 'DIIS_X: Resetting kOptim!'
    write(6,*) '        Caused by inconsistent B matrix values.'
  else if (Case1) then
    write(6,*) 'DIIS_X: Resetting BFGS depth!'
    write(6,*) '        Too large B matrix values.'
  else if (Case2) then
    write(6,*) 'DIIS_X: Resetting kOptim!'
    write(6,*) '        Caused by energies and gradients which are inconsistent with a convex energy functional.'
  end if
  if (Case1) then
    ! The BFGS update is probably to blame. Reset the update depth.
    Iter_Start = Iter
    IterSO = 1
  else
    !if (Case2) then
    write(6,*) 'kOptim=',kOptim,'-> kOptim=',1
    kOptim = 1
    Iter_Start = Iter
    IterSO = 1
    !else
    !  write(6,*) 'kOptim=',kOptim,'-> kOptim=',kOptim-1
    !  kOptim = kOptim-1
    !  Iter_Start = Iter_Start+1
    !  IterSO = IterSO-1
    !end if
  end if
  call mma_deallocate(Err2)
  call mma_deallocate(Err1)
  call mma_deallocate(Bij)
  Go To 100
end if

if (kOptim /= 1) then
  do i=1,kOptim-1
    if (delta*sqrt(Bij(i,i)) > sqrt(Bij(kOptim,kOptim))) then
      write(6,*) 'DIIS_X: Reduction of the subspace dimension due to numerical imbalance of the values in the B-Matrix'
      write(6,*) 'kOptim=',kOptim,'-> kOptim=',kOptim-1
      kOptim = kOptim-1
      Iter_Start = Iter_Start+1
      IterSO = IterSO-1
      call mma_deallocate(Err2)
      call mma_deallocate(Err1)
      call mma_deallocate(Bij)
      Go To 100
    end if
  end do
end if

! Deallocate memory for error vectors & gradient
!call mma_deallocate(Err4)
!call mma_deallocate(Err3)
call mma_deallocate(Err2)
call mma_deallocate(Err1)

#ifdef _DEBUGPRINT_
write(6,*) '   Calculation of the norms in Diis :'
Fmt = '(6f16.8)'
Text = 'B-matrix squared in Diis :'
call RecPrt(Text,Fmt,Bij,nBij,nBij)
write(6,*)
write(6,*)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Here, the DIIS coeffs for DIIS only or QNR/DIIS are
! computed, either with C1DIIS or C2DIIS
! (-> stored in vector CInter)
!                                                                      *
!***********************************************************************
!                                                                      *
if (.not. c1Diis) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !     C2DIIS case                                                    *
  !                                                                    *
  !       References:                                                  *
  !       H. Sellers, Int. J. Quantum Chem. 45, 31-41(1993).           *
  !       doi:10.1002/qua.560450106                                    *
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (QNRStp) then
    AccCon = 'QNRc2DIIS'
  else
    AccCon = 'c2DIIS   '
  end if

  ! Form a unit eigenvector matrix
  call mma_allocate(EVector,kOptim,kOptim,Label='EVector')
  call mma_allocate(EValue,kOptim,Label='EValue')

  call dcopy_(kOptim**2,[Zero],0,EVector,1)
  call dcopy_(kOptim,[One],0,EVector,kOptim+1)

  ! Form a triangular B-matrix

  ij = 1
  do i=1,kOptim
    call dcopy_(i,Bij(i,1),nBij,BijTri(ij),1)
    ij = ij+i
  end do

# ifdef _DEBUGPRINT_
  Fmt = '(5g25.15)'
  Text = 'B-matrix before Jacobi :'
  call TriPrt(Text,Fmt,BijTri,kOptim)
  Text = 'EigenVectors before Jacobi :'
  call RecPrt(Text,Fmt,EVector,kOptim,kOptim)
  write(6,*)
  write(6,*)
# endif

  ! Diagonalize B-matrix

  EMax = Zero
  do i=1,kOptim*(kOptim+1)/2
    EMax = max(EMax,abs(BijTri(i)))
  end do
  do i=1,kOptim*(kOptim+1)/2
    if (abs(BijTri(i)) < EMax*1.0D-14) BijTri(i) = Zero
  end do

  call mma_allocate(Scratch,kOptim**2,Label='Scratch')

  Dummy = Zero
  iDum = 0
  call Diag_Driver('V','A','L',kOptim,BijTri,Scratch,kOptim,Dummy,Dummy,iDum,iDum,EValue,EVector,kOptim,1,0,'J',nFound,iErr)

  call mma_deallocate(Scratch)
  call dCopy_(kOptim*(kOptim+1)/2,[Zero],0,BijTri,1)

  iDiag = 0
  do i=1,kOptim
    iDiag = iDiag+i
    BijTri(iDiag) = EValue(i)
  end do

# ifdef _DEBUGPRINT_
  Fmt = '(5g25.15)'
  Text = 'B-matrix after Jacobi :'
  call TriPrt(Text,Fmt,BijTri,kOptim)
  Text = 'EigenValues :'
  call RecPrt(Text,Fmt,EValue,1,kOptim)
  Text = 'EigenVectors :'
  call RecPrt(Text,Fmt,EVector,kOptim,kOptim)
  write(6,*)
  write(6,*)
# endif

  ! Renormalize the eigenvectors to the C1-DIIS format

# ifdef _DEBUGPRINT_
  write(6,*) ' Normalization constants :'
# endif

  do iVec=1,kOptim
    Alpha = Zero
    do i=1,kOptim
      Alpha = Alpha+EVector(i,iVec)
    end do

#   ifdef _DEBUGPRINT_
    Fmt = '(A7,i2,A4,f16.8)'
    write(6,Fmt) ' Alpha(',iVec,') = ',Alpha
#   endif

    EVector(:,iVec) = EVector(:,iVec)/Alpha
  end do

  do kVec=1,kOptim
    ee1 = Zero
    do iVec=1,kOptim
      do jVec=1,kOptim
        ee1 = ee1+EVector(iVec,kVec)*EVector(jVec,kVec)*Bij(iVec,jVec)
      end do
    end do
    !write(6,*) 'EValue(kVec),ee1:',EValue(kVec),ee1
    EValue(kVec) = ee1
  end do

# ifdef _DEBUGPRINT_
  Fmt = '(6es16.8)'
  Text = 'B-matrix after scaling :'
  call TriPrt(Text,Fmt,BijTri,kOptim)
  Text = 'EigenValues after scaling :'
  call RecPrt(Text,Fmt,EValue,1,kOptim)
  Text = 'EigenVectors after scaling :'
  call RecPrt(Text,Fmt,EVector,kOptim,kOptim)
  write(6,*)
  write(6,*)
# endif

  ! Select a vector.
  ee1 = 1.0D+72
  DD1 = 1.0D+72
# ifdef _DEBUGPRINT_
  cDotV = 1.0D+72
# endif
  ipBst = -99999999
  call mma_allocate(Scratch,mOV,label='Scratch')
  do iVec=1,kOptim

    ! Pick up eigenvalue (ee2) and the norm of the eigenvector (c2).

    ee2 = EValue(iVec)
    c2 = DDot_(kOptim,EVector(1,iVec),1,EVector(1,iVec),1)
    if (QNRStp) then
      call dcopy_(kOptim,EVector(:,iVec),1,CInter(:,1),1)
      if (nD == 2) call DCopy_(nCI,CInter(:,1),1,CInter(:,2),1)
      call OptClc_x(CInter,nCI,nD,Scratch,mOV,Ind,MxOptm,kOptim,kOV,LLx,DD)
    else
      DD = Zero
    end if
#   ifdef _DEBUGPRINT_
    write(6,*) '<e|e>=',ee2
    write(6,*) 'c**2=',c2
    write(6,*) 'DD=',DD
#   endif

    ! Reject if coefficients are too large (linear dep.).

    if ((sqrt(c2) > ThrCff) .and. (ee2 < 1.0D-5) .and. (.not. QNRStp)) then
#     ifdef _DEBUGPRINT_
      Fmt = '(A,i2,5x,g12.6)'
      Text = '|c| is too large,     iVec, |c| = '
      write(6,Fmt) Text(1:36),iVec,sqrt(c2)
#     endif
      cycle
    end if

    ! Keep the best candidate

    if (QNRStp) then
      if (ee2*DD < ee1*DD1) then
        !if (DD < DD1) then
        ee1 = ee2
        ipBst = iVec
        DD1 = DD
#       ifdef _DEBUGPRINT_
        cDotV = c2
#       endif
      end if
#     ifdef _NOT_USED_
      if ((ee2/ee1 > f1) .and. (ee2/ee1 < f2)) then
        if (DD < DD1) then
          ee1 = ee2
          ipBst = iVec
          DD1 = DD
#         ifdef _DEBUGPRINT_
          cDotV = c2
#         endif
        end if
      else if (ee2 < ee1) then
        ee1 = ee2
        ipBst = iVec
        DD1 = DD
#       ifdef _DEBUGPRINT_
        cDotV = c2
#       endif
      end if
#     endif
    else
      if (ee2 < ee1) then
        ! New vector lower eigenvalue.
        ee1 = ee2
        ipBst = iVec
#       ifdef _DEBUGPRINT_
        cDotV = c2
#       endif
      end if
    end if

  end do
  call mma_deallocate(Scratch)

  if ((ipBst < 1) .or. (ipBst > kOptim)) then
    write(6,*) ' No proper solution found in C2-DIIS !'
    Fmt = '(6es16.8)'
    Text = 'EigenValues :'
    call RecPrt(Text,Fmt,EValue,1,kOptim)
    Text = 'EigenVectors :'
    call RecPrt(Text,Fmt,EVector,kOptim,kOptim)
    call Quit_OnConvError()
  end if
  call dcopy_(kOptim,EVector(1,ipBst),1,CInter(1,1),1)

# ifdef _DEBUGPRINT_
  write(6,*)
  write(6,*) ' Selected root :',ipBst
  write(6,'(A,f16.8)') '  c**2 =         ',cDotV
# endif

  call mma_deallocate(EValue)
  call mma_deallocate(EVector)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !     C1DIIS                                                         *
  !                                                                    *
  !       References:                                                  *
  !       P. Csaszar and P. Pulay, J. Mol. Struc., 114, 31-34 (1984).  *
  !       doi:10.1016/S0022-2860(84)87198-7                            *
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (QNRStp) then
    AccCon = 'QNRc1DIIS'
  else
    AccCon = 'c1DIIS   '
  end if

  ! Set up the missing part of the matrix in Eq. (5) and the
  ! vector on the RHS in the same equation. Note the sign change!

  do i=1,kOptim
    Bij(kOptim+1,i) = -One ! note sign change
    Bij(i,kOptim+1) = -One ! note sign change
    GDiis(i) = Zero
  end do
  Bij(kOptim+1,kOptim+1) = Zero
  GDiis(kOptim+1) = -One  ! note sign change

# ifdef _DEBUGPRINT_
  write(6,*) ' B matrix in DIIS_e:'
  do i=1,kOptim+1
    write(6,'(7f16.8)') (Bij(i,j),j=1,kOptim+1),GDiis(i)
  end do
# endif

  ! Condition the B matrix

  B11 = sqrt(Bij(1,1)*Bij(kOptim,kOptim))
  do i=1,kOptim
    do j=1,kOptim
      Bij(i,j) = Bij(i,j)/B11
    end do
  end do

  ! Solve for the coefficients, solve the equations.

  call Gauss(kOptim+1,nBij,Bij,CInter(1,1),GDiis)

  ! Normalize sum of interpolation coefficients

  Fact = Zero
  do i=1,kOptim
    Fact = Fact+CInter(i,1)
  end do

  Fact = One/Fact
  do i=1,kOptim
    CInter(i,1) = Fact*CInter(i,1)
  end do

  ! Make sure new density gets a weight
  call C_Adjust(CInter(1,1),kOptim,0.05d0)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if

call mma_deallocate(Bij)

! Temporary fix for UHF.

if (nD == 2) call DCopy_(nCI,CInter(1,1),1,CInter(1,2),1)
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
Fmt = '(6es16.8)'
Text = 'The solution vector :'
call RecPrt(Text,Fmt,CInter(1,1),1,kOptim)
#endif

call Timing(Cpu2,Tim1,Tim2,Tim3)
TimFld(6) = TimFld(6)+(Cpu2-Cpu1)

return

end subroutine DIIS_x
