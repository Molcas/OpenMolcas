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

use Index_Functions, only: nTri_Elem
use InfSCF, only: AccCon, C1DIIS, Energy, Iter, Iter_Start, IterSO, kOptim, kOV, mOV, MxOptm, TimFld
use LnkLst, only: LLx
use Interfaces_SCF, only: OptClc_X
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Ten
#ifdef _NEW_
use Constants, only: Half
#endif
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nD, nCI
real(kind=wp), intent(inout) :: CInter(nCI,nD)
logical(kind=iwp), intent(in) :: QNRstp
integer(kind=iwp), intent(out) :: Ind(MxOptm)
integer(kind=iwp) :: i, iDiag, iDum, iErr, ij, ipBst, iVec, j, kVec, nBij, nFound
real(kind=wp) :: Alpha, B11, Bii_Min, c2, cpu1, cpu2, DD, DD1, Dummy, E_Min, E_Min_G, ee1, ee2, EMax, Fact, tim1, tim2, tim3
#ifdef _DEBUGPRINT_
real(kind=wp) :: cDotV
#endif
logical(kind=iwp) :: Case1, Case2, Case3
character(len=80) :: Frmt, Text
real(kind=wp), allocatable :: Bij(:,:), BijTri(:), Err1(:), Err2(:), EValue(:), EVector(:,:), GDiis(:), Scratch(:)
real(kind=wp), parameter :: CThr = 0.05_wp, delta = 1.0e-4_wp, delta_E = 1.0e-4_wp, Fact_Decline = 15.0_wp, ThrCff = Ten
real(kind=wp), external :: DDot_

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

call Timing(Cpu1,Tim1,Tim2,Tim3)
do
  Case1 = .false.
  Case2 = .false.
  Case3 = .false.

  ! Select from the kOptim last iterations

  do i=1,kOptim
    Ind(i) = iter-kOptim+i
  end do
# ifdef _DEBUGPRINT_
  !write(u6,*) 'Iter, Iter_Start=',Iter,Iter_Start
  write(u6,*) 'kOptim=',kOptim
  write(u6,*) 'Ind(i):',(Ind(i),i=1,kOptim)
# endif

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
  Bij(:,:) = Zero

  ! Compute norms, <e_i|e_j>

# ifdef _DEBUGPRINT_
  write(u6,*) 'kOV(:)=',kOV
  write(u6,*) 'mOV   =',mOV
  call RecPrt('Energy',' ',Energy,1,iter)
# endif
  E_Min_G = Zero
  E_Min = Zero
  Bii_min = 1.0e99_wp
  do i=1,kOptim
    call ErrV(mOV,Ind(i),QNRStp,Err1)
    !if (QNRStp) call ErrV(mOV,Ind(i),.false.,Err3)
#   ifdef _DEBUGPRINT_
    call NrmClc(Err1,mOV,'Diis  ','Err(i) ')
#   endif
    do j=1,i-1

      call ErrV(mOV,Ind(j),QNRStp,Err2)
      !if (QNRStp) call ErrV(mOV,Ind(j),.false.,Err4)
#     ifdef _DEBUGPRINT_
      call NrmClc(Err2,mOV,'Diis  ','Err(j)  ')
#     endif
      if (QNRStp) then
#       ifdef _NEW_
        !Bij(i,j) = Half*real(nD,kind=wp)*(DDot_(mOV,Err1,1,Err4,1)+DDot_(mOV,Err3,1,Err2,1))
#       else
        Bij(i,j) = real(nD,kind=wp)*DDot_(mOV,Err1,1,Err2,1)
#       endif
      else
        Bij(i,j) = real(nD,kind=wp)*DDot_(mOV,Err1,1,Err2,1)
      end if
      Bij(j,i) = Bij(i,j)
    end do
    if (QNRStp) then
#     ifdef _NEW_
      !Bij(i,i) = real(nD,kind=wp)*DDot_(mOV,Err1,1,Err3,1)
#     else
      Bij(i,i) = real(nD,kind=wp)*DDot_(mOV,Err1,1,Err1,1)
#     endif
    else
      Bij(i,i) = real(nD,kind=wp)*DDot_(mOV,Err1,1,Err1,1)
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
    if (Bij(i,i) < 1.0e-6_wp) cycle
    if (Fact_Decline*Bij(i,i) < Bij(i+1,i+1)) Case3 = .true.
  end do
  if (Energy(Ind(i)) >= E_min) Case3 = .false.

  if (qNRStp .and. (Case1 .or. Case2 .or. Case3)) then
#   ifdef _DEBUGPRINT_
    write(u6,*) 'Case1=',Case1
    write(u6,*) 'Case2=',Case2
    write(u6,*) 'Case3=',Case3
    write(u6,*) '   RESETTING kOptim!!!!'
    write(u6,*) '   Calculation of the norms in Diis :'
    Frmt = '(6(G0.12,2x))'
    Text = 'B-matrix squared in Diis :'
    call RecPrt(Text,Frmt,Bij,nBij,nBij)
    write(u6,'(A,2(G0.9,2x))') 'Bij(i,i),      Bii_Min=',Bij(i,i),Bii_Min
    write(u6,'(A,2F16.6)') 'Energy(Ind(i)),E_Min_G=',Energy(Ind(i)),E_Min_G
    write(u6,*) Energy(Ind(i)),E_Min_g
    write(u6,*)

#   endif
    ! Rest the depth of the DIIS and the BFGS update.
    if (Case3) then
      write(u6,*) 'DIIS_X: Resetting kOptim!'
      write(u6,*) '        Caused by inconsistent B matrix values.'
    else if (Case1) then
      write(u6,*) 'DIIS_X: Resetting BFGS depth!'
      write(u6,*) '        Too large B matrix values.'
    else if (Case2) then
      write(u6,*) 'DIIS_X: Resetting kOptim!'
      write(u6,*) '        Caused by energies and gradients which are inconsistent with a convex energy functional.'
    end if
    if (Case1) then
      ! The BFGS update is probably to blame. Reset the update depth.
      Iter_Start = Iter
      IterSO = 1
    else
      !if (Case2) then
      write(u6,*) 'kOptim=',kOptim,'-> kOptim=',1
      kOptim = 1
      Iter_Start = Iter
      IterSO = 1
      !else
      !  write(u6,*) 'kOptim=',kOptim,'-> kOptim=',kOptim-1
      !  kOptim = kOptim-1
      !  Iter_Start = Iter_Start+1
      !  IterSO = IterSO-1
      !end if
    end if
    call mma_deallocate(Err2)
    call mma_deallocate(Err1)
    call mma_deallocate(Bij)
  else
    do i=1,kOptim-1
      if (delta*sqrt(Bij(i,i)) > sqrt(Bij(kOptim,kOptim))) then
        write(u6,*) 'DIIS_X: Reduction of the subspace dimension due to numerical imbalance of the values in the B-Matrix'
        write(u6,*) 'kOptim=',kOptim,'-> kOptim=',kOptim-1
        kOptim = kOptim-1
        Iter_Start = Iter_Start+1
        IterSO = IterSO-1
        call mma_deallocate(Err2)
        call mma_deallocate(Err1)
        call mma_deallocate(Bij)
        exit
      end if
    end do
    if (i >= kOptim) exit
  end if

end do

! Deallocate memory for error vectors & gradient
!call mma_deallocate(Err4)
!call mma_deallocate(Err3)
call mma_deallocate(Err2)
call mma_deallocate(Err1)

#ifdef _DEBUGPRINT_
write(u6,*) '   Calculation of the norms in Diis :'
Frmt = '(6f16.8)'
Text = 'B-matrix squared in Diis :'
call RecPrt(Text,Frmt,Bij,nBij,nBij)
write(u6,*)
write(u6,*)
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

  call unitmat(EVector,kOptim)

  ! Form a triangular B-matrix

  call mma_allocate(BijTri,nTri_Elem(MxOptm),Label='BijTri')
  ij = 1
  do i=1,kOptim
    BijTri(ij:ij+i-1) = Bij(i,1:i)
    ij = ij+i
  end do

# ifdef _DEBUGPRINT_
  Frmt = '(5g25.15)'
  Text = 'B-matrix before Jacobi :'
  call TriPrt(Text,Frmt,BijTri,kOptim)
  Text = 'EigenVectors before Jacobi :'
  call RecPrt(Text,Frmt,EVector,kOptim,kOptim)
  write(u6,*)
  write(u6,*)
# endif

  ! Diagonalize B-matrix

  EMax = Zero
  do i=1,nTri_Elem(kOptim)
    EMax = max(EMax,abs(BijTri(i)))
  end do
  do i=1,nTri_Elem(kOptim)
    if (abs(BijTri(i)) < EMax*1.0e-14_wp) BijTri(i) = Zero
  end do

  call mma_allocate(Scratch,kOptim**2,Label='Scratch')

  Dummy = Zero
  iDum = 0
  call Diag_Driver('V','A','L',kOptim,BijTri,Scratch,kOptim,Dummy,Dummy,iDum,iDum,EValue,EVector,kOptim,1,0,'J',nFound,iErr)

  call mma_deallocate(Scratch)
  BijTri(1:nTri_Elem(kOptim)) = Zero

  iDiag = 0
  do i=1,kOptim
    iDiag = iDiag+i
    BijTri(iDiag) = EValue(i)
  end do

# ifdef _DEBUGPRINT_
  Frmt = '(5g25.15)'
  Text = 'B-matrix after Jacobi :'
  call TriPrt(Text,Frmt,BijTri,kOptim)
  Text = 'EigenValues :'
  call RecPrt(Text,Frmt,EValue,1,kOptim)
  Text = 'EigenVectors :'
  call RecPrt(Text,Frmt,EVector,kOptim,kOptim)
  write(u6,*)
  write(u6,*)
# endif

  ! Renormalize the eigenvectors to the C1-DIIS format

# ifdef _DEBUGPRINT_
  write(u6,*) ' Normalization constants :'
# endif

  do iVec=1,kOptim
    Alpha = Zero
    do i=1,kOptim
      Alpha = Alpha+EVector(i,iVec)
    end do
    Alpha = sum(EVector(:,iVec))

#   ifdef _DEBUGPRINT_
    Frmt = '(A7,i2,A4,f16.8)'
    write(u6,Frmt) ' Alpha(',iVec,') = ',Alpha
#   endif

    EVector(:,iVec) = EVector(:,iVec)/Alpha
  end do

  do kVec=1,kOptim
    ee1 = Zero
    do iVec=1,kOptim
      ee1 = ee1+EVector(iVec,kVec)*sum(EVector(:,kVec)*Bij(iVec,1:kOptim))
    end do
    !write(u6,*) 'EValue(kVec),ee1:',EValue(kVec),ee1
    EValue(kVec) = ee1
  end do

# ifdef _DEBUGPRINT_
  Frmt = '(6es16.8)'
  Text = 'B-matrix after scaling :'
  call TriPrt(Text,Frmt,BijTri,kOptim)
  Text = 'EigenValues after scaling :'
  call RecPrt(Text,Frmt,EValue,1,kOptim)
  Text = 'EigenVectors after scaling :'
  call RecPrt(Text,Frmt,EVector,kOptim,kOptim)
  write(u6,*)
  write(u6,*)
# endif
  call mma_deallocate(BijTri)

  ! Select a vector.
  ee1 = 1.0e72_wp
  DD1 = 1.0e72_wp
# ifdef _DEBUGPRINT_
  cDotV = 1.0e72_wp
# endif
  ipBst = -99999999
  call mma_allocate(Scratch,mOV,label='Scratch')
  do iVec=1,kOptim

    ! Pick up eigenvalue (ee2) and the norm of the eigenvector (c2).

    ee2 = EValue(iVec)
    c2 = DDot_(kOptim,EVector(1,iVec),1,EVector(1,iVec),1)
    if (QNRStp) then
      CInter(1:kOptim,1) = EVector(1:kOptim,iVec)
      if (nD == 2) CInter(:,2) = CInter(:,1)
      call OptClc_x(CInter,nCI,nD,Scratch,mOV,Ind,MxOptm,kOptim,kOV,LLx,DD)
    else
      DD = Zero
    end if
#   ifdef _DEBUGPRINT_
    write(u6,*) '<e|e>=',ee2
    write(u6,*) 'c**2=',c2
    write(u6,*) 'DD=',DD
#   endif

    ! Reject if coefficients are too large (linear dep.).

    if ((sqrt(c2) > ThrCff) .and. (ee2 < 1.0e-5_wp) .and. (.not. QNRStp)) then
#     ifdef _DEBUGPRINT_
      Frmt = '(A,i2,5x,g12.6)'
      Text = '|c| is too large,     iVec, |c| = '
      write(u6,Frmt) Text(1:36),iVec,sqrt(c2)
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
      if ((ee2/ee1 > Half) .and. (ee2/ee1 < Two)) then
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
    write(u6,*) ' No proper solution found in C2-DIIS !'
    Frmt = '(6es16.8)'
    Text = 'EigenValues :'
    call RecPrt(Text,Frmt,EValue,1,kOptim)
    Text = 'EigenVectors :'
    call RecPrt(Text,Frmt,EVector,kOptim,kOptim)
    call Quit_OnConvError()
  end if
  CInter(1:kOptim,1) = EVector(:,ipBst)

# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) ' Selected root :',ipBst
  write(u6,'(A,f16.8)') '  c**2 =         ',cDotV
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

  call mma_allocate(GDiis,MxOptm+1,Label='GDiis')
  Bij(kOptim+1,1:kOptim) = -One ! note sign change
  Bij(1:kOptim,kOptim+1) = -One ! note sign change
  Bij(kOptim+1,kOptim+1) = Zero
  GDiis(1:kOptim) = Zero
  GDiis(kOptim+1) = -One  ! note sign change

# ifdef _DEBUGPRINT_
  write(u6,*) ' B matrix in DIIS_e:'
  do i=1,kOptim+1
    write(u6,'(7f16.8)') (Bij(i,j),j=1,kOptim+1),GDiis(i)
  end do
# endif

  ! Condition the B matrix

  B11 = sqrt(Bij(1,1)*Bij(kOptim,kOptim))
  Bij(1:kOptim,1:kOptim) = Bij(1:kOptim,1:kOptim)/B11

  ! Solve for the coefficients, solve the equations.

  call Gauss(kOptim+1,nBij,Bij,CInter(:,1),GDiis)
  call mma_deallocate(GDiis)

  ! Normalize sum of interpolation coefficients

  Fact = One/sum(CInter(1:kOptim,1))
  CInter(1:kOptim,1) = Fact*CInter(1:kOptim,1)

  ! Make sure new density gets a weight
  if (CInter(kOptim,1) < CThr) then
    Fact = (One-CThr)/(One-CInter(kOptim,1))
    CInter(1:kOptim-1,1) = Fact*CInter(1:kOptim-1,1)
    CInter(kOptim,1) = CThr
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if

call mma_deallocate(Bij)

! Temporary fix for UHF.

if (nD == 2) CInter(:,2) = CInter(:,1)
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
Frmt = '(6es16.8)'
Text = 'The solution vector :'
call RecPrt(Text,Frmt,CInter(1,1),1,kOptim)
#endif

call Timing(Cpu2,Tim1,Tim2,Tim3)
TimFld(6) = TimFld(6)+(Cpu2-Cpu1)

return

end subroutine DIIS_x
