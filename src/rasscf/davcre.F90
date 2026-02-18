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

subroutine DAVCRE(C,HC,HH,CC,E,HD,SC,Q,QQ,S,SXSEL,NROOT,ITMAX,NDIM,ITERSX,NSXS)
! RASSCF program: version IBM-3090: SX section
!
! Using the Davidson method, in the multiple root version suggested
! by B.LIU. This routine finds the NROOT lowest eigenvalues and
! eigenvectors of a secular problem of dimension NDIM.
! This is a vectorized version adapted to run optimally on IBM 3090
! It is completely CPU bound (except maybe for the construction of
! the sigma vector). Multiple calls to this subroutine may therefore
! be necessary for cases where more iterations are needed than
! is allowed by the core space requirements to store all C and
! sigma vectors.
!
! Externals: HMAT (set up the Davidson H-matrix HH)
!            COVLP (calculate the overlap between two CI vectors)
!            Jacob (full diagonalization routine)
!            DGEMM and other ESSL routines for matrix operations
! Parameters: C  SuperCI-vectors
!             HC SuperCI Sigma vectors
!             HH Davidson's H-matrix
!             CC      "     eigenvectors
!             E       "     eigenvalues
!             HD diagonal elements of SuperCI
!             Q  the Davidson update vectors
!             QQ the norm of all Q-vectors
!             SC scratch area
!
! ********** IBM-3090 Release 88 09 08 *****

use Index_Functions, only: nTri_Elem
use fciqmc, only: DoNECI
use wadr, only: DIA, PA, SXN
use PrintLevel, only: DEBUG, INSANE
use output_ras, only: IPRLOC, RC_SX
use RASDim, only: MxSXIt
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: NROOT, ITMAX, NDIM, ITERSX, NSXS
real(kind=wp) :: C((NROOT+NSXS)*NROOT*(ITMAX+1)), HC((NROOT+NSXS)*NROOT*ITMAX), HH((ITMAX*NROOT)*(ITMAX*NROOT+1)), &
                 CC((ITMAX*NROOT)**2), E((ITMAX*NROOT)), HD(NROOT+NSXS), SC((NROOT+NSXS)), Q((NROOT+NSXS)*(NROOT+1)), QQ(NROOT), &
                 S(ITMAX*NROOT**2)
character(len=*) :: SXSEL
integer(kind=iwp) :: i, iConvA, iConvL, iConvQ, ii, ij, iPass, iPrLev, iSel, iST, iSTC, iSTQ, j, ji, jST, k, kDimH, Length, nCR, &
                     nDimH, nDimH2, nST, nTotDC, nTrial
real(kind=wp) :: ASQ, Ei, ENO, Ovl, QNorm, SWAP, XMX, XNorm, XX
character(len=4) :: IOUTW, IOUTX
real(kind=wp), allocatable :: C1(:), C2(:), X(:)
real(kind=wp), parameter :: THRA = 1.0e-13_wp, THRLD1 = 1.0e-8_wp, THRLD2 = 5.0e-14_wp, THRQ = 1.0e-7_wp, THRZ = 1.0e-6_wp
real(kind=wp), external :: DDot_
#include "warnings.h"

! Local print level (if any)
IPRLEV = IPRLOC(1)
if (IPRLEV >= DEBUG) then
  write(u6,*) ' Entering DAVCRE'
  write(u6,*) 'Super-CI diagonalization. Max iterations: ',ITMAX
end if

! Set values at loop head and starting guess of the vectors C

Rc_SX = 0
NTRIAL = NROOT
NCR = NROOT*NDIM
C(1:NCR) = Zero
II = 0
do I=1,NROOT
  C(II+I) = One
  II = II+NDIM
end do

! memory allocation for COVLP

call mma_allocate(C1,NSXS,Label='C1')
call mma_allocate(C2,NSXS,Label='C2')
call mma_allocate(X,NSXS,Label='X')

! Begin Davidson iterations
! Start by setting up the Davidson HH-matrix
! The dimension of this matrix, NDIMH, is updated in HMAT

NDIMH = 0

!*************************************************************
!               Head of Super-CI iteration loop:
!*************************************************************
do ITERSX=1,ITMAX
  call HMAT(C,HC,HH,HD,NDIM,NDIMH,NTRIAL)
  KDIMH = nTri_Elem(NDIMH)

  ! Echo the HH-matrix before diagonalizing it.

  HH(KDIMH+1:2*KDIMH) = HH(1:KDIMH)

  if (IPRLEV >= DEBUG) then
    write(u6,*) ' Davidson H-matrix in iteration ',ITERSX
    write(u6,*) ' Davidson H-matrix triangular of size =  ',NDIMH
    write(u6,'(1x,8F14.6)') (HH(I),I=1,KDIMH)
  end if
  E(1) = HH(1)
  CC(1) = One
  if (NDIMH > 1) then

    ! set eigenvector array to identity before JACO call

    NDIMH2 = NDIMH**2
    CC(1:NDIMH2) = Zero
    II = -NDIMH
    do I=1,NDIMH
      II = II+NDIMH+1
      CC(II) = One
    end do

    call Jacob(HH(KDIMH+1),CC,NDIMH,NDIMH)
    if (SXSEL == 'LOWEST') then
      ! Root selection here assumes picking the lowest root(s).
      call JACORD(HH(KDIMH+1),CC,NDIMH,NDIMH)
      II = 0
      do I=1,NROOT
        II = II+I
        E(I) = HH(II+KDIMH)
      end do
    else
      ! Root selection by maximum overlap.
      ISEL = 1
      XMX = abs(CC(1))
      do I=2,NDIMH
        XX = abs(CC(1+NDIMH*(I-1)))
        if (XX > XMX) then
          ISEL = I
          XMX = XX
        end if
      end do
      if (ISEL /= I) then
        do I=1,NDIMH
          SWAP = CC(I+NDIMH*(ISEL-1))
          CC(I+NDIMH*(ISEL-1)) = CC(I)
          CC(I) = -SWAP
        end do
        SWAP = HH(KDIMH+nTri_Elem(ISEL))
        HH(KDIMH+nTri_Elem(ISEL)) = HH(KDIMH+1)
        HH(KDIMH+1) = SWAP
      end if
    end if

    if (IPRLEV >= DEBUG) then
      write(u6,*) ' Eigenvalues:'
      write(u6,'(1X,8F14.6)') (E(I),I=1,NROOT)
      write(u6,*) ' Eigenvectors:'
      NST = 0
      do I=1,NROOT
        write(u6,'(1X,10F11.6)') (CC(J+NST),J=1,NDIMH)
        NST = NST+NDIMH
      end do
    end if
  end if

  ! Now perform the Davidson update
  ! First form for each root I, a vector Q(I), where
  ! Q(I,K)=sum(J,M) CC(I,J,M)*(HC(M,J,K)-E(I)*C(M,J,K))
  ! here J runs over iterations and M over the number of trial
  ! vectors in each iterations E(I) is the current estimate of the
  ! energy of the I:th root, C(M,J) is the M:th trial vector used in
  ! the J:th iteration, HC being the corresponding sigma vector. CC
  ! is the matrix of eigenvectors of the Davidson Hamiltonian.

  ! STEP 1: form the matrix CCE=-CC*E (CCE in array SC)
  JI = 0
  do I=1,NROOT
    EI = E(I)
    do J=1,NDIMH
      JI = JI+1
      S(JI) = -CC(JI)*EI
    end do
  end do

  ! Step 2: form the matrix Q=HC*CC
  call DGEMM_('N','N',NDIM,NROOT,NDIMH,One,HC,NDIM,CC,NDIMH,Zero,Q(NDIM+1),NDIM)
  ! Step 3: add the contribution C*CCE
  call DGEMM_('N','N',NDIM,NROOT,NDIMH,One,C,NDIM,S,NDIMH,One,Q(NDIM+1),NDIM)

  ! Note that the NDIM first positions in Q are untouched so far
  ! the vector will be moved later
  if (IPRLEV >= INSANE) then
    write(u6,*) ' The residual vectors Q of size:',NDIM
    write(u6,'(1x,8F14.10)') (Q(NDIM+I),I=1,NDIM)
  end if

  ! Check the norms of the Q-vectors for convergence, and obtain
  ! bounds to the eigenvalues. E(I) is an upper bound to the I:th
  ! eigenvalue, from McDonald's theorem, while E(I)-NORM(Q(I))
  ! is a lower bound to the eigenvalue. This is the Weinstein lower
  ! bound formula in a multiple root form.

  ICONVQ = 0
  IST = 1+NDIM
  do I=1,NROOT
    call COVLP(Q(IST),Q(IST),DIA,PA,SXN,C1,C2,X,QQ(I))
    IST = IST+NDIM
  end do

  do I=1,NROOT
    QNORM = sqrt(QQ(I))
    if (QNORM < THRQ) ICONVQ = ICONVQ+1
    EI = E(I)
    ENO = EI-QNORM
    if (IPRLEV >= DEBUG) then
      if (NROOT > 1) then
        if (I > 1) then
          write(u6,'(20X,F16.8,A,I2,A,F16.8)') ENO,' <  SX energy ',I,'  < ',EI
        else
          write(u6,'(1X,I2,A,4X,F16.8,A,I2,A,F16.8)') ITERSX,' SX iteration ',ENO,' <  SX energy ',I,'  < ',EI
        end if
      else
        write(u6,'(1X,I2,A,4X,F16.8,A,I2,A,F16.8)') ITERSX,' SX iteration ',ENO,' <  SX energy ',I,'  < ',EI
      end if
      write(u6,*) '  Norm of Q:',QNORM
    end if
    IST = IST+NDIM
  end do
  !PAM00 End of replacement.

  ! Reset ICONVQ unless all roots are converged

  if (ICONVQ < NROOT) ICONVQ = 0

  ! Check the expansion vectors for convergence. This is done
  ! by computing the sum of the squares of CC(I,J), for each
  ! root I, and for the J trial vectors of the current iteration

  NST = NDIMH-NTRIAL+1
  ICONVA = 0
  do I=1,NROOT
    ASQ = DDOT_(NTRIAL,CC(NST),1,CC(NST),1)
    if (ASQ < THRA) ICONVA = ICONVA+1
    NST = NST+NDIMH
  end do

  ! Reset ICONVA unless all roots are converged.

  if (ICONVA < NROOT) ICONVA = 0

  ! Branch out of iteration loop if convergence has been
  ! achieved on either the energy or the Davidson expansion vectors.

  if ((ICONVQ /= 0) .or. (ICONVA /= 0)) then
    IOUTW = 'y   '
    if (NROOT > 1) IOUTW = 'ies '
    IOUTX = '    '
    if (NROOT > 1) IOUTX = 's   '
    if (IPRLEV >= DEBUG) then
      if (ICONVQ /= 0) write(u6,*) ' Convergence on CI energ'//IOUTW
      if (ICONVA /= 0) write(u6,*) ' Convergence on Davidson expansion vector'//IOUTX
    end if
    exit
  end if

  ! Calculate the D-vectors as D(I,K)=Q(I,K)/(E(I)-HD(K,K))

  IST = 1
  do I=1,NROOT
    EI = E(I)
    do K=1,NDIM
      SC(K) = EI-HD(K)
      if (abs(SC(K)) < THRZ) SC(K) = One
    end do
    Q(IST:IST+NDIM-1) = Q(IST+NDIM:IST+2*NDIM-1)/SC(1:NDIM)
    IST = IST+NDIM
  end do
  ! Remove any unwanted components. These are signalled by
  ! huge elements of SX hamiltonian diagonal (set in SXHAM).
  do I=1,NDIM
    if (HD(I) > 1.0e20_wp) then
      do J=1,NROOT
        Q(I+NDIM*(J-1)) = Zero
      end do
    end if
  end do
  if (IPRLEV >= INSANE) then
    write(u6,*) ' The correction vectors D(stored in Q):'
    write(u6,'(1x,8F14.10)') (Q(I),I=1,NDIM)
  end if

  ! Q now contains the Davidson correction vectors. These will now
  ! be orthogonalized to all old vectors and then to each other.
  !PAM01 The ON is done in two passes. A vector smaller than THRLD
  ! after the first ON pass is considered to contribute nothing worthwhile
  ! to the Davidson procedure, and is discarded. Else, it is normalized.
  ! If it is not still normalized to within an accuracy THRLD after the
  ! second pass, then again it is discarded: it means that the rescaling
  ! after pass 1, together with numerical noise amplification in
  ! COVLP, is beginning to be noticable. The latter will of course happen
  ! when some orbital rotation(s) are almost redundant.

  IPASS = 0
  ICONVL = 0
  NTOTDC = NROOT

  do

    ! First form the overlap matrix
    IJ = 0
    ISTQ = 1
    do I=1,NTOTDC
      ISTC = 1
      do J=1,NDIMH
        IJ = IJ+1
        call COVLP(Q(ISTQ),C(ISTC),DIA,PA,SXN,C1,C2,X,OVL)
        S(IJ) = -OVL
        ISTC = ISTC+NDIM
      end do
      ISTQ = ISTQ+NDIM
    end do
    call DGEMM_('N','N',NDIM,NTOTDC,NDIMH,One,C,NDIM,S,NDIMH,One,Q,NDIM)

    ! The Q-vectors are now orthogonal to all C-vectors.
    if (IPRLEV >= INSANE) then
      write(u6,*) '  D vector orthogonal to all C vectors:'
      write(u6,'(1x,8F14.10)') (Q(I),I=1,NDIM)
    end if
    ! Now orthogonalize them to one another, rejecting those which
    ! appear with too small a norm

    IST = 1
    NTRIAL = 0

    ! Long loop over NTOTDC vectors:
    do I=1,NTOTDC
      if ((I /= 1) .and. (NTRIAL /= 0)) then

        ! Orthogonalize this vector to the preceding Q's of this iteration

        JST = 1
        do J=1,NTRIAL
          call COVLP(Q(IST),Q(JST),DIA,PA,SXN,C1,C2,X,OVL)
          S(J) = -OVL
          JST = JST+NDIM
        end do
        !call DGEMX(NDIM,NTRIAL,One,Q,NDIM,S,1,Q(IST),1)
        call DGEMV_('N',NDIM,NTRIAL,One,Q,NDIM,S,1,One,Q(IST),1)
        if (IPRLEV >= INSANE) then
          write(u6,*) '  Q vector orthogonal to preceding Q vectors:'
          write(u6,'(1x,8F14.10)') (Q(k),k=1,NDIM)
        end if
      end if

      ! Normalize this vector and move to trial set if norm large enough

      call COVLP(Q(IST),Q(IST),DIA,PA,SXN,C1,C2,X,XNORM)
      ! Due to large noise amplification in COVLP, the squared-norm
      ! can actually come out as a negative number.
      ! Acceptable, only if it is very close to zero. Else, quit.
      if (XNORM < -1.0e-9_wp) then
        write(u6,*)
        write(u6,*) '      *** Error in subroutine DAVCRE ***'
        write(u6,*) ' The squared norm of a trial vector has been'
        write(u6,*) ' computed to be negative:'
        write(u6,*) '      XNORM=',XNORM
        write(u6,*) ' This is possible only for some severe malfunction'
        write(u6,*) ' of the rasscf program. Please issue a bug report.'
        write(u6,*)
        if (.not. DoNECI) then
          call Quit(_RC_GENERAL_ERROR_)
        else
          write(u6,*) ' non positive-semi definite matrix occurred.'
          write(u6,*) ' Calculation will continue. '
          write(u6,*) ' Divergent results might occur'
          write(u6,*) ' Tests for possible solution on the way...'
        end if
      end if
      XNORM = sqrt(max(Zero,XNORM))
      if (IPRLEV >= INSANE) write(u6,'(1X,A,I3,A,I3,A,ES16.8)') 'Pass ',IPASS,' New orthogonal vector ',I,' has norm ',XNORM

      !PAM01 Two different treatments, depending on if this is first or
      ! second orthonormalization pass:
      if (IPASS == 0) then
        ! First pass:
        if ((ITERSX == 1) .or. (XNORM > THRLD1)) then
          ISTQ = NTRIAL*NDIM+1
          NTRIAL = NTRIAL+1
          XNORM = One/(XNORM+1.0e-24_wp)
          !PAM01 Note that ISTQ can be (and is!) the same as IST:
          Q(ISTQ:ISTQ+NDIM-1) = XNORM*Q(IST:IST+NDIM-1)
        end if
      else
        ! Second pass: Demand accurate normalization, else we know that
        ! poor independence, maybe with strong rounding-error amplification
        ! in COVLP, has began to erode the orthonormalization of the
        ! basis vectors, hence the integrity of the Davidson procedure.
        if ((ITERSX == 1) .or. (abs(XNORM-One) < THRLD2)) then
          ISTQ = NTRIAL*NDIM+1
          NTRIAL = NTRIAL+1
          XNORM = One/(XNORM+1.0e-24_wp)
          !PAM01 Note that ISTQ can be (and is!) the same as IST:
          Q(ISTQ:ISTQ+NDIM-1) = XNORM*Q(IST:IST+NDIM-1)
        end if
      end if
      IST = IST+NDIM

    ! End over long loop over NTOTDC vectors I=1..NTOTDC
    end do

    ! NTRIAL new orthogonal vectors have now been formed. if NTRIAL
    ! equals zero there are no new linearly independent trial vectors
    ! and we branch out. If NTRIAL does not equal zero make a second
    ! orthonormalization pass

    if (NTRIAL == 0) ICONVL = 1
    NTOTDC = NTRIAL
    IPASS = IPASS+1
    if ((IPASS == 2) .or. (NTOTDC == 0)) exit
  end do
  if (IPRLEV >= INSANE) then
    write(u6,*) ' The correction vectors, after ON:'
    write(u6,'(1x,8F14.10)') (Q(I),I=1,NDIM*NTRIAL)
  end if

  ! Check if the set of vectors has become linearly dependent

  if (ICONVL == 1) then
    if (IPRLEV >= DEBUG) write(u6,*) ' Trial vector set has become linearly dependent'
    exit
  end if

  ! Move the new orthogonal vectors to C

  NST = 1+NDIMH*NDIM
  LENGTH = NTRIAL*NDIM
  C(NST:NST+LENGTH-1) = Q(1:LENGTH)

  if (IPRLEV >= DEBUG) write(u6,'(1X,A,I2,A)') ' Adding ',NTRIAL,' new vectors.'
  if (IPRLEV >= INSANE) then
    IST = NDIMH*NDIM
    do I=1,NTRIAL
      write(u6,'(1X,A,I2,A)') ' New vector ',I,' is:'
      write(u6,'(1X,10F11.6)') (C(IST+K),K=1,NDIM)
      IST = IST+NDIM
    end do
  end if
end do

if (ITERSX > ITMAX) then
  ! At this point, the calculation has neither converged, nor
  ! produced linearly dependent trial vectors.

  if (IPRLEV >= DEBUG) then
    if (ITMAX < MXSXIT) then
      write(u6,*) ' Super-CI not converged. Max SX iter increased.'
    else
      write(u6,*) ' Super-CI not converged.'
    end if
  end if
  ITMAX = min(ITMAX+2,MXSXIT)
  Rc_SX = 16
end if

! Compute CI vectors for all roots
!
! CI NEW (K) = sum(J,M)CC(I,J M)*C M J (K)
!
! Here J runs over the iterations, and M over the
! number of trial vectors used in each iteration J

call DGEMM_('N','N',NDIM,NROOT,NDIMH,One,C,NDIM,CC,NDIMH,Zero,Q,NDIM)
if (IPRLEV >= INSANE) then
  write(u6,*) ' Unnormalized final CI vectors (in Q):'
  write(u6,'(1x,8F14.10)') (Q(I),I=1,NDIM)
end if

! Normalize the final CI vectors

IST = 1
do I=1,NROOT
  call COVLP(Q(IST),Q(IST),DIA,PA,SXN,C1,C2,X,XNORM)
  XNORM = One/XNORM
  C(IST:IST+NDIM-1) = XNORM*Q(IST:IST+NDIM-1)
  IST = IST+NDIM
end do

! Print vector if desired

if (IPRLEV >= INSANE) then
  IST = 0
  do I=1,NROOT
    write(u6,*) ' SX-CI vector for root ',I
    write(u6,'(1X,8F8.4)') (C(IST+K),K=1,NDIM)
    IST = IST+NDIM
  end do
end if

! End of diagonalization
! Free memory for COVLP

call mma_deallocate(C1)
call mma_deallocate(C2)
call mma_deallocate(X)

end subroutine DAVCRE
