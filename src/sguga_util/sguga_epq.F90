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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!               2026, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
! E_pq|Psi> in the split graph (SG_Epq_Psi) and the operators built on it
! (E_pq E_qp diagonal, CI transformation SG_TRACI).
! (Default public, as SGUGA_Core; the interface is the module SGUGA.)
module SGUGA_Epq

use SGUGA_Core

#include "intent.fh"

implicit none

contains

subroutine SG_Epq_Psi(SGS,CIS,EXS,IP,IQ,CPQ,ISYCI,CI,SGM,CI_Tag)

  use Symmetry_Info, only: Mul
  use Constants, only: Zero, One
  use Definitions, only: wp, iwp

#include "intent.fh"

  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(inout), target :: EXS
  integer(kind=iwp), intent(in) :: IP, IQ, ISYCI
  real(kind=wp), intent(in) :: CPQ, CI(*)
  real(kind=wp), intent(_OUT_) :: SGM(*)
! CI_Tag: the tag of the vector CI (SG_New_Tag); the intermediates of a
! previous call with the same tag may be reused. SG_No_Reuse: no reuse.
! Without it the vector is recognized by its fingerprints (CI_Fingerprint,
! a pass over the vector).
  integer(kind=iwp), intent(in), optional :: CI_Tag
  integer(kind=iwp) :: I, IC, ICS, INDEO, IOC, IOLW, IOUW, IPPOW, IPSHFT, ISGSTA, ISTA, ISYDC, ISYDSG, ISYP, ISYPQ, ISYQ, ISYSGM, &
                       ISYUC, ISYUSG, J, JC, JSTA, LICP, LLW, LUW, MVSGM, NCP, NDWNC, NDWNSG, NS1, NUPC, &
                       NUPSG, NCP1, NCP2, MV, MVX, nCSFs, LICP1, LICP2, JQ, JP
  real(kind=wp) :: X

! declarations to facilitate the reuse of sigma vectors if possible (the key
! is kept in EXS, next to SGTMP)
  logical(kind=iwp) :: Reuse_Sigma, Swap, Test1, Test2, Moved, Tagged, Same_Key
  integer(kind=iwp) :: iTag
  real(kind=wp) :: CI_ID(2)
  integer(kind=iwp) :: iOff, jOff

!*********************************************************************
!  GIVEN ACTIVE LEVEL INDICES IP AND IQ, AND INPUT CI ARRAYS
!  CI AND SGM, THIS ROUTINE ADDS TO SGM THE RESULT OF ACTING ON
!  CI WITH THE NUMBER CPQ TIMES THE EXCITATION OPERATOR E(IP,IQ).
!  THE ADDITIONAL ENTRIES IN THE PARAMETER LIST ARE TABLES THAT
!  WERE PREPARED BY GINIT AND ITS SUBROUTINES.
!*********************************************************************

  nCSFs = CIS%NCSF(ISYCI)
! SYMMETRY OF ORBITALS:
  ISYP = SGS%ISm(IP)
  ISYQ = SGS%ISm(IQ)
  ISYPQ = Mul(ISYP,ISYQ)
! SYMMETRY OF SIGMA ARRAY:
  ISYSGM = Mul(ISYPQ,ISYCI)

  if (IQ /= IP) then

! Unified code for the case P/=Q
    Swap = IQ > IP
    Test1 = merge(IQ <= SGS%MidLev,IP <= SGS%MidLev,Swap)
    Test2 = merge(SGS%MidLev < IP,SGS%MidLev < IQ,Swap)

    if (Test1) then

      ! EXCITING CASE, IQ<IP<=MIDLEV.
      INDEO = merge(2*SGS%nLev+(IQ*(IQ-1))/2+IP,2*SGS%nLev+(IP*(IP-1))/2+IQ,Swap)
      do MVSGM=1,CIS%nMidV
        do ISYUSG=1,SGS%nSym

          NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
          if (NS1 == 0) cycle
          ISYDSG = Mul(ISYUSG,ISYSGM)           ! compute the lower symmetry
          ISYDC = Mul(ISYPQ,ISYDSG)             ! compute the symmetry of the sigma vector
          NDWNC = CIS%NOW(2,ISYDC,MVSGM)
          if (NDWNC == 0) cycle
          ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM) ! get the off-set to the sigma vector block
          NUPSG = CIS%NOW(1,ISYUSG,MVSGM)           ! number of upper half-walks by symmetry and midvertex.
          IOC = CIS%IOCSF(ISYUSG,MVSGM,ISYCI)       ! get the off-set to the CI vector block
          NCP = merge(EXS%NOCP(INDEO,ISYDSG,MVSGM),EXS%NOCP(INDEO,ISYDC,MVSGM),Swap)
          if (NCP == 0) cycle
          NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)

          ! get the off-set to the block of compressed coupling coefficients.
          LICP = merge(EXS%IOCP(INDEO,ISYDSG,MVSGM),EXS%IOCP(INDEO,ISYDC,MVSGM),Swap)
          ! CASE IS: LOWER HALF, EXCITE:
          call sort_icoup_block(EXS%ICOUP(1,LICP+1),NCP,Swap,Moved)
          ! the walk index (EXS%IOrd) of a reordered block must follow
          if (Moved .and. allocated(EXS%IOrd)) call Index_Coup_Block(EXS,LICP,NCP)
          call Apply_col(CPQ,NUPSG,NDWNC,CI(IOC+1),NDWNSG,SGM(ISGSTA+1),NCP,EXS%ICOUP(1,LICP+1),swap=Swap)

        end do
      end do
    else if (Test2) then
      ! EXCITING CASE, MIDLEV<IQ<IP
      INDEO = merge(2*SGS%nLev+(IQ*(IQ-1))/2+IP,2*SGS%nLev+(IP*(IP-1))/2+IQ,Swap)
      do MVSGM=1,CIS%nMidV
        do ISYUSG=1,SGS%nSym
          NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
          if (NS1 == 0) cycle
          ISYUC = Mul(ISYPQ,ISYUSG)
          NUPC = CIS%NOW(1,ISYUC,MVSGM)
          if (NUPC == 0) cycle
          ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
          ISYDSG = Mul(ISYUSG,ISYSGM)
          NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
          IOC = CIS%IOCSF(ISYUC,MVSGM,ISYCI)
          NCP = merge(EXS%NOCP(INDEO,ISYUSG,MVSGM),EXS%NOCP(INDEO,ISYUC,MVSGM),Swap)
          if (NCP == 0) cycle
          LICP = merge(EXS%IOCP(INDEO,ISYUSG,MVSGM),EXS%IOCP(INDEO,ISYUC,MVSGM),Swap)
          ! CASE IS: UPPER HALF, EXCITE:
          call apply_row(CPQ,NDWNSG,NUPC,CI(IOC+1),NUPSG,SGM(ISGSTA+1),NCP,EXS%ICOUP(1,LICP+1),swap=Swap)
        end do
      end do
    else
      ! EXCITING CASE, IQ<=MIDLEV<IP
      !
      !   The canonical case, p>q
      !   sigma(Jl(m),Ju(m)) = Sum(Il(m'),Iu(m')) C(Il(m'),Iu(m'))
      !                       x A^p(Ju(m),Iu(m')) x A^q(Ju(m),Iu(m'))
      iOff = 1
      JP = merge(IQ,IP,SWAP)
      JQ = merge(IP,IQ,SWAP)
      ! The intermediates depend on the CI vector, on the upper-half orbital JP
      ! and on the symmetry of the lower-half orbital JQ.
      ! A tagged vector is recognized by its tag (the fingerprints only as the
      ! sanity check of SG_Checks), an untagged one by its fingerprints.
      Reuse_Sigma = .false.
      CI_ID(:) = Zero
      iTag = 0
      if (present(CI_Tag)) iTag = CI_Tag
      Tagged = iTag > 0
      if (EXS%Reuse_SGTMP .and. (iTag >= 0)) then
        Same_Key = EXS%Save_Valid .and. (EXS%Save_Owner == 1) .and. (EXS%Save_Swap .eqv. Swap) .and. (EXS%Save_U == JP) .and. &
                   (EXS%Save_LSym == SGS%ISm(JQ)) .and. (EXS%Save_CISym == ISYCI) .and. (EXS%Save_nCSF == nCSFs) .and. &
                   (EXS%Save_Tag == iTag)
        if (Tagged) then
          Reuse_Sigma = Same_Key
          if (SG_Checks) then
            CI_ID(:) = CI_Fingerprint(CI,nCSFs)
            if (Reuse_Sigma .and. any(EXS%Save_ID(:) /= CI_ID(:))) then
              write(u6,*) 'SG_Epq_Psi: the same tag for different CI vectors, tag ',iTag
              call Abend()
            end if
          end if
        else
          CI_ID(:) = CI_Fingerprint(CI,nCSFs)
          Reuse_Sigma = Same_Key .and. all(EXS%Save_ID(:) == CI_ID(:))
        end if
      end if
      ! Loop over all mid vertex
      do MVSGM=1,CIS%nMidV  ! m
        ! loop over delta b = -1 and +1
        do MV=1,2
          ! Find the connecting vertex, MVX=m' , in the case it exist.
          !
          ! In the the canonical case, p>q, we have the matrix element
          !
          !   A^p(Ju(m),Iu(m')), m'>m
          !
          !   Accessed INDEO,Gamma(Iu),m'
          MVX = merge(EXS%MVR(MVSGM,MV),EXS%MVL(MVSGM,MV),SWAP)
          if (MVX == 0) cycle
          if ((Swap .and. (MVX < MVSGM)) .or. ((.not. Swap) .and. (MVX > MVSGM))) then
            write(u6,*) 'SG_Epq_Psi: midvertex partner in the wrong order (MV, MVSGM, MVX, Swap):',MV,MVSGM,MVX,Swap
            call Abend()
          end if

          do ISYUSG=1,SGS%nSym
            NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
            if (NS1 == 0) cycle
            ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
            NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
            ISYDSG = Mul(ISYUSG,ISYSGM)
            ISYUC = merge(Mul(ISYQ,ISYUSG),Mul(ISYP,ISYUSG),Swap)
            ISYDC = merge(Mul(ISYP,ISYDSG),Mul(ISYQ,ISYDSG),Swap)

            NUPC = CIS%NOW(1,ISYUC,MVX)
            if (NUPC == 0) cycle
            NDWNC = CIS%NOW(2,ISYDC,MVX)
            if (NDWNC == 0) cycle

!       INDEO =merge(JP, SGS%nLev+JP, MV==1)
            INDEO = (MV-1)*SGS%nLev+JP
            NCP1 = merge(EXS%NOCP(INDEO,ISYUSG,MVSGM),EXS%NOCP(INDEO,ISYUC,MVX),Swap)
            if (NCP1 == 0) cycle

            ! CASE IS: UPPER HALF, EXCITE:
            LICP1 = merge(EXS%IOCP(INDEO,ISYUSG,MVSGM),EXS%IOCP(INDEO,ISYUC,MVX),Swap)
            IOC = CIS%IOCSF(ISYUC,MVX,ISYCI)

            ! IN CASE OF REUSE COMPUTE THE TEMPORARY SIGMA VECTOR REGARDLESS OF THE NCP2 VALUE.
            ! Put it at the offset jOff
            jOff = iOff
            if (EXS%Reuse_SGTMP .and. (.not. Reuse_Sigma)) then
              call Build_Crossing_Temporary( &
                NDWNC,NUPC,CI(IOC+1),NUPSG,EXS%SGTMP(jOff), &
                NCP1,EXS%ICOUP(1,LICP1+1),Swap)
            end if

            ! Update the pointer in SGTMP to the next slot.
            if (EXS%Reuse_SGTMP) then
              iOff = iOff+NUPSG*NDWNC
              if (iOff-1 > size(EXS%SGTMP)) then
                write(u6,*) 'iOff-1>SIZE(EXS%SGTMP) -- 1'
                write(u6,*) iOff-1,size(EXS%SGTMP)
                call Abend()
              end if
            end if

!       INDEO = merge(JQ, SGS%nLev+JQ, MV==1)
            INDEO = (MV-1)*SGS%nLev+JQ
            NCP2 = merge(EXS%NOCP(INDEO,ISYDSG,MVSGM),EXS%NOCP(INDEO,ISYDC,MVX),Swap)
            if (NCP2 == 0) cycle

            ! In case of no reuse compute the sigma vector here.
            if (.not. EXS%Reuse_SGTMP) then
              jOff = 1
              call Build_Crossing_Temporary( &
                NDWNC,NUPC,CI(IOC+1),NUPSG,EXS%SGTMP(jOff), &
                NCP1,EXS%ICOUP(1,LICP1+1),Swap)
            end if

            ! CASE IS: LOWER HALF, EXCITE:
            NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
            LICP2 = merge(EXS%IOCP(INDEO,ISYDSG,MVSGM),EXS%IOCP(INDEO,ISYDC,MVX),Swap)
            call Contract_Crossing_Temporary( &
              CPQ,NUPSG,NDWNC,EXS%SGTMP(jOff),NDWNSG, &
              SGM(ISGSTA+1),NCP2,EXS%ICOUP(1,LICP2+1),Swap)

          end do
        end do
      end do

      if (EXS%Reuse_SGTMP) then
        if (.not. Reuse_Sigma) then
          EXS%Save_Valid = .true.
          EXS%Save_Owner = 1
          EXS%Save_Swap = Swap
          EXS%Save_U = JP
          EXS%Save_LSym = SGS%ISm(JQ)
          EXS%Save_CISym = ISYCI
          EXS%Save_nCSF = nCSFs
          EXS%Save_Tag = iTag
          EXS%Save_ID(:) = CI_ID(:)
          ! a vector used only once leaves nothing to reuse
          if (iTag < 0) EXS%Save_Valid = .false.
        end if
      else
        ! SGTMP was overwritten
        EXS%Save_Valid = .false.
      end if

    end if

  else

    ! THEN THIS IS A SPECIAL CASE.
    if (IP > SGS%MidLev) then

      ! SPECIAL CASE: WEIGHT OPERATOR, IP=IQ.
      ! IP=IQ>MIDLEV
      do MVSGM=1,CIS%nMidV
        do ISYUSG=1,SGS%nSym
          NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
          if (NS1 == 0) cycle
          ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
          ISYDSG = Mul(ISYUSG,ISYSGM)
          NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
          IOUW = CIS%IOW(1,ISYUSG,MVSGM)
          IPSHFT = 2*(IP-1-SGS%MidLev)
          LUW = 1+IOUW-CIS%nIpWlk+IPSHFT/(2*CSFG_nPack)
          IPSHFT = mod(IPSHFT,2*CSFG_nPack)
          IPPOW = 2**IPSHFT
          do I=1,NUPSG
            IC = CIS%ICase(LUW+I*CIS%nIpWlk)
            ICS = mod(IC/IPPOW,4)
            if (ICS == 0) cycle
            X = CPQ*real((1+ICS)/2,kind=wp)
            ISTA = ISGSTA+I
            call DAXPY_(NDWNSG,X,CI(ISTA),NUPSG,SGM(ISTA),NUPSG)
          end do
        end do
      end do

    else

      ! SPECIAL CASE: WEIGHT OPERATOR, IP=IQ.
      ! IP=IQ < MIDLEV.
      do MVSGM=1,CIS%nMidV
        do ISYUSG=1,SGS%nSym
          NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
          if (NS1 == 0) cycle
          ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
          ISYDSG = Mul(ISYUSG,ISYSGM)
          NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
          IOLW = CIS%IOW(2,ISYDSG,MVSGM)
          IPSHFT = 2*(IP-1)
          LLW = 1+IOLW-CIS%nIpWlk+IPSHFT/(2*CSFG_nPack)
          IPSHFT = mod(IPSHFT,2*CSFG_nPack)
          IPPOW = 2**IPSHFT
          do J=1,NDWNSG
            JC = CIS%ICase(LLW+J*CIS%nIpWlk)
            ICS = mod(JC/IPPOW,4)
            if (ICS == 0) cycle
            X = CPQ*real((1+ICS)/2,kind=wp)
            JSTA = ISGSTA+NUPSG*(J-1)+1
            SGM(JSTA:JSTA+NUPSG-1) = SGM(JSTA:JSTA+NUPSG-1)+X*CI(JSTA:JSTA+NUPSG-1)
          end do
        end do
      end do

    end if
  end if

contains

  ! Two weighted checksums of the CI vector, with distinct, non-repeating
  ! weights 1+frac(i*w), so that vectors with equal component sums (e.g.
  ! E_pq e_I = c e_J) are told apart.
  function CI_Fingerprint(C,n) result(F)
    integer(kind=iwp), intent(in) :: n
    real(kind=wp), intent(in) :: C(n)
    real(kind=wp) :: F(2)
    real(kind=wp), parameter :: W1 = 0.6180339887498949_wp, W2 = 0.4142135623730951_wp
    real(kind=wp) :: X
    integer(kind=iwp) :: k
    F(:) = Zero
    do k=1,n
      X = real(k,kind=wp)
      F(1) = F(1)+C(k)*(One+(X*W1-aint(X*W1)))
      F(2) = F(2)+C(k)*(One+(X*W2-aint(X*W2)))
    end do
  end function CI_Fingerprint

  subroutine sort_icoup_block(ICOUP,NCP,swap,Moved)

    integer(kind=iwp), intent(in) :: NCP
    integer(kind=iwp), intent(inout) :: ICOUP(3,NCP)
    logical(kind=iwp), intent(in) :: swap
    logical(kind=iwp), intent(out) :: Moved
    integer(kind=iwp) :: i, j, temp(3)

    Moved = .false.
    do i=2,NCP
      temp = ICOUP(:,i)
      j = i-1
      do while (j >= 1)
        if (swap) then
          if (ICOUP(1,j) <= temp(1)) exit
        else
          if (ICOUP(2,j) <= temp(2)) exit
        end if
        ICOUP(:,j+1) = ICOUP(:,j)
        j = j-1
        Moved = .true.
      end do
      ICOUP(:,j+1) = temp
    end do
  end subroutine sort_icoup_block

! Reusable construction stage for a crossing excitation operator.
  subroutine Build_Crossing_Temporary( &
    NDWN,NUPC,CI,NUPSG,TMP,NCP,ICOUP,Swap)
    integer(kind=iwp), intent(in) :: NDWN, NUPC, NUPSG, NCP
    integer(kind=iwp), intent(in) :: ICOUP(3,NCP)
    real(kind=wp), intent(in) :: CI(NUPC,NDWN)
    real(kind=wp), intent(out) :: TMP(NUPSG,NDWN)
    logical(kind=iwp), intent(in) :: Swap
    TMP = Zero
    call Apply_row(One,NDWN,NUPC,CI,NUPSG,TMP,NCP,ICOUP,swap=Swap)
  end subroutine Build_Crossing_Temporary

! Reusable contraction stage for a crossing excitation operator.
  subroutine Contract_Crossing_Temporary( &
    Scale,NUP,NDWNC,TMP,NDWNSG,SIGMA,NCP,ICOUP,Swap)
    integer(kind=iwp), intent(in) :: NUP, NDWNC, NDWNSG, NCP
    integer(kind=iwp), intent(in) :: ICOUP(3,NCP)
    real(kind=wp), intent(in) :: Scale, TMP(NUP,NDWNC)
    real(kind=wp), intent(inout) :: SIGMA(NUP,NDWNSG)
    logical(kind=iwp), intent(in) :: Swap
    call Apply_col(Scale,NUP,NDWNC,TMP,NDWNSG,SIGMA,NCP,ICOUP,swap=Swap)
  end subroutine Contract_Crossing_Temporary

  subroutine apply_col(CPQ,NUP,NDWNC,CI,NDWNSG,SIGMA,NCP,ICOUP,swap)

    ! SIGMA(:,I2) = SIGMA(:,I2) + CPQ*V*CI(:,I1) for every coupling (I1,I2,V)
    ! of the lower half (columns of the CSF block). One contiguous column update
    ! per coupling: the blocks are small, so gathering the columns for a GEMM
    ! cost more (copies, temporaries, call overhead) than it gained.
    integer(kind=iwp), intent(in) :: NUP, NDWNC, NDWNSG, NCP, ICOUP(3,NCP)
    real(kind=wp), intent(in) :: CPQ, CI(NUP,NDWNC)
    real(kind=wp), intent(inout) :: SIGMA(NUP,NDWNSG)
    logical(kind=iwp), intent(in) :: swap
    integer(kind=iwp) :: I1, I2, ICP, iFrom, iTo
    real(kind=wp) :: X

    ! with swap the roles of the walk indices are exchanged
    iTo = merge(1,2,swap)
    iFrom = merge(2,1,swap)
    do ICP=1,NCP
      I2 = ICOUP(iTo,ICP)
      I1 = ICOUP(iFrom,ICP)
      X = CPQ*EXS%VTab(ICOUP(3,ICP))
      SIGMA(:,I2) = SIGMA(:,I2)+X*CI(:,I1)
    end do

  end subroutine apply_col

  subroutine apply_row(CPQ,NDWN,NUPC,CI,NUPSG,SIGMA,NCP,ICOUP,swap)

    ! SIGMA(I2,:) = SIGMA(I2,:) + CPQ*V*CI(I1,:) for every coupling (I1,I2,V)
    ! of the upper half (rows of the CSF block). Columns outside, couplings
    ! inside: every access stays within one column of CI and of SIGMA (no
    ! strided rows, no tile copies).
    integer(kind=iwp), intent(in) :: NDWN, NUPC, NUPSG, NCP, ICOUP(3,NCP)
    real(kind=wp), intent(in) :: CPQ, CI(NUPC,NDWN)
    real(kind=wp), intent(inout) :: SIGMA(NUPSG,NDWN)
    logical(kind=iwp), intent(in) :: swap
    integer(kind=iwp) :: ICP, IDWN
    integer(kind=iwp), pointer :: I1LIST(:), I2LIST(:)
    real(kind=wp), pointer :: XLIST(:)

    XLIST => EXS%XLIST
    I1LIST => EXS%I1LIST
    I2LIST => EXS%I2LIST

    if (swap) then
      do ICP=1,NCP
        I1list(ICP) = ICOUP(2,ICP)
        I2list(ICP) = ICOUP(1,ICP)
        Xlist(ICP) = CPQ*EXS%VTab(ICOUP(3,ICP))
      end do
    else
      do ICP=1,NCP
        I1list(ICP) = ICOUP(1,ICP)
        I2list(ICP) = ICOUP(2,ICP)
        Xlist(ICP) = CPQ*EXS%VTab(ICOUP(3,ICP))
      end do
    end if

    do IDWN=1,NDWN
      do ICP=1,NCP
        SIGMA(I2list(ICP),IDWN) = SIGMA(I2list(ICP),IDWN)+Xlist(ICP)*CI(I1list(ICP),IDWN)
      end do
    end do

    nullify(XLIST,I1LIST,I2LIST)

  end subroutine apply_row

end subroutine SG_Epq_Psi

subroutine SG_Apply_Ak_Psi(iState,LSM,T,k,Scale,CI,Sigma)

  use Constants, only: Zero, One
  use stdalloc, only: mma_allocate, mma_deallocate

  implicit none

  integer(kind=iwp), intent(In) :: iState
  integer(kind=iwp), intent(In) :: LSM
  integer(kind=iwp), intent(In) :: k

  real(kind=wp), intent(In) :: T(:,:)
  real(kind=wp), intent(In) :: Scale
  real(kind=wp), intent(In) :: CI(:)

  real(kind=wp), intent(InOut) :: Sigma(:)

  integer(kind=iwp) :: p
  integer(kind=iwp) :: nLev

  real(kind=wp) :: CPK
  real(kind=wp), parameter :: Thr = 1.0e-14_wp

  real(kind=wp), allocatable :: Tmp(:)
  logical(kind=iwp) :: Save_Reuse

  Save_reuse = EXS(istate)%ReUse_SGTMP
  EXS(istate)%ReUse_SGTMP = .false.

  nLev = size(T,1)

  call mma_allocate(Tmp,size(CI),Label='Tmp')

  do p=1,nLev

!
!   Symmetry constraint
!

    if (SGS(iState)%ISM(p) /= SGS(iState)%ISM(k)) cycle

!
!   CPK = Tpk - delta_pk
!

    CPK = T(p,k)

    if (p == k) CPK = CPK-One

    CPK = Scale*CPK

    if (abs(CPK) <= Thr) cycle

    Tmp(:) = Zero

    call SG_Epq_Psi(SGS(iState),CIS(iState),EXS(iState), &
                    p,k,CPK,LSM,CI,Tmp)

    Sigma(:) = Sigma(:)+Tmp(:)

  end do

  call mma_deallocate(Tmp)

  EXS(istate)%ReUse_SGTMP = Save_reuse
end subroutine SG_Apply_Ak_Psi

subroutine SG_EpqEqp_Diag(SGS,CIS,EXS,IP,IQ,CPQ,ISYCI,CI,DIAG)

  use Symmetry_Info, only: Mul
  use Constants, only: Zero, One
  use Definitions, only: wp, iwp, u6

#include "intent.fh"

  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(inout), target :: EXS
  integer(kind=iwp), intent(in) :: IP, IQ, ISYCI
  real(kind=wp), intent(in) :: CPQ
  real(kind=wp), intent(inout) :: CI(*)
  real(kind=wp), intent(_OUT_) :: DIAG(*)
! todo: cleanup
!integer(kind=iwp) :: INDEO, IOC, ISGSTA, ISYDC, ISYDSG, ISYP, ISYPQ, ISYQ, ISYSGM, &
  integer(kind=iwp) :: INDEO, IOC, ISYDC, ISYDSG, ISYP, ISYPQ, ISYQ, ISYSGM, &
                       ISYUC, ISYUSG, LICP, MVSGM, NCP, NDWNC, NDWNSG, NS1, NTMP, NUPC, &
                       NUPSG, NCP1, NCP2, MV, MVX, nCSFs, LICP1, LICP2

! declarations to facilitate the reuse of the intermediates in SGTMP (the key
! is kept in EXS; here they depend only on the couplings, not on CI)
  logical(kind=iwp) :: Reuse_Sigma, Swap, Test1, Test2
  integer(kind=iwp) :: iOff, jOff, JP, JQ

!*********************************************************************
!  GIVEN ACTIVE LEVEL INDICES IP AND IQ, AND INPUT CI ARRAYS
!  CI AND SGM, THIS ROUTINE ADDS TO SGM THE RESULT OF ACTING ON
!  CI WITH THE NUMBER CPQ TIMES THE EXCITATION OPERATOR E(IP,IQ).
!  THE ADDITIONAL ENTRIES IN THE PARAMETER LIST ARE TABLES THAT
!  WERE PREPARED BY GINIT AND ITS SUBROUTINES.
!*********************************************************************

  nCSFs = CIS%NCSF(ISYCI)
! SYMMETRY OF ORBITALS:
  ISYP = SGS%ISm(IP)
  ISYQ = SGS%ISm(IQ)
  ISYPQ = Mul(ISYP,ISYQ)
! SYMMETRY OF SIGMA ARRAY:
  ISYSGM = Mul(ISYPQ,ISYCI)

  Swap = IQ > IP
  Test1 = merge(IQ <= SGS%MidLev,IP <= SGS%MidLev,Swap)
  Test2 = merge(SGS%MidLev < IP,SGS%MidLev < IQ,Swap)

  if (IQ /= IP) then
    Swap = IQ > IP
    Test1 = merge(IQ <= SGS%MidLev,IP <= SGS%MidLev,Swap)
    Test2 = merge(SGS%MidLev < IP,SGS%MidLev < IQ,Swap)
    if (Test1) then
      ! EXCITING CASE, IQ<IP<=MIDLEV.
      INDEO = merge(2*SGS%nLev+(IQ*(IQ-1))/2+IP,2*SGS%nLev+(IP*(IP-1))/2+IQ,Swap)
      do MVSGM=1,CIS%nMidV
        do ISYUSG=1,SGS%nSym

          NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
          if (NS1 == 0) cycle
          ISYDSG = Mul(ISYUSG,ISYSGM)           ! compute the lower symmetry
          ISYDC = Mul(ISYPQ,ISYDSG)             ! compute the symmetry of the sigma vector
          NDWNC = CIS%NOW(2,ISYDC,MVSGM)
          if (NDWNC == 0) cycle
!       ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM) ! get the off-set to the sigma vector block
          NUPSG = CIS%NOW(1,ISYUSG,MVSGM)           ! number of upper half-walks by symmetry and midvertex.
          IOC = CIS%IOCSF(ISYUSG,MVSGM,ISYCI)       ! get the off-set to the CI vector block
          NCP = merge(EXS%NOCP(INDEO,ISYDSG,MVSGM),EXS%NOCP(INDEO,ISYDC,MVSGM),SWAP)
          if (NCP == 0) cycle
          NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)

          ! off-set to the block of compressed coupling coefficients.
          LICP = merge(EXS%IOCP(INDEO,ISYDSG,MVSGM),EXS%IOCP(INDEO,ISYDC,MVSGM),SWAP)
          ! CASE IS: LOWER HALF, EXCITE:

!       call sort_icoup_block(EXS%ICOUP(1,LICP+1),NCP,Swap)
!       call Apply_col(CPQ,NUPSG,NDWNC,CI(IOC+1),NDWNSG,SGM(ISGSTA+1),NCP,EXS%ICOUP(1,LICP+1),swap=Swap)
          call Apply_colX(CPQ,NUPSG,NDWNC,DIAG(IOC+1),NCP,EXS%ICOUP(1,LICP+1),swap=Swap)

        end do
      end do
    else if (Test2) then
      ! EXCITING CASE, MIDLEV<IQ<IP
      INDEO = merge(2*SGS%nLev+(IQ*(IQ-1))/2+IP,2*SGS%nLev+(IP*(IP-1))/2+IQ,SWAP)
      do MVSGM=1,CIS%nMidV
        do ISYUSG=1,SGS%nSym
          NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
          if (NS1 == 0) cycle
          ISYUC = Mul(ISYPQ,ISYUSG)
          NUPC = CIS%NOW(1,ISYUC,MVSGM)
          if (NUPC == 0) cycle
!       ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
          ISYDSG = Mul(ISYUSG,ISYSGM)
          NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
          IOC = CIS%IOCSF(ISYUC,MVSGM,ISYCI)
          NCP = merge(EXS%NOCP(INDEO,ISYUSG,MVSGM),EXS%NOCP(INDEO,ISYUC,MVSGM),SWAP)
          if (NCP == 0) cycle
          LICP = merge(EXS%IOCP(INDEO,ISYUSG,MVSGM),EXS%IOCP(INDEO,ISYUC,MVSGM),SWAP)
          ! CASE IS: UPPER HALF, EXCITE:
!       call apply_row(CPQ,NDWNSG,NUPC,CI(IOC+1),NUPSG,SGM(ISGSTA+1),NCP,EXS%ICOUP(1,LICP+1),swap=Swap)
          call apply_rowX(CPQ,NDWNSG,NUPC,DIAG(IOC+1),NCP,EXS%ICOUP(1,LICP+1),swap=Swap)
        end do
      end do
    else

      ! EXCITING CASE, IQ<=MIDLEV<IP

      iOff = 1
      JP = merge(IQ,IP,SWAP)
      JQ = merge(IP,IQ,SWAP)
      Reuse_Sigma = EXS%Reuse_SGTMP .and. EXS%Save_Valid .and. (EXS%Save_Owner == 2) .and. (EXS%Save_Swap .eqv. Swap) .and. &
                    (EXS%Save_U == JP) .and. (EXS%Save_LSym == SGS%ISm(JQ)) .and. (EXS%Save_CISym == ISYCI)

      do MVSGM=1,CIS%nMidV
        do MV=1,2
          MVX = merge(EXS%MVR(MVSGM,MV),EXS%MVL(MVSGM,MV),SWAP)
          if (MVX == 0) cycle
          do ISYUSG=1,SGS%nSym
            NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
            if (NS1 == 0) cycle
!         ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
            NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
            ISYDSG = Mul(ISYUSG,ISYSGM)
            ISYUC = merge(Mul(ISYQ,ISYUSG),Mul(ISYP,ISYUSG),SWAP)
            ISYDC = merge(Mul(ISYP,ISYDSG),Mul(ISYQ,ISYDSG),SWAP)

            NUPC = CIS%NOW(1,ISYUC,MVX)
            if (NUPC == 0) cycle
            NDWNC = CIS%NOW(2,ISYDC,MVX)
            if (NDWNC == 0) cycle

            INDEO = merge(merge(IQ,SGS%nLev+IQ,MV == 1),merge(IP,SGS%nLev+IP,MV == 1),SWAP)
            NCP1 = merge(EXS%NOCP(INDEO,ISYUSG,MVSGM),EXS%NOCP(INDEO,ISYUC,MVX),SWAP)
            if (NCP1 == 0) cycle

            ! CASE IS: UPPER HALF, EXCITE:
            LICP1 = merge(EXS%IOCP(INDEO,ISYUSG,MVSGM),EXS%IOCP(INDEO,ISYUC,MVX),Swap)
            IOC = CIS%IOCSF(ISYUC,MVX,ISYCI)

            ! IN CASE OF REUSE COMPUTE THE TEMPORARY SIGMA VECTOR REGARDLESS OF THE NCP2 VALUE.
            if (EXS%Reuse_SGTMP .and. (.not. Reuse_Sigma)) then
              NTMP = NUPC
              EXS%SGTMP(iOff:iOff+NTMP-1) = Zero
!            call Apply_row(One,NDWNC,NUPC,CI(IOC+1),NUPSG,EXS%SGTMP(iOff),NCP1,EXS%ICOUP(1,LICP1+1),swap=Swap)
              call Apply_rowY(One,NUPC,EXS%SGTMP(iOff),NCP1,EXS%ICOUP(1,LICP1+1),swap=Swap)
            end if

            jOff = iOff
            if (EXS%Reuse_SGTMP) then
              iOff = iOff+NUPC
              if (iOff-1 > size(EXS%SGTMP)) then
                write(u6,*) 'iOff-1>SIZE(EXS%SGTMP) -- 1'
                write(u6,*) iOff-1,size(EXS%SGTMP)
                call Abend()
              end if
            end if

            INDEO = merge(merge(IP,SGS%nLev+IP,MV == 1),merge(IQ,SGS%nLev+IQ,MV == 1),SWAP)
            NCP2 = merge(EXS%NOCP(INDEO,ISYDSG,MVSGM),EXS%NOCP(INDEO,ISYDC,MVX),SWAP)
            if (NCP2 == 0) cycle

            if (.not. EXS%Reuse_SGTMP) then
              NTMP = NUPC
              EXS%SGTMP(jOff:jOff+NTMP-1) = Zero
!            call Apply_row(One,NDWNC,NUPC,CI(IOC+1),NUPSG,EXS%SGTMP(jOff),NCP1,EXS%ICOUP(1,LICP1+1),swap=Swap)
              call Apply_rowY(One,NUPC,EXS%SGTMP(jOff),NCP1,EXS%ICOUP(1,LICP1+1),swap=Swap)
            end if

            ! CASE IS: LOWER HALF, EXCITE:
            NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
            LICP2 = merge(EXS%IOCP(INDEO,ISYDSG,MVSGM),EXS%IOCP(INDEO,ISYDC,MVX),SWAP)
!         call Apply_col(CPQ,NUPSG,NDWNC,EXS%SGTMP(jOff),NDWNSG,SGM(ISGSTA+1),NCP2,EXS%ICOUP(1,LICP2+1),swap=Swap)
            call Apply_colY(CPQ,NUPC,EXS%SGTMP(jOff),NDWNC,DIAG(IOC+1),NCP2,EXS%ICOUP(1,LICP2+1),swap=Swap)

          end do
        end do
      end do
      if (EXS%Reuse_SGTMP) then
        if (.not. Reuse_Sigma) then
          EXS%Save_Valid = .true.
          EXS%Save_Owner = 2
          EXS%Save_Swap = Swap
          EXS%Save_U = JP
          EXS%Save_LSym = SGS%ISm(JQ)
          EXS%Save_CISym = ISYCI
          EXS%Save_nCSF = nCSFs
        end if
      else
        ! SGTMP was overwritten
        EXS%Save_Valid = .false.
      end if
    end if
  else
    CI(1:nCSFs) = One
    call SG_Epq_Psi(SGS,CIS,EXS,IP,IQ,CPQ,ISYCI,CI,DIAG)
    DIAG(1:nCSFs) = DIAG(1:nCSFs)**2
  end if

contains

  subroutine apply_colX(CPQ,NUP,NDWNSG,SIGMA,NCP,ICOUP,swap)
    integer(kind=iwp), intent(in) :: NUP, NDWNSG, NCP, ICOUP(3,NCP)
    real(kind=wp), intent(in) :: CPQ
    real(kind=wp), intent(inout) :: SIGMA(NUP,NDWNSG)
    logical(kind=iwp), intent(in) :: swap   ! =.false. (EXC1), =.true. (DEX1)

    integer(kind=iwp) :: ICP, JLFT, JRGT, IUP, I1, idx
    real(kind=wp) :: X
    real(kind=wp), pointer :: VTAB(:)

    vTab => EXS%VTab
    do ICP=1,NCP
      JLFT = ICOUP(1,ICP)
      JRGT = ICOUP(2,ICP)

      I1 = merge(JRGT,JLFT,Swap)

      idx = ICOUP(3,ICP)
      X = (CPQ*VTab(idx))**2
!$OMP SIMD
      do IUP=1,NUP
        SIGMA(IUP,I1) = SIGMA(IUP,I1)+X
      end do
    end do
    vTab => null()
  end subroutine apply_colX

  subroutine apply_rowX(CPQ,NDWN,NUPSG,SIGMA,NCP,ICOUP,swap)
    integer(kind=iwp), intent(in) :: NDWN, NUPSG, NCP, ICOUP(3,NCP)
    real(kind=wp), intent(in) :: CPQ
    real(kind=wp), intent(inout) :: SIGMA(NUPSG,NDWN)
    logical(kind=iwp), intent(in) :: swap   ! =.false. (EXC2), =.true. (DEX2)

    integer(kind=iwp) :: ICP, ILFT, IRGT, IDWN, I1, idx
    real(kind=wp) :: X
    real(kind=wp), pointer :: VTAB(:)

    vTab => EXS%VTab
    do ICP=1,NCP
      ILFT = ICOUP(1,ICP)
      IRGT = ICOUP(2,ICP)

      I1 = merge(IRGT,ILFT,Swap)

      idx = ICOUP(3,ICP)
      X = (CPQ*VTab(idx))**2
!$OMP SIMD
      do IDWN=1,NDWN
        SIGMA(I1,IDWN) = SIGMA(I1,IDWN)+X
      end do
    end do
    vTab => null()
  end subroutine apply_rowX

  subroutine apply_rowY(CPQ,NUPSG,UP,NCP,ICOUP,swap)
    integer(kind=iwp), intent(in) :: NUPSG, NCP, ICOUP(3,NCP)
    real(kind=wp), intent(in) :: CPQ
    real(kind=wp), intent(inout) :: UP(NUPSG)
    logical(kind=iwp), intent(in) :: swap   ! =.false. (EXC2), =.true. (DEX2)

! todo: cleanup
! integer(kind=iwp) :: ICP, ILFT, IRGT, I1, I2, idx
    integer(kind=iwp) :: ICP, ILFT, IRGT, I1, idx
    real(kind=wp) :: X
    real(kind=wp), pointer :: VTAB(:)

    vTab => EXS%VTab
    do ICP=1,NCP
      ILFT = ICOUP(1,ICP)
      IRGT = ICOUP(2,ICP)

      I1 = merge(IRGT,ILFT,SWAP)
! I2=MERGE(ILFT,IRGT,SWAP)

      idx = ICOUP(3,ICP)
      X = (CPQ*VTab(idx))**2
      UP(I1) = UP(I1)+X
    end do
    vTab => null()
  end subroutine apply_rowY

  subroutine apply_colY(CPQ,NUP,Up,NDWNSG,SIGMA,NCP,ICOUP,swap)
    integer(kind=iwp), intent(in) :: NUP, NDWNSG, NCP, ICOUP(3,NCP)
    real(kind=wp), intent(in) :: CPQ, Up(NUP)
    real(kind=wp), intent(inout) :: SIGMA(NUP,NDWNSG)
    logical(kind=iwp), intent(in) :: swap   ! =.false. (EXC1), =.true. (DEX1)

! todo: cleanup
! integer(kind=iwp) :: ICP, JLFT, JRGT, IUP, I1, I2, idx
    integer(kind=iwp) :: ICP, JLFT, JRGT, IUP, I1, idx
    real(kind=wp) :: X
    real(kind=wp), pointer :: VTAB(:)

    vTab => EXS%VTab
    do ICP=1,NCP
      JLFT = ICOUP(1,ICP)
      JRGT = ICOUP(2,ICP)

      I1 = merge(JRGT,JLFT,Swap)
! I2=MERGE(JLFT,JRGT,Swap)

      idx = ICOUP(3,ICP)
      X = (CPQ*VTab(idx))**2
!$OMP SIMD
      do IUP=1,NUP
        SIGMA(IUP,I1) = SIGMA(IUP,I1)+X*Up(IUP)
      end do
    end do
    vTab => null()
  end subroutine apply_colY

end subroutine SG_EpqEqp_Diag

subroutine SG_TRACI(iState,LSM,T,nT,CI,nCSF)

  use molcas, only: MxAct
  use Constants, only: Zero, Half, Three
  use stdalloc, only: mma_allocate, mma_deallocate

  implicit none

  integer(kind=iwp), intent(In) :: iState
  integer(kind=iwp), intent(In) :: LSM
  integer(kind=iwp), intent(In) :: nT
  integer(kind=iwp), intent(In) :: nCSF

  real(kind=wp), intent(In), target :: T(nT)
  real(kind=wp), intent(InOut) :: CI(nCSF)

  integer(kind=iwp) :: k, iOff, iSym, nAct, p, kLocal, pLocal
  integer(kind=iwp) :: nLev, Level(MxAct)

  real(kind=wp) :: tkk, X

  real(kind=wp), allocatable :: Tmp(:)
  real(kind=wp), allocatable :: T2(:,:)

  nLev = SGS(iState)%NLEV

! Call RecPrt('SG_TRACI; CI',' ',CI,1,nCSF)
  call mma_allocate(T2,nLev,nLev,Label='T2')
! Call RecPrt('T',' ',T,1,nT)

  T2(:,:) = Zero
  iOff = 0

  do iSym=1,SGS(iState)%nSym

    nAct = 0

    do p=1,nLev
      if (SGS(iState)%ISM(p) == iSym) then
        nAct = nAct+1
        Level(nAct) = p
      end if
    end do

    do kLocal=1,nAct
      k = Level(kLocal)

      do pLocal=1,nAct
        p = Level(pLocal)
        T2(p,k) = T(iOff+pLocal+(kLocal-1)*nAct)
      end do
    end do
    iOff = iOff+nAct*nAct
  end do
! Call RecPrt('T2',' ',T2,nLev,nLev)
!
! Work array
!
  call mma_allocate(Tmp,nCSF,Label='tmp')
!
! Apply row operators sequentially.
!

  do k=1,nLev

!    Write (6,*) 'ik=',k
    TMP(:) = Zero

    ! tmp = 1/2 A_k CI
    call SG_Apply_Ak_Psi(iState,LSM,T2,k,Half,CI,TMP)
!    Call RecPrt('TMP',' ',TMP,1,nCSF)

    tkk = T2(k,k)
    X = Three-tkk
!    Write (6,*) 'Ckk, X=',TKK, X
    CI(:) = CI(:)+X*TMP(:)
!    Call RecPrt('CI',' ',CI,1,nCSF)

    call SG_Apply_Ak_Psi(iState,LSM,T2,k,One,TMP,CI)
!    Call RecPrt('CI',' ',CI,1,nCSF)

  end do

!
! Cleanup
!
  call mma_deallocate(Tmp)
  call mma_deallocate(T2)

end subroutine SG_TRACI

end module SGUGA_Epq
