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

!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*
subroutine SG_Epq_Psi(SGS,CIS,EXS,IP,IQ,CPQ,ISYCI,CI,SGM)

use Symmetry_Info, only: Mul
use sguga, only: CIStruct, EXStruct, SGStruct, nPack
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp), intent(in) :: IP, IQ, ISYCI
real(kind=wp), intent(in) :: CPQ, CI(*)
real(kind=wp), intent(_OUT_) :: SGM(*)
integer(kind=iwp) :: I, IC, ICS, INDEO, IOC, IOLW, IOUW, IPPOW, IPSHFT, ISGSTA, ISTA, ISYDC, ISYDSG, ISYP, ISYPQ, ISYQ, ISYSGM, &
                     ISYUC, ISYUSG, J, JC, JSTA, LICP, LLW, LUW, MV1, MV2, MV4, MV5, MVSGM, NCP, NDWNC, NDWNSG, NS1, NTMP, NUPC, &
                     NUPSG, NCP1, NCP2, MV, MVX
real(kind=wp) :: X

!***********************************************************************
!  GIVEN ACTIVE LEVEL INDICES IP AND IQ, AND INPUT CI ARRAYS
!  CI AND SGM, THIS ROUTINE ADDS TO SGM THE RESULT OF ACTING ON
!  CI WITH THE NUMBER CPQ TIMES THE EXCITATION OPERATOR E(IP,IQ).
!  THE ADDITIONAL ENTRIES IN THE PARAMETER LIST ARE TABLES THAT
!  WERE PREPARED BY GINIT AND ITS SUBROUTINES.
!***********************************************************************

! SYMMETRY OF ORBITALS:
ISYP = SGS%ISm(IP)
ISYQ = SGS%ISm(IQ)
ISYPQ = Mul(ISYP,ISYQ)
! SYMMETRY OF SIGMA ARRAY:
ISYSGM = Mul(ISYPQ,ISYCI)
! SYMMETRY OF SIGMA ARRAY:
ISYSGM = Mul(ISYPQ,ISYCI)

if (IQ < IP) then
  ! THEN THIS IS AN EXCITING OPERATOR.
  if (IP <= SGS%MidLev) then

    ! EXCITING CASE, IQ<IP<=MIDLEV.
    INDEO = 2*SGS%nLev+(IP*(IP-1))/2+IQ
    do MVSGM=1,CIS%nMidV
      do ISYUSG=1,SGS%nSym

        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)  ! Number of  CSFs for a given upper symmetry, midvertex, and total symmetry
        if (NS1 == 0) cycle
        ISYDSG = Mul(ISYUSG,ISYSGM)           ! compute the lower symmetry
        ISYDC = Mul(ISYPQ,ISYDSG)             ! compute the symmetry of the sigma vector
        NDWNC = CIS%NOW(2,ISYDC,MVSGM)        ! number of lower half-walks by symmetry and midvertex.
        if (NDWNC == 0) cycle
        ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM) ! get the off-set to the sigma vector block
        NUPSG = CIS%NOW(1,ISYUSG,MVSGM)           ! number of upper half-walks by symmetry and midvertex.
        IOC = CIS%IOCSF(ISYUSG,MVSGM,ISYCI)       ! get the off-set to the CI vector block
        NCP = EXS%NOCP(INDEO,ISYDC,MVSGM)         ! get the number of compressed coupling coefficients in the block
        if (NCP > 0) then
          LICP = EXS%IOCP(INDEO,ISYDC,MVSGM)      ! get the off-set to the block of compressed coupling coefficients.
          ! CASE IS: LOWER HALF, EXCITE:
          call EXC1(CPQ,NUPSG,CI(IOC+1),SGM(ISGSTA+1),NCP,EXS%ICOUP(1,LICP+1))
        end if

      end do
    end do

  else if (SGS%MidLev < IQ) then

    ! EXCITING CASE, MIDLEV<IQ<IP
    INDEO = 2*SGS%nLev+(IP*(IP-1))/2+IQ
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
        NCP = EXS%NOCP(INDEO,ISYUC,MVSGM)
        if (NCP == 0) cycle
        LICP = EXS%IOCP(INDEO,ISYUC,MVSGM)
        ! CASE IS: UPPER HALF, EXCITE:
        call EXC2(CPQ,NDWNSG,NUPC,CI(IOC+1),NUPSG,SGM(ISGSTA+1),NCP,EXS%ICOUP(1,LICP+1))
      end do
    end do

  else

    ! EXCITING CASE, IQ<=MIDLEV<IP
    do MVSGM=1,CIS%nMidV
      do MV = 1, 2
        MVX = EXS%MVL(MVSGM,MV)
        if (MVX == 0) cycle
        do ISYUSG=1,SGS%nSym
          NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
          if (NS1 == 0) cycle
          ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
          ISYDSG = Mul(ISYUSG,ISYSGM)
          ISYUC = Mul(ISYP,ISYUSG)
          ISYDC = Mul(ISYQ,ISYDSG)

          NUPC = CIS%NOW(1,ISYUC,MVX)
          if (NUPC == 0) cycle
          NDWNC = CIS%NOW(2,ISYDC,MVX)
          if (NDWNC == 0) cycle

          if (MV==1) then
            INDEO = IP
          else
            INDEO = SGS%nLev+IP
          end if
          NCP1 = EXS%NOCP(INDEO,ISYUC,MVX)
          if (NCP1 == 0) cycle
          NTMP = NUPSG*NDWNC
          EXS%SGTMP(1:NTMP) = Zero
          LICP = EXS%IOCP(INDEO,ISYUC,MVX)
          IOC = CIS%IOCSF(ISYUC,MVX,ISYCI)
          if (MV==1) then
            INDEO = IQ
          else
            INDEO = SGS%nLev+IQ
          end if
          NCP2 = EXS%NOCP(INDEO,ISYDC,MVX)
          if (NCP2 == 0) cycle
          ! CASE IS: UPPER HALF, EXCITE:
          call EXC2(CPQ,NDWNC,NUPC,CI(IOC+1),NUPSG,EXS%SGTMP,NCP1,EXS%ICOUP(1,LICP+1))
          LICP = EXS%IOCP(INDEO,ISYDC,MVX)
          ! CASE IS: LOWER HALF, EXCITE:
          call EXC1(One,NUPSG,EXS%SGTMP,SGM(ISGSTA+1),NCP2,EXS%ICOUP(1,LICP+1))

        end do
      end do
    end do

  end if
else if (IP < IQ) then
  ! THEN THIS IS A DEEXCITING OPERATOR.
  if (IQ <= SGS%MidLev) then

    ! DEEXCITING OPERATOR, IP<IQ<=MIDLEV.
    INDEO = 2*SGS%nLev+(IQ*(IQ-1))/2+IP
    do MVSGM=1,CIS%nMidV
      do ISYUSG=1,SGS%nSym
        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
        if (NS1 == 0) cycle
        ISYDSG = Mul(ISYUSG,ISYSGM)
        ISYDC = Mul(ISYPQ,ISYDSG)
        NDWNC = CIS%NOW(2,ISYDC,MVSGM)
        if (NDWNC == 0) cycle
        ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
        NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
        IOC = CIS%IOCSF(ISYUSG,MVSGM,ISYCI)
        NCP = EXS%NOCP(INDEO,ISYDSG,MVSGM)
        if (NCP == 0) cycle
        LICP = EXS%IOCP(INDEO,ISYDSG,MVSGM)
        ! CASE IS: LOWER HALF, DEEXCITE:
        call DEX1(CPQ,NUPSG,CI(IOC+1),SGM(ISGSTA+1),NCP,EXS%ICOUP(1,LICP+1))
      end do
    end do

  else if (SGS%MidLev < IP) then

    ! DEEXCITING OPERATOR, MIDLEV<IP<IQ
    INDEO = 2*SGS%nLev+(IQ*(IQ-1))/2+IP
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
        NCP = EXS%NOCP(INDEO,ISYUSG,MVSGM)
        if (NCP == 0) cycle
        LICP = EXS%IOCP(INDEO,ISYUSG,MVSGM)
        ! CASE IS: UPPER HALF, DEEXCITE:
        call DEX2(CPQ,NDWNSG,NUPC,CI(IOC+1),NUPSG,SGM(ISGSTA+1),NCP,EXS%ICOUP(1,LICP+1))
      end do
    end do

  else

    ! DEEXCITING CASE, IP<=MIDLEV<IQ.
    do MVSGM=1,CIS%nMidV
      do MV = 1, 2
         MVX = EXS%MVR(MVSGM,MV)
         if (MVX == 0) cycle

         do ISYUSG=1,SGS%nSym
           NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
           if (NS1 == 0) cycle
           ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
           NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
           ISYDSG = Mul(ISYUSG,ISYSGM)
           ISYUC = Mul(ISYQ,ISYUSG)
           ISYDC = Mul(ISYP,ISYDSG)

           NUPC = CIS%NOW(1,ISYUC,MVX)
           if (NUPC == 0) cycle
           NDWNC = CIS%NOW(2,ISYDC,MVX)
           if (NDWNC == 0) cycle

           if (MV==1) then
             INDEO = IQ
           else
             INDEO = SGS%nLev+IQ
           endif
           NCP1 = EXS%NOCP(INDEO,ISYUSG,MVSGM)
           if (NCP1 == 0) cycle
           NTMP = NUPSG*NDWNC
           EXS%SGTMP(1:NTMP) = Zero
           LICP = EXS%IOCP(INDEO,ISYUSG,MVSGM)
           IOC = CIS%IOCSF(ISYUC,MVX,ISYCI)
           if (MV==1) then
             INDEO = IP
           else
             INDEO = SGS%nLev+IP
           endif
           NCP2 = EXS%NOCP(INDEO,ISYDSG,MVSGM)
           if (NCP2 == 0) cycle
           ! CASE IS: UPPER HALF, DEEXCITE:
           call DEX2(CPQ,NDWNC,NUPC,CI(IOC+1),NUPSG,EXS%SGTMP,NCP1,EXS%ICOUP(1,LICP+1))
           LICP = EXS%IOCP(INDEO,ISYDSG,MVSGM)
           ! CASE IS: LOWER HALF, DEEXCITE:
           call DEX1(One,NUPSG,EXS%SGTMP,SGM(ISGSTA+1),NCP2,EXS%ICOUP(1,LICP+1))

        end do

      end do
    end do

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
        LUW = 1+IOUW-CIS%nIpWlk+IPSHFT/(2*nPack)
        IPSHFT = mod(IPSHFT,2*nPack)
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
        LLW = 1+IOLW-CIS%nIpWlk+IPSHFT/(2*nPack)
        IPSHFT = mod(IPSHFT,2*nPack)
        IPPOW = 2**IPSHFT
        do J=1,NDWNSG
          JC = CIS%ICase(LLW+J*CIS%nIpWlk)
          ICS = mod(JC/IPPOW,4)
          if (ICS == 0) cycle
          X = CPQ*real((1+ICS)/2,kind=wp)
          JSTA = ISGSTA + NUPSG*(J-1) + 1
          SGM(JSTA:JSTA+NUPSG-1) = SGM(JSTA:JSTA+NUPSG-1)+X*CI(JSTA:JSTA+NUPSG-1)
        end do
      end do
    end do

  end if
end if

contains

subroutine EXC1(CPQ,NUP,A,B,NCP,ICOUP)

  integer(kind=iwp), intent(in) :: NUP, NCP, ICOUP(3,NCP)
  real(kind=wp), intent(in) :: CPQ, A(NUP,*)
  real(kind=wp), intent(inout) :: B(NUP,*)
  integer(kind=iwp) :: ICP, JLFT, JRGT, IUP
  real(kind=wp) :: X

  ! CASE: ADD EPQ*A TO B, WHERE Q<P<=MIDLEV AND A AND B ARE SINGLE
  ! MATRIX BLOCKS.
  ! ONLY LOWER WALKS AFFECTED => COLUMN OPERATION.
  ! P>Q, EXCITING OPERATOR  => A IS  LEFTHAND WALK.
  ! USE VECTOR ROUTINE CALL IF LONG ENOUGH:
  if (NUP > 20) then
    do ICP=1,NCP
      JLFT = ICOUP(1,ICP)
      JRGT = ICOUP(2,ICP)
      X = CPQ*EXS%VTab(ICOUP(3,ICP))
      call DAXPY_(NUP,X,A(:,JLFT),1,B(:,JRGT),1)
    end do
  else
    do ICP=1,NCP
      JLFT = ICOUP(1,ICP)
      JRGT = ICOUP(2,ICP)
      X = CPQ*EXS%VTab(ICOUP(3,ICP))
!!$OMP SIMD
      do IUP=1,NUP
         B(IUP,JRGT) = B(IUP,JRGT)+X*A(IUP,JLFT)
      end do
    end do
  end if

end subroutine EXC1

subroutine EXC2(CPQ,NDWN,NUPA,A,NUPB,B,NCP,ICOUP)

  integer(kind=iwp), intent(in) :: NDWN, NUPA, NUPB, NCP, ICOUP(3,NCP)
  real(kind=wp), intent(in) :: CPQ, A(NUPA,NDWN)
  real(kind=wp), intent(inout) :: B(NUPB,NDWN)
  integer(kind=iwp) :: ICP, ILFT, IRGT, IDWN
  real(kind=wp) :: X

  ! CASE: ADD EPQ*A TO B, WHERE MIDLEV<Q<P AND A AND B ARE SINGLE
  ! MATRIX BLOCKS.
  ! ONLY UPPER WALKS AFFECTED =>    ROW OPERATION.
  ! Q<P, EXCITING OPERATOR  => A IS  LEFTHAND WALK.
  if (NDWN > 20) then
    do ICP=1,NCP
      ILFT = ICOUP(1,ICP)
      IRGT = ICOUP(2,ICP)
      X = CPQ*EXS%VTab(ICOUP(3,ICP))
      call DAXPY_(NDWN,X,A(ILFT,1),NUPA,B(IRGT,1),NUPB)
    end do
  else
    do ICP=1,NCP
      ILFT = ICOUP(1,ICP)
      IRGT = ICOUP(2,ICP)
      X = CPQ*EXS%VTab(ICOUP(3,ICP))
!!$OMP SIMD
      do IDWN=1,NDWN
         B(IRGT,IDWN) = B(IRGT,IDWN)+X*A(ILFT,IDWN)
      end do
    end do
  end if

end subroutine EXC2

subroutine DEX1(CPQ,NUP,A,B,NCP,ICOUP)

  integer(kind=iwp), intent(in) :: NUP, NCP, ICOUP(3,NCP)
  real(kind=wp), intent(in) :: CPQ, A(NUP,*)
  real(kind=wp), intent(inout) :: B(NUP,*)
  integer(kind=iwp) :: ICP, JLFT, JRGT, IUP
  real(kind=wp) :: X

  ! CASE: ADD EPQ*A TO B, WHERE P<Q<=MIDLEV AND A AND B ARE SINGLE
  ! MATRIX BLOCKS.
  ! ONLY LOWER WALKS AFFECTED => COLUMN OPERATION.
  ! P<Q, DEEXCITING OPERATOR  => A IS RIGHTHAND WALK.
  ! CHOSE DAXPY IF LARGE VECTOR LENGTH:
  if (NUP > 20) then
    do ICP=1,NCP
      JLFT = ICOUP(1,ICP)
      JRGT = ICOUP(2,ICP)
      X = CPQ*EXS%VTab(ICOUP(3,ICP))
      call DAXPY_(NUP,X,A(:,JRGT),1,B(:,JLFT),1)
    end do
  else
    do ICP=1,NCP
      JLFT = ICOUP(1,ICP)
      JRGT = ICOUP(2,ICP)
      X = CPQ*EXS%VTab(ICOUP(3,ICP))
!!$OMP SIMD
      do IUP=1,NUP
         B(IUP,JLFT) = B(IUP,JLFT)+X*A(IUP,JRGT)
      end do
    end do
  end if

end subroutine DEX1

subroutine DEX2(CPQ,NDWN,NUPA,A,NUPB,B,NCP,ICOUP)

  integer(kind=iwp), intent(in) :: NDWN, NUPA, NUPB, NCP, ICOUP(3,NCP)
  real(kind=wp), intent(in) :: CPQ, A(NUPA,NDWN)
  real(kind=wp), intent(inout) :: B(NUPB,NDWN)
  integer(kind=iwp) :: ICP, ILFT, IRGT, IDWN
  real(kind=wp) :: X

  ! CASE: ADD EPQ*A TO B, WHERE MIDLEV<P<Q AND A AND B ARE SINGLE
  ! MATRIX BLOCKS.
  ! ONLY UPPER WALKS AFFECTED =>    ROW OPERATION.
  ! P<Q, DEEXCITING OPERATOR  => A IS RIGHTHAND WALK.
  ! CHOSE DAXPY IF LARGE LENGTH.
  if (NDWN > 20) then
    do ICP=1,NCP
      ILFT = ICOUP(1,ICP)
      IRGT = ICOUP(2,ICP)
      X = CPQ*EXS%VTab(ICOUP(3,ICP))
      call DAXPY_(NDWN,X,A(IRGT,1),NUPA,B(ILFT,1),NUPB)
    end do
  else
    do ICP=1,NCP
      ILFT = ICOUP(1,ICP)
      IRGT = ICOUP(2,ICP)
      X = CPQ*EXS%VTab(ICOUP(3,ICP))
!!$OMP SIMD
      do IDWN=1,NDWN
         B(ILFT,IDWN) = B(ILFT,IDWN)+X*A(IRGT,IDWN)
      end do
    end do
  end if

end subroutine DEX2

end subroutine SG_Epq_Psi
