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
type(EXStruct), intent(inout), target :: EXS
integer(kind=iwp), intent(in) :: IP, IQ, ISYCI
real(kind=wp), intent(in) :: CPQ, CI(*)
real(kind=wp), intent(_OUT_) :: SGM(*)
integer(kind=iwp) :: I, IC, ICS, INDEO, IOC, IOLW, IOUW, IPPOW, IPSHFT, ISGSTA, ISTA, ISYDC, ISYDSG, ISYP, ISYPQ, ISYQ, ISYSGM, &
                     ISYUC, ISYUSG, J, JC, JSTA, LICP, LLW, LUW, MVSGM, NCP, NDWNC, NDWNSG, NS1, NTMP, NUPC, &
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

if (IQ < IP) then
  ! THEN THIS IS AN EXCITING OPERATOR.
  if (IP <= SGS%MidLev) then

    ! EXCITING CASE, IQ<IP<=MIDLEV.
    INDEO = 2*SGS%nLev+(IP*(IP-1))/2+IQ
    do MVSGM=1,CIS%nMidV
      do ISYUSG=1,SGS%nSym

        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM) ; if (NS1 == 0) cycle
        ISYDSG = Mul(ISYUSG,ISYSGM)           ! compute the lower symmetry
        ISYDC = Mul(ISYPQ,ISYDSG)             ! compute the symmetry of the sigma vector
        NDWNC = CIS%NOW(2,ISYDC,MVSGM)       ; if (NDWNC == 0) cycle
        ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM) ! get the off-set to the sigma vector block
        NUPSG = CIS%NOW(1,ISYUSG,MVSGM)           ! number of upper half-walks by symmetry and midvertex.
        IOC = CIS%IOCSF(ISYUSG,MVSGM,ISYCI)       ! get the off-set to the CI vector block
        NCP = EXS%NOCP(INDEO,ISYDC,MVSGM)    ; if (NCP == 0) cycle
        NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
        LICP = EXS%IOCP(INDEO,ISYDC,MVSGM)      ! get the off-set to the block of compressed coupling coefficients.
        ! CASE IS: LOWER HALF, EXCITE:
        call sort_icoup_block(EXS%ICOUP(1,LICP+1),NCP,.false.)
        call Apply_col(CPQ,NUPSG,NDWNC,CI(IOC+1),NDWNSG,SGM(ISGSTA+1),NCP,EXS%ICOUP(1,LICP+1),swap=.false.)

      end do
    end do

  else if (SGS%MidLev < IQ) then

    ! EXCITING CASE, MIDLEV<IQ<IP
    INDEO = 2*SGS%nLev+(IP*(IP-1))/2+IQ
    do MVSGM=1,CIS%nMidV
      do ISYUSG=1,SGS%nSym
        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)  ; if (NS1 == 0) cycle
        ISYUC = Mul(ISYPQ,ISYUSG)
        NUPC = CIS%NOW(1,ISYUC,MVSGM)         ; if (NUPC == 0) cycle
        ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
        NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
        ISYDSG = Mul(ISYUSG,ISYSGM)
        NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
        IOC = CIS%IOCSF(ISYUC,MVSGM,ISYCI)
        NCP = EXS%NOCP(INDEO,ISYUC,MVSGM)     ; if (NCP == 0) cycle
        LICP = EXS%IOCP(INDEO,ISYUC,MVSGM)
        ! CASE IS: UPPER HALF, EXCITE:
        call apply_row(CPQ,NDWNSG,NUPC,CI(IOC+1),NUPSG,SGM(ISGSTA+1),NCP,EXS%ICOUP(1,LICP+1),swap=.false.)
      end do
    end do

  else

    ! EXCITING CASE, IQ<=MIDLEV<IP
    do MVSGM=1,CIS%nMidV
      do MV = 1, 2
        MVX = EXS%MVL(MVSGM,MV) ; if (MVX == 0) cycle
        do ISYUSG=1,SGS%nSym
          NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM) ; if (NS1 == 0) cycle
          ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
          ISYDSG = Mul(ISYUSG,ISYSGM)
          ISYUC = Mul(ISYP,ISYUSG)
          ISYDC = Mul(ISYQ,ISYDSG)

          NUPC = CIS%NOW(1,ISYUC,MVX)         ; if (NUPC == 0) cycle
          NDWNC = CIS%NOW(2,ISYDC,MVX)        ; if (NDWNC == 0) cycle

          INDEO = merge(IP, SGS%nLev+IP, MV==1)

          NCP1 = EXS%NOCP(INDEO,ISYUC,MVX)  ; if (NCP1 == 0) cycle
          NTMP = NUPSG*NDWNC
          EXS%SGTMP(1:NTMP) = Zero
          LICP = EXS%IOCP(INDEO,ISYUC,MVX)
          IOC = CIS%IOCSF(ISYUC,MVX,ISYCI)

          INDEO = merge(IQ, SGS%nLev+IQ, MV==1)

          NCP2 = EXS%NOCP(INDEO,ISYDC,MVX) ; if (NCP2 == 0) cycle
          ! CASE IS: UPPER HALF, EXCITE:
          call Apply_row(CPQ,NDWNC,NUPC,CI(IOC+1),NUPSG,EXS%SGTMP,NCP1,EXS%ICOUP(1,LICP+1),swap=.false.)
          NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
          LICP = EXS%IOCP(INDEO,ISYDC,MVX)
          ! CASE IS: LOWER HALF, EXCITE:
          call Apply_col(One,NUPSG,NDWNC,EXS%SGTMP,NDWNSG,SGM(ISGSTA+1),NCP2,EXS%ICOUP(1,LICP+1),swap=.false.)

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
        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM) ; if (NS1 == 0) cycle
        ISYDSG = Mul(ISYUSG,ISYSGM)
        ISYDC = Mul(ISYPQ,ISYDSG)
        NDWNC = CIS%NOW(2,ISYDC,MVSGM)       ; if (NDWNC == 0) cycle
        ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
        NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
        IOC = CIS%IOCSF(ISYUSG,MVSGM,ISYCI)
        NCP = EXS%NOCP(INDEO,ISYDSG,MVSGM)   ; if (NCP == 0) cycle
        NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
        LICP = EXS%IOCP(INDEO,ISYDSG,MVSGM)
        ! CASE IS: LOWER HALF, DEEXCITE:
        call sort_icoup_block(EXS%ICOUP(1,LICP+1),NCP,.true.)
        call Apply_col(CPQ,NUPSG,NDWNC,CI(IOC+1),NDWNSG,SGM(ISGSTA+1),NCP,EXS%ICOUP(1,LICP+1),swap=.True.)
      end do
    end do

  else if (SGS%MidLev < IP) then

    ! DEEXCITING OPERATOR, MIDLEV<IP<IQ
    INDEO = 2*SGS%nLev+(IQ*(IQ-1))/2+IP
    do MVSGM=1,CIS%nMidV
      do ISYUSG=1,SGS%nSym
        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)  ; if (NS1 == 0) cycle
        ISYUC = Mul(ISYPQ,ISYUSG)
        NUPC = CIS%NOW(1,ISYUC,MVSGM)         ; if (NUPC == 0) cycle
        ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
        NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
        ISYDSG = Mul(ISYUSG,ISYSGM)
        NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
        IOC = CIS%IOCSF(ISYUC,MVSGM,ISYCI)
        NCP = EXS%NOCP(INDEO,ISYUSG,MVSGM)    ; if (NCP == 0) cycle
        LICP = EXS%IOCP(INDEO,ISYUSG,MVSGM)
        ! CASE IS: UPPER HALF, DEEXCITE:
        call Apply_row(CPQ,NDWNSG,NUPC,CI(IOC+1),NUPSG,SGM(ISGSTA+1),NCP,EXS%ICOUP(1,LICP+1),swap=.true.)
      end do
    end do

  else

    ! DEEXCITING CASE, IP<=MIDLEV<IQ.
    do MVSGM=1,CIS%nMidV
      do MV = 1, 2
         MVX = EXS%MVR(MVSGM,MV) ; if (MVX == 0) cycle

         do ISYUSG=1,SGS%nSym
           NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM) ; if (NS1 == 0) cycle
           ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
           NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
           ISYDSG = Mul(ISYUSG,ISYSGM)
           ISYUC = Mul(ISYQ,ISYUSG)
           ISYDC = Mul(ISYP,ISYDSG)

           NUPC = CIS%NOW(1,ISYUC,MVX)          ; if (NUPC == 0) cycle
           NDWNC = CIS%NOW(2,ISYDC,MVX)         ; if (NDWNC == 0) cycle

           INDEO = merge(IQ, SGS%nLev+IQ, MV==1)

           NCP1 = EXS%NOCP(INDEO,ISYUSG,MVSGM)  ; if (NCP1 == 0) cycle
           NTMP = NUPSG*NDWNC
           EXS%SGTMP(1:NTMP) = Zero
           LICP = EXS%IOCP(INDEO,ISYUSG,MVSGM)
           IOC = CIS%IOCSF(ISYUC,MVX,ISYCI)

           INDEO = merge(IP, SGS%nLev+IP, MV==1)

           NCP2 = EXS%NOCP(INDEO,ISYDSG,MVSGM) ; if (NCP2 == 0) cycle
           ! CASE IS: UPPER HALF, DEEXCITE:
           call Apply_row(CPQ,NDWNC,NUPC,CI(IOC+1),NUPSG,EXS%SGTMP,NCP1,EXS%ICOUP(1,LICP+1),swap=.true.)
           NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
           LICP = EXS%IOCP(INDEO,ISYDSG,MVSGM)
           ! CASE IS: LOWER HALF, DEEXCITE:
           call Apply_col(One,NUPSG,NDWNC,EXS%SGTMP,NDWNSG,SGM(ISGSTA+1),NCP2,EXS%ICOUP(1,LICP+1),swap=.true.)

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
        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM) ; if (NS1 == 0) cycle
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
          ICS = mod(IC/IPPOW,4) ; if (ICS == 0) cycle
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
        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM) ; if (NS1 == 0) cycle
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
          ICS = mod(JC/IPPOW,4) ; if (ICS == 0) cycle
          X = CPQ*real((1+ICS)/2,kind=wp)
          JSTA = ISGSTA + NUPSG*(J-1) + 1
          SGM(JSTA:JSTA+NUPSG-1) = SGM(JSTA:JSTA+NUPSG-1)+X*CI(JSTA:JSTA+NUPSG-1)
        end do
      end do
    end do

  end if
end if

contains



subroutine sort_icoup_block(ICOUP, NCP, swap)
  use Definitions, only: iwp
  implicit none
  integer(kind=iwp), intent(inout) :: ICOUP(3,NCP)
  integer(kind=iwp), intent(in) :: NCP
  logical(kind=iwp), intent(in) :: swap
  integer(kind=iwp) :: i, j
  integer(kind=iwp) :: temp(3)

  do i = 2, NCP
     temp = ICOUP(:,i)
     j = i-1
     do while (j >= 1)
        if (swap) then
           if (ICOUP(1,j) <= temp(1)) exit
        else
           if (ICOUP(2,j) <= temp(2)) exit
        end if
        ICOUP(:,j+1) = ICOUP(:,j)
        j = j - 1
     end do
     ICOUP(:,j+1) = temp
  end do
end subroutine sort_icoup_block

subroutine apply_col(CPQ, NUP, NDWNC, CI, NDWNSG, SIGMA, NCP, ICOUP, swap)
  use Definitions, only: wp, iwp
  implicit none

  integer(kind=iwp), intent(in) :: NUP, NDWNC, NDWNSG, NCP, ICOUP(3,NCP)
  real(kind=wp), intent(in) :: CPQ, CI(NUP,NDWNC)
  real(kind=wp), intent(inout) :: SIGMA(NUP,NDWNSG)
  logical(kind=iwp), intent(in) :: swap

  integer(kind=iwp) :: ICP, start, finish, k, nk, i, I2, ICP2
  integer(kind=iwp) :: I1
  ! Tune these parameters, KBLOCK=15,32,..., THR_GEMV=6,8,12,
  integer(kind=iwp), parameter :: KBLOCK=16, THR_GEMV=8

  integer(kind=iwp) :: I1list(KBLOCK)
  real(kind=wp) :: Xlist(KBLOCK)
  real(kind=wp) :: ABLOCK(NUP,KBLOCK)
  real(kind=wp) :: TEMP(NUP,1)

  real(kind=wp), pointer :: VTAB(:)
  real(kind=wp) :: X

  VTAB => EXS%VTab

  if (swap) then
    ICP = 1
    do while (ICP <= NCP)

      I2 = ICOUP(1,ICP)
      start = ICP
      finish = ICP

      do while (finish <= NCP .and. ICOUP(1,finish) == I2)
        finish = finish + 1
      end do

      do k = start, finish-1, KBLOCK

        nk = min(KBLOCK, finish-k)

        if (nk < THR_GEMV) then
           ! fallback scalar kernel
           do i=1,nk
             I1 = ICOUP(2,k+i-1)
             X = CPQ * VTAB(ICOUP(3,k+i-1))
!$OMP SIMD
             do ICP2=1,NUP
               SIGMA(ICP2,I2)=SIGMA(ICP2,I2)+X*CI(ICP2,I1)
             end do
           end do
        else
           do i = 1, nk
             I1list(i) = ICOUP(2, k+i-1)
             Xlist(i)  = CPQ * VTAB(ICOUP(3, k+i-1))
           end do

           do i = 1, nk
             ABLOCK(:,i) = CI(:, I1list(i))
           end do

           call DGEMM_('N','N', NUP, 1, nk, 1.0_wp, ABLOCK, NUP, Xlist, nk, 0.0_wp, TEMP, NUP)

           SIGMA(:,I2) = SIGMA(:,I2) + TEMP(:,1)
        end if

      end do

      ICP = finish
    end do

  else
    ICP = 1
    do while (ICP <= NCP)

      I2 = ICOUP(2,ICP)
      start = ICP
      finish = ICP

      do while (finish <= NCP .and. ICOUP(2,finish) == I2)
        finish = finish + 1
      end do

      do k = start, finish-1, KBLOCK

        nk = min(KBLOCK, finish-k)

        if (nk < THR_GEMV) then
           do i=1,nk
             I1 = ICOUP(1,k+i-1)
             X = CPQ * VTAB(ICOUP(3,k+i-1))
!$OMP SIMD
             do ICP2=1,NUP
               SIGMA(ICP2,I2)=SIGMA(ICP2,I2)+X*CI(ICP2,I1)
             end do
           end do
        else
           do i = 1, nk
             I1list(i) = ICOUP(1, k+i-1)
             Xlist(i)  = CPQ * VTAB(ICOUP(3, k+i-1))
           end do

           do i = 1, nk
             ABLOCK(:,i) = CI(:, I1list(i))
           end do

           call DGEMM_('N','N', NUP, 1, nk, 1.0_wp, ABLOCK, NUP, Xlist, nk, 0.0_wp, TEMP, NUP)

           SIGMA(:,I2) = SIGMA(:,I2) + TEMP(:,1)
        end if

      end do

      ICP = finish
    end do

  end if

  VTAB => null()
end subroutine apply_col

subroutine apply_row(CPQ, NDWN, NUPC, CI, NUPSG, SIGMA, NCP, ICOUP, swap)
  integer(kind=iwp), intent(in) :: NDWN, NUPC, NUPSG, NCP, ICOUP(3,NCP)
  real(kind=wp), intent(in) :: CPQ, CI(NUPC,NDWN)
  real(kind=wp), intent(inout) :: SIGMA(NUPSG,NDWN)
  logical(kind=iwp), intent(in) :: swap

  integer(kind=iwp) :: ICP, IDWN
  real(kind=wp), pointer :: VTAB(:), XLIST(:)
  integer(kind=iwp), pointer ::  I1LIST(:), I2LIST(:)

  VTAB => EXS%VTab
  XLIST => EXS%XLIST
  I1LIST => EXS%I1LIST
  I2LIST => EXS%I2LIST

  if (swap) then
    do ICP=1,NCP
      I1list(ICP)=ICOUP(2,ICP)
      I2list(ICP)=ICOUP(1,ICP)
      Xlist(ICP)=CPQ*VTAB(ICOUP(3,ICP))
    end do
  else
    do ICP=1,NCP
      I1list(ICP)=ICOUP(1,ICP)
      I2list(ICP)=ICOUP(2,ICP)
      Xlist(ICP)=CPQ*VTAB(ICOUP(3,ICP))
    end do
  end if

  do ICP=1,NCP
!$OMP SIMD
    do IDWN=1,NDWN
      SIGMA(I2list(ICP),IDWN)=SIGMA(I2list(ICP),IDWN)+Xlist(ICP)*CI(I1list(ICP),IDWN)
    end do
  end do

  VTAB => null()
end subroutine apply_row

end subroutine SG_Epq_Psi
