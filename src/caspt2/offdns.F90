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
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine OFFDNS(ISYM1,ICASE1,ISYM2,ICASE2,X1,nX1,X2,nX2,DPT2,nDPT2,Y,nY,LIST,MLIST)
! Compute off-diagonal contributions to a trans density matrix.
! Sub-diagonal blocks only. If a density matrix is required,
! i.e. both wave functions equal, use symmetry. Else, call
! again with wave functions interchanged.

use Symmetry_Info, only: Mul
use EQSOLV, only: IFCOUP, LLIST, NLIST
use Sigma_data, only: IFTEST, INCF1, INCF2, INCX1, INCX2, INCX3, INCY1, INCY2, INCY3, LEN1, LEN2, NLST1, NLST2, VAL1, VAL2
use caspt2_module, only: NAGEB, NAGTB, NASH, NASUP, NIGEJ, NIGTJ, NISH, NISUP, NORB, NSSH, NSYM, NTGEU, NTGTU, NTUV
use Constants, only: One, Two, Three, Six, Half, OneHalf
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ISYM1, ICASE1, ISYM2, ICASE2, nX1, nX2, nY, nDPT2, mLIST, LIST(mList)
real(kind=wp), intent(inout) :: X1(nX1), X2(nX2), DPT2(nDPT2), Y(nY)
integer(kind=iwp) :: ICD, ICEM, ICEP, ICGM, ICGP, IDIA, IDIT, IDJB, IDOFF, IDTA, IJSYM, IMLTOP, IOFCD(8,8), IOFCEM(8,8), &
                     IOFCEP(8,8), IOFCGM(8,8), IOFCGP(8,8), IOFDIA(8), IOFDIT(8), IOFDTA(8), IOXIA, ISYM, ISYM12, ISYMA, ISYMAB, &
                     ISYMB, ISYMI, ISYMIJ, ISYMJ, IX, IXIA, IXTA, IXTI, IY, JSYM, KOD, LLST1, LLST2, NA, NAS1, NAS2, NB, NI, NIS1, &
                     NJ, NO, NT, NU
real(kind=wp), parameter :: SQR2 = sqrt(Two), SQR3 = sqrt(Three), SQR32 = sqrt(OneHalf), SQRI2 = sqrt(Half), SQRI6 = One/sqrt(Six)

NA = 0 ! dummy initialize

IFTEST = 0
IMLTOP = 2
KOD = IFCOUP(ICASE2,ICASE1)
if (KOD == 0) return

ISYM12 = Mul(ISYM1,ISYM2)
NAS1 = NASUP(ISYM1,ICASE1)
NIS1 = NISUP(ISYM1,ICASE1)
NAS2 = NASUP(ISYM2,ICASE2)
! Set up various offset arrays:
IDOFF = 0
do ISYM=1,NSYM
  NI = NISH(ISYM)
  NA = NASH(ISYM)
  NO = NORB(ISYM)
  IOFDIT(ISYM) = IDOFF+NO*NI
  IOFDIA(ISYM) = IOFDIT(ISYM)+NO*NA
  IOFDTA(ISYM) = IOFDIA(ISYM)+NI
  IDOFF = IDOFF+NO**2
  ICD = 0
  ICEP = 0
  ICEM = 0
  ICGP = 0
  ICGM = 0
  do JSYM=1,NSYM
    IJSYM = Mul(ISYM,JSYM)
    IOFCD(ISYM,JSYM) = ICD
    IOFCEP(ISYM,JSYM) = ICEP
    IOFCEM(ISYM,JSYM) = ICEM
    IOFCGP(ISYM,JSYM) = ICGP
    IOFCGM(ISYM,JSYM) = ICGM
    ICD = ICD+NSSH(JSYM)*NISH(IJSYM)
    ICEP = ICEP+NSSH(JSYM)*NIGEJ(IJSYM)
    ICEM = ICEM+NSSH(JSYM)*NIGTJ(IJSYM)
    ICGP = ICGP+NISH(JSYM)*NAGEB(IJSYM)
    ICGM = ICGM+NISH(JSYM)*NAGTB(IJSYM)
  end do
end do

select case (KOD)
  ! -----------------------------------------------
  case (1)
    ! A&BP One-el
    NLST1 = NLIST(ISYM1,ISYM2,12)
    NLST2 = NLIST(ISYM1,ISYM2,14)
    if (NLST1*NLST2 /= 0) then
      LLST1 = LLIST(ISYM1,ISYM2,12)
      VAL1(1) = One
      VAL1(2) = Two
      LLST2 = LLIST(ISYM1,ISYM2,14)
      VAL2(1) = One
      VAL2(2) = SQR2
      IXTI = 1
      INCX1 = 1
      INCX2 = NASH(ISYM1)
      IDIT = 1+IOFDIT(ISYM12)
      INCF1 = NORB(ISYM12)
      INCF2 = 1
      IY = 1
      INCY1 = 1
      INCY2 = NTGEU(ISYM2)
      call MLTSCA(IMLTOP,LIST(LLST1),NLST1,LIST(LLST2),NLST2,X1(IXTI),size(X1(IXTI:)),DPT2(IDIT),size(DPT2(IDIT:)),Y(IY), &
                  size(Y(IY:)))
    end if

    ! A&BP Two-el
    NLST1 = NLIST(ISYM1,ISYM2,3)
    NLST2 = NLIST(ISYM1,ISYM2,14)
    if (NLST1*NLST2 /= 0) then
      LLST1 = LLIST(ISYM1,ISYM2,3)
      VAL1(1) = -One
      VAL1(2) = -Two
      LLST2 = LLIST(ISYM1,ISYM2,14)
      VAL2(1) = One
      VAL2(2) = SQR2
      IX = 1
      INCX1 = 1
      INCX2 = NTUV(ISYM1)
      IDIT = 1+IOFDIT(ISYM12)
      INCF1 = NORB(ISYM12)
      INCF2 = 1
      IY = 1
      INCY1 = 1
      INCY2 = NTGEU(ISYM2)
      call MLTSCA(IMLTOP,LIST(LLST1),NLST1,LIST(LLST2),NLST2,X2(IX),size(X2(IX:)),DPT2(IDIT),size(DPT2(IDIT:)),Y(IY),size(Y(IY:)))
    end if
    ! -----------------------------------------------
  case (2)
    ! A&BM One-el
    NLST1 = NLIST(ISYM1,ISYM2,13)
    NLST2 = NLIST(ISYM1,ISYM2,15)
    if (NLST1*NLST2 /= 0) then
      LLST1 = LLIST(ISYM1,ISYM2,13)
      VAL1(1) = Three
      VAL1(2) = -Three
      LLST2 = LLIST(ISYM1,ISYM2,15)
      ! Original:
      !VAL2(1) = -One
      !VAL2(2) = One
      ! Fix for sign error noted by Takeshi, May 2015:
      VAL2(1) = One
      VAL2(2) = -One
      IXTI = 1
      INCX1 = 1
      INCX2 = NASH(ISYM1)
      IDIT = 1+IOFDIT(ISYM12)
      INCF1 = NORB(ISYM12)
      INCF2 = 1
      IY = 1
      INCY1 = 1
      INCY2 = NTGTU(ISYM2)
      call MLTSCA(IMLTOP,LIST(LLST1),NLST1,LIST(LLST2),NLST2,X1(IXTI),size(X1(IXTI:)),DPT2(IDIT),size(DPT2(IDIT:)),Y(IY), &
                  size(Y(IY:)))
    end if

    ! A&BM Two-el
    NLST1 = NLIST(ISYM1,ISYM2,4)
    NLST2 = NLIST(ISYM1,ISYM2,15)
    if (NLST1*NLST2 /= 0) then
      LLST1 = LLIST(ISYM1,ISYM2,4)
      VAL1(1) = -One
      VAL1(2) = One
      LLST2 = LLIST(ISYM1,ISYM2,15)
      ! Original:
      !VAL2(1) = -One
      !VAL2(2) = One
      ! Fix for sign error noted by Takeshi, May 2015:
      VAL2(1) = One
      VAL2(2) = -One
      IX = 1
      INCX1 = 1
      INCX2 = NTUV(ISYM1)
      IDIT = 1+IOFDIT(ISYM12)
      INCF1 = NORB(ISYM12)
      INCF2 = 1
      IY = 1
      INCY1 = 1
      INCY2 = NTGTU(ISYM2)
      call MLTSCA(IMLTOP,LIST(LLST1),NLST1,LIST(LLST2),NLST2,X2(IX),size(X2(IX:)),DPT2(IDIT),size(DPT2(IDIT:)),Y(IY),size(Y(IY:)))
    end if
    ! -----------------------------------------------
  case (3)
    ! A&D  Two-el
    LLST1 = LLIST(ISYM1,ISYM2,1)
    NLST1 = NLIST(ISYM1,ISYM2,1)
    if (NLST1 /= 0) then
      VAL1(1) = One
      VAL1(2) = One
      IX = 1
      INCX1 = 1
      INCX2 = NAS1
      IDTA = 1+IOFDTA(ISYM12)
      INCF1 = 1
      INCF2 = NORB(ISYM12)
      IY = 1+NAS2*IOFCD(ISYM2,ISYM12)
      INCY1 = 1
      INCY2 = NAS2
      INCY3 = NAS2*NISH(ISYM1)
      LEN1 = NISH(ISYM1)
      LEN2 = NSSH(ISYM12)
      call MLTMV(IMLTOP,LIST(LLST1),NLST1,X2(IX),size(X2(IX:)),DPT2(IDTA),size(DPT2(IDTA:)),Y(IY),size(Y(IY:)))
    end if
    ! -----------------------------------------------
  case (4)
    ! A&EP One-el
    NT = NASH(ISYM1)
    if ((ISYM2 == ISYM1) .and. (NT /= 0)) then
      do ISYMIJ=1,NSYM
        ISYMA = Mul(ISYMIJ,ISYM1)
        NA = NSSH(ISYMA)
        NLST1 = NLIST(ISYM1,ISYMIJ,14)
        if (NA*NLST1 /= 0) then
          LLST1 = LLIST(ISYM1,ISYMIJ,14)
          VAL1(1) = One
          VAL1(2) = SQR2
          IXTI = 1
          INCX1 = NT
          INCX2 = 1
          IDIA = 1+IOFDIA(ISYMA)
          INCF1 = 1
          INCF2 = NORB(ISYMA)
          IY = 1+NAS2*IOFCEP(ISYM2,ISYMA)
          INCY1 = NT*NA
          INCY2 = 1
          INCY3 = NT
          LEN1 = NT
          LEN2 = NA
          call MLTMV(IMLTOP,LIST(LLST1),NLST1,X1(IXTI),size(X1(IXTI:)),DPT2(IDIA),size(DPT2(IDIA:)),Y(IY),size(Y(IY:)))
        end if
      end do
    end if
    ! -----------------------------------------------
  case (5)
    ! A&EM One-el
    NT = NASH(ISYM1)
    if ((ISYM2 == ISYM1) .and. (NT /= 0)) then
      do ISYMIJ=1,NSYM
        ISYMA = Mul(ISYMIJ,ISYM1)
        NA = NSSH(ISYMA)
        NLST1 = NLIST(ISYM1,ISYMIJ,15)
        if (NA*NLST1 /= 0) then
          LLST1 = LLIST(ISYM1,ISYMIJ,15)
          VAL1(1) = -SQR3
          VAL1(2) = SQR3
          IXTI = 1
          INCX1 = NT
          INCX2 = 1
          IDIA = 1+IOFDIA(ISYMA)
          INCF1 = 1
          INCF2 = NORB(ISYMA)
          IY = 1+NAS2*IOFCEM(ISYM2,ISYMA)
          INCY1 = NT*NA
          INCY2 = 1
          INCY3 = NT
          LEN1 = NT
          LEN2 = NA
          call MLTMV(IMLTOP,LIST(LLST1),NLST1,X1(IXTI),size(X1(IXTI:)),DPT2(IDIA),size(DPT2(IDIA:)),Y(IY),size(Y(IY:)))
        end if
      end do
    end if
    ! -----------------------------------------------
  case (6)
    ! BP&EP Two-el
    NA = NSSH(ISYM12)
    NLST1 = NLIST(ISYM1,ISYM2,9)
    if (NA*NLST1 /= 0) then
      LLST1 = LLIST(ISYM1,ISYM2,9)
      VAL1(1) = SQRI2
      VAL1(2) = SQRI2
      IX = 1
      INCX1 = 1
      INCX2 = NAS1
      IDTA = 1+IOFDTA(ISYM12)
      INCF1 = 1
      INCF2 = NORB(ISYM12)
      IY = 1+NAS2*IOFCEP(ISYM2,ISYM12)
      INCY1 = 1
      INCY2 = NAS2*NA
      INCY3 = NAS2
      LEN1 = NIS1
      LEN2 = NA
      call MLTMV(IMLTOP,LIST(LLST1),NLST1,X2(IX),size(X2(IX:)),DPT2(IDTA),size(DPT2(IDTA:)),Y(IY),size(Y(IY:)))
    end if
    ! -----------------------------------------------
  case (7)
    ! BM&EM Two-el
    NA = NSSH(ISYM12)
    NLST1 = NLIST(ISYM1,ISYM2,10)
    if (NA*NLST1 /= 0) then
      LLST1 = LLIST(ISYM1,ISYM2,10)
      ! Original:
      !VAL1(1) = -SQRI6
      !VAL1(2) = SQRI6
      ! Fix for sign error noted by Takeshi, May 2015:
      VAL1(1) = SQRI6
      VAL1(2) = -SQRI6
      IX = 1
      INCX1 = 1
      INCX2 = NAS1
      IDTA = 1+IOFDTA(ISYM12)
      INCF1 = 1
      INCF2 = NORB(ISYM12)
      IY = 1+NAS2*IOFCEM(ISYM2,ISYM12)
      INCY1 = 1
      INCY2 = NAS2*NA
      INCY3 = NAS2
      LEN1 = NIS1
      LEN2 = NA
      call MLTMV(IMLTOP,LIST(LLST1),NLST1,X2(IX),size(X2(IX:)),DPT2(IDTA),size(DPT2(IDTA:)),Y(IY),size(Y(IY:)))
    end if
    ! -----------------------------------------------
  case (8)
    ! C&D  One-el
    NI = NISH(ISYM12)
    NLST1 = NLIST(ISYM1,ISYM2,11)
    if (NI*NLST1 /= 0) then
      LLST1 = LLIST(ISYM1,ISYM2,11)
      VAL1(1) = Two
      VAL1(2) = One
      IXTA = 1
      INCX1 = 1
      INCX2 = NASH(ISYM1)
      IDIT = 1+IOFDIT(ISYM12)
      INCF1 = NORB(ISYM12)
      INCF2 = 1
      IY = 1+NAS2*IOFCD(ISYM2,ISYM1)
      INCY1 = 1
      INCY2 = NAS2*NI
      INCY3 = NAS2
      LEN1 = NSSH(ISYM1)
      LEN2 = NI
      call MLTMV(IMLTOP,LIST(LLST1),NLST1,X1(IXTA),size(X1(IXTA:)),DPT2(IDIT),size(DPT2(IDIT:)),Y(IY),size(Y(IY:)))
    end if

    ! C&D  Two-el
    NI = NISH(ISYM12)
    NLST1 = NLIST(ISYM1,ISYM2,2)
    if (NI*NLST1 /= 0) then
      LLST1 = LLIST(ISYM1,ISYM2,2)
      VAL1(1) = -One
      VAL1(2) = -One
      IX = 1
      INCX1 = 1
      INCX2 = NAS1
      IDIT = 1+IOFDIT(ISYM12)
      INCF1 = NORB(ISYM12)
      INCF2 = 1
      IY = 1+NAS2*IOFCD(ISYM2,ISYM1)
      INCY1 = 1
      INCY2 = NAS2*NI
      INCY3 = NAS2
      LEN1 = NSSH(ISYM1)
      LEN2 = NI
      call MLTMV(IMLTOP,LIST(LLST1),NLST1,X2(IX),size(X2(IX:)),DPT2(IDIT),size(DPT2(IDIT:)),Y(IY),size(Y(IY:)))
    end if
    ! -----------------------------------------------
  case (9)
    ! C&FP One-el
    NLST1 = NLIST(ISYM1,ISYM2,12)
    NLST2 = NLIST(ISYM1,ISYM2,16)
    if (NLST1 /= 0) then
      LLST1 = LLIST(ISYM1,ISYM2,12)
      VAL1(1) = -One
      VAL1(2) = -Two
      LLST2 = LLIST(ISYM1,ISYM2,16)
      VAL2(1) = One
      VAL2(2) = SQR2
      IXTA = 1
      INCX1 = 1
      INCX2 = NASH(ISYM1)
      IDTA = 1+IOFDTA(ISYM12)
      INCF1 = 1
      INCF2 = NORB(ISYM12)
      IY = 1
      INCY1 = 1
      INCY2 = NTGEU(ISYM2)
      call MLTSCA(IMLTOP,LIST(LLST1),NLST1,LIST(LLST2),NLST2,X1(IXTA),size(X1(IXTA:)),DPT2(IDTA),size(DPT2(IDTA:)),Y(IY), &
                  size(Y(IY:)))
    end if

    ! C&FP Two-el
    NLST1 = NLIST(ISYM1,ISYM2,5)
    NLST2 = NLIST(ISYM1,ISYM2,16)
    if (NLST1*NLST2 /= 0) then
      LLST1 = LLIST(ISYM1,ISYM2,5)
      VAL1(1) = One
      VAL1(2) = Two
      LLST2 = LLIST(ISYM1,ISYM2,16)
      VAL2(1) = One
      VAL2(2) = SQR2
      IX = 1
      INCX1 = 1
      INCX2 = NTUV(ISYM1)
      IDTA = 1+IOFDTA(ISYM12)
      INCF1 = 1
      INCF2 = NORB(ISYM12)
      IY = 1
      INCY1 = 1
      INCY2 = NTGEU(ISYM2)
      call MLTSCA(IMLTOP,LIST(LLST1),NLST1,LIST(LLST2),NLST2,X2(IX),size(X2(IX:)),DPT2(IDTA),size(DPT2(IDTA:)),Y(IY),size(Y(IY:)))
    end if
    ! -----------------------------------------------
  case (10)
    ! C&FM One-el
    NLST1 = NLIST(ISYM1,ISYM2,13)
    NLST2 = NLIST(ISYM1,ISYM2,17)
    if (NLST1*NLST2 /= 0) then
      LLST1 = LLIST(ISYM1,ISYM2,13)
      VAL1(1) = -One
      VAL1(2) = One
      LLST2 = LLIST(ISYM1,ISYM2,17)
      VAL2(1) = One
      VAL2(2) = -One
      IXTA = 1
      INCX1 = 1
      INCX2 = NASH(ISYM1)
      IDTA = 1+IOFDTA(ISYM12)
      INCF1 = 1
      INCF2 = NORB(ISYM12)
      IY = 1
      INCY1 = 1
      INCY2 = NTGTU(ISYM2)
      call MLTSCA(IMLTOP,LIST(LLST1),NLST1,LIST(LLST2),NLST2,X1(IXTA),size(X1(IXTA:)),DPT2(IDTA),size(DPT2(IDTA:)),Y(IY), &
                  size(Y(IY:)))
    end if

    ! C&FM Two-el
    NLST1 = NLIST(ISYM1,ISYM2,6)
    NLST2 = NLIST(ISYM1,ISYM2,17)
    if (NLST1*NLST2 /= 0) then
      LLST1 = LLIST(ISYM1,ISYM2,6)
      VAL1(1) = -One
      VAL1(2) = One
      LLST2 = LLIST(ISYM1,ISYM2,17)
      VAL2(1) = One
      VAL2(2) = -One
      IX = 1
      INCX1 = 1
      INCX2 = NTUV(ISYM1)
      IDTA = 1+IOFDTA(ISYM12)
      INCF1 = 1
      INCF2 = NORB(ISYM12)
      IY = 1
      INCY1 = 1
      INCY2 = NTGTU(ISYM2)
      call MLTSCA(IMLTOP,LIST(LLST1),NLST1,LIST(LLST2),NLST2,X2(IX),size(X2(IX:)),DPT2(IDTA),size(DPT2(IDTA:)),Y(IY),size(Y(IY:)))
    end if
    ! -----------------------------------------------
  case (11)
    ! C&GP One-el
    NT = NASH(ISYM1)
    if ((ISYM2 == ISYM1) .and. (NT /= 0)) then
      do ISYMAB=1,NSYM
        ISYMI = Mul(ISYMAB,ISYM1)
        NI = NISH(ISYMI)
        NLST1 = NLIST(ISYM1,ISYMAB,16)
        if (NI*NLST1 /= 0) then
          LLST1 = LLIST(ISYM1,ISYMAB,16)
          VAL1(1) = SQRI2
          VAL1(2) = One
          IXTA = 1
          INCX1 = NT
          INCX2 = 1
          IDIA = 1+IOFDIA(ISYMI)
          INCF1 = NORB(ISYMI)
          INCF2 = 1
          IY = 1+NAS2*IOFCGP(ISYM2,ISYMI)
          INCY1 = NT*NI
          INCY2 = 1
          INCY3 = NT
          LEN1 = NT
          LEN2 = NI
          call MLTMV(IMLTOP,LIST(LLST1),NLST1,X1(IXTA),size(X1(IXTA:)),DPT2(IDIA),size(DPT2(IDIA:)),Y(IY),size(Y(IY:)))
        end if
      end do
    end if
    ! -----------------------------------------------
  case (12)
    ! C&GM One-el
    NT = NASH(ISYM1)
    if ((ISYM2 == ISYM1) .and. (NT /= 0)) then
      do ISYMAB=1,NSYM
        ISYMI = Mul(ISYMAB,ISYM1)
        NI = NISH(ISYMI)
        NLST1 = NLIST(ISYM1,ISYMAB,17)
        if (NI*NLST1 /= 0) then
          LLST1 = LLIST(ISYM1,ISYMAB,17)
          VAL1(1) = SQR32
          VAL1(2) = -SQR32
          IXTA = 1
          INCX1 = NT
          INCX2 = 1
          IDIA = 1+IOFDIA(ISYMI)
          INCF1 = NORB(ISYMI)
          INCF2 = 1
          IY = 1+NAS2*IOFCGM(ISYM2,ISYMI)
          INCY1 = NT*NI
          INCY2 = 1
          INCY3 = NT
          LEN1 = NT
          LEN2 = NI
          call MLTMV(IMLTOP,LIST(LLST1),NLST1,X1(IXTA),size(X1(IXTA:)),DPT2(IDIA),size(DPT2(IDIA:)),Y(IY),size(Y(IY:)))
        end if
      end do
    end if
    ! -----------------------------------------------
  case (13)
    ! D&EP One-el
    NT = NASH(ISYM2)
    if ((ISYM1 == 1) .and. (NT /= 0)) then
      IOXIA = 0
      do ISYMI=1,NSYM
        NI = NISH(ISYMI)
        NA = NSSH(ISYMI)
        ISYMIJ = Mul(ISYMI,ISYM2)
        LLST1 = LLIST(ISYMI,ISYMIJ,14)
        NLST1 = NLIST(ISYMI,ISYMIJ,14)
        if (NI*NA*NLST1 /= 0) then
          VAL1(1) = SQRI2
          VAL1(2) = One
          IXIA = IOXIA+1
          INCX1 = 1
          INCX2 = NI
          IDIT = 1+IOFDIT(ISYM2)
          INCF1 = 1
          INCF2 = NORB(ISYM2)
          IY = 1+NAS2*IOFCEP(ISYM2,ISYMI)
          INCY1 = NT*NA
          INCY2 = NT
          INCY3 = 1
          LEN1 = NA
          LEN2 = NT
          call MLTMV(IMLTOP,LIST(LLST1),NLST1,X1(IXIA),size(X1(IXIA:)),DPT2(IDIT),size(DPT2(IDIT:)),Y(IY),size(Y(IY:)))
        end if
        IOXIA = IOXIA+NI*NA
      end do
    end if

    ! D&EP Two-el
    NLST1 = NLIST(ISYM1,ISYM2,7)
    NU = NASH(ISYM2)
    if (NLST1*NU /= 0) then
      LLST1 = LLIST(ISYM1,ISYM2,7)
      VAL1(1) = -One
      VAL1(2) = -One
      do ISYMA=1,NSYM
        NA = NSSH(ISYMA)
        ISYMI = Mul(ISYMA,ISYM1)
        NI = NISH(ISYMI)
        ISYMIJ = Mul(ISYMI,ISYM12)
        LLST2 = LLIST(ISYMI,ISYMIJ,14)
        NLST2 = NLIST(ISYMI,ISYMIJ,14)
        if (NA*NI*NLST2 /= 0) then
          VAL2(1) = SQRI2
          VAL2(2) = One
          IX = 1+NAS1*IOFCD(ISYM1,ISYMA)
          INCX1 = 1
          INCX2 = NAS1
          INCX3 = NAS1*NI
          IDIT = 1+IOFDIT(ISYM12)
          INCF1 = NORB(ISYM12)
          INCF2 = 1
          IY = 1+NU*IOFCEP(ISYM2,ISYMA)
          INCY1 = 1
          INCY2 = NU*NA
          INCY3 = NU
          LEN1 = NA
          call MLTDXP(IMLTOP,LIST(LLST1),LIST(LLST2),X2(IX:),size(X2(IX:)),DPT2(IDIT:),size(DPT2(IDIT:)),Y(IY:),size(Y(IY:)))
        end if
      end do
    end if
    ! -----------------------------------------------
  case (14)
    ! D&EM One-el
    NT = NASH(ISYM2)
    if ((ISYM1 == 1) .and. (NT /= 0)) then
      IOXIA = 0
      do ISYMI=1,NSYM
        NI = NISH(ISYMI)
        NA = NSSH(ISYMI)
        ISYMIJ = Mul(ISYMI,ISYM2)
        LLST1 = LLIST(ISYMI,ISYMIJ,15)
        NLST1 = NLIST(ISYMI,ISYMIJ,15)
        if (NI*NLST1 /= 0) then
          VAL1(1) = SQR32
          VAL1(2) = -SQR32
          IXIA = IOXIA+1
          INCX1 = 1
          INCX2 = NI
          IDIT = 1+IOFDIT(ISYM2)
          INCF1 = 1
          INCF2 = NORB(ISYM2)
          IY = 1+NAS2*IOFCEM(ISYM2,ISYMI)
          INCY1 = NT*NA
          INCY2 = NT
          INCY3 = 1
          LEN1 = NA
          LEN2 = NT
          call MLTMV(IMLTOP,LIST(LLST1),NLST1,X1(IXIA),size(X1(IXIA:)),DPT2(IDIT),size(DPT2(IDIT:)),Y(IY),size(Y(IY:)))
        end if
        IOXIA = IOXIA+NI*NA
      end do
    end if

    ! D&EM Two-el
    NLST1 = NLIST(ISYM1,ISYM2,7)
    NU = NASH(ISYM2)
    if (NLST1*NU /= 0) then
      LLST1 = LLIST(ISYM1,ISYM2,7)
      VAL1(1) = -One
      VAL1(2) = One
      do ISYMA=1,NSYM
        NA = NSSH(ISYMA)
        ISYMI = Mul(ISYMA,ISYM1)
        NI = NISH(ISYMI)
        ISYMIJ = Mul(ISYMI,ISYM12)
        NLST2 = NLIST(ISYMI,ISYMIJ,15)
        if (NA*NI*NLST2 /= 0) then
          LLST2 = LLIST(ISYMI,ISYMIJ,15)
          VAL2(1) = SQRI6
          VAL2(2) = -SQRI6
          IX = 1+NAS1*IOFCD(ISYM1,ISYMA)
          INCX1 = 1
          INCX2 = NAS1
          INCX3 = NAS1*NI
          IDIT = 1+IOFDIT(ISYM12)
          INCF1 = NORB(ISYM12)
          INCF2 = 1
          IY = 1+NU*IOFCEM(ISYM2,ISYMA)
          INCY1 = 1
          INCY2 = NU*NA
          INCY3 = NU
          LEN1 = NA
          call MLTDXP(IMLTOP,LIST(LLST1),LIST(LLST2),X2(IX:),size(X2(IX:)),DPT2(IDIT:),size(DPT2(IDIT:)),Y(IY:),size(Y(IY:)))
        end if
      end do
    end if
    ! -----------------------------------------------
  case (15)
    ! D&GP Two-el
    NLST1 = NLIST(ISYM1,ISYM2,8)
    NU = NASH(ISYM2)
    if (NLST1*NU /= 0) then
      LLST1 = LLIST(ISYM1,ISYM2,8)
      VAL1(1) = One
      VAL1(2) = One
      do ISYMA=1,NSYM
        ISYMI = Mul(ISYMA,ISYM1)
        NI = NISH(ISYMI)
        ISYMAB = Mul(ISYMA,ISYM12)
        NLST2 = NLIST(ISYMA,ISYMAB,16)
        if (NI*NLST2 /= 0) then
          LLST2 = LLIST(ISYMA,ISYMAB,16)
          VAL2(1) = SQRI2
          VAL2(2) = One
          IX = 1+NAS1*IOFCD(ISYM1,ISYMA)
          INCX1 = 1
          INCX2 = NAS1*NI
          INCX3 = NAS1
          IDTA = 1+IOFDTA(ISYM12)
          INCF1 = 1
          INCF2 = NORB(ISYM12)
          IY = 1+NU*IOFCGP(ISYM2,ISYMI)
          INCY1 = 1
          INCY2 = NU*NI
          INCY3 = NU
          LEN1 = NI
          call MLTDXP(IMLTOP,LIST(LLST1),LIST(LLST2),X2(IX:),size(X2(IX:)),DPT2(IDTA:),size(DPT2(IDTA:)),Y(IY:),size(Y(IY:)))
        end if
      end do
    end if
    ! -----------------------------------------------
  case (16)
    ! D&GM Two-el
    NLST1 = NLIST(ISYM1,ISYM2,8)
    NU = NASH(ISYM2)
    if (NLST1*NU /= 0) then
      LLST1 = LLIST(ISYM1,ISYM2,8)
      VAL1(1) = -One
      VAL1(2) = One

      do ISYMA=1,NSYM
        ISYMI = Mul(ISYMA,ISYM1)
        NI = NISH(ISYMI)
        ISYMAB = Mul(ISYMA,ISYM12)
        NLST2 = NLIST(ISYMA,ISYMAB,17)
        if (NI*NLST2 /= 0) then
          LLST2 = LLIST(ISYMA,ISYMAB,17)
          VAL2(1) = SQRI6
          VAL2(2) = -SQRI6
          IX = 1+NAS1*IOFCD(ISYM1,ISYMA)
          INCX1 = 1
          INCX2 = NAS1*NI
          INCX3 = NAS1
          IDTA = 1+IOFDTA(ISYM12)
          INCF1 = 1
          INCF2 = NORB(ISYM12)
          IY = 1+NU*IOFCGM(ISYM2,ISYMI)
          INCY1 = 1
          INCY2 = NU*NI
          INCY3 = NU
          LEN1 = NI
          call MLTDXP(IMLTOP,LIST(LLST1),LIST(LLST2),X2(IX:),size(X2(IX:)),DPT2(IDTA:),size(DPT2(IDTA:)),Y(IY:),size(Y(IY:)))
        end if
      end do
    end if
    ! -----------------------------------------------
  case (17)
    ! EP&HP Two-el
    NA = NSSH(ISYM12)
    NLST1 = NLIST(ISYM12,ISYM2,16)
    if (NA*NLST1 /= 0) then
      LLST1 = LLIST(ISYM12,ISYM2,16)
      VAL1(1) = SQRI2
      VAL1(2) = One
      IX = 1+NAS1*IOFCEP(ISYM1,ISYM12)
      INCX1 = NAS1
      INCX2 = 1
      INCX3 = NAS1*NA
      IDTA = 1+IOFDTA(ISYM1)
      INCF1 = NORB(ISYM1)
      INCF2 = 1
      IY = 1
      INCY1 = 1
      INCY2 = NAS2
      LEN1 = NAS1
      LEN2 = NIGEJ(ISYM2)
      call MLTR1(IMLTOP,LIST(LLST1),X2(IX:),size(X2(IX:)),DPT2(IDTA:),size(DPT2(IDTA:)),Y(IY:),size(Y(IY:)))
    end if
    ! -----------------------------------------------
  case (18)
    ! EM&HM Two-el
    NA = NSSH(ISYM12)
    NLST1 = NLIST(ISYM12,ISYM2,17)
    if (NA*NLST1 /= 0) then
      LLST1 = LLIST(ISYM12,ISYM2,17)
      VAL1(1) = SQRI2
      VAL1(2) = -SQRI2
      IX = 1+NAS1*IOFCEM(ISYM1,ISYM12)
      INCX1 = NAS1
      INCX2 = 1
      INCX3 = NAS1*NA
      IDTA = 1+IOFDTA(ISYM1)
      INCF1 = NORB(ISYM1)
      INCF2 = 1
      IY = 1
      INCY1 = 1
      INCY2 = NAS2
      LEN1 = NAS1
      LEN2 = NIGTJ(ISYM2)
      call MLTR1(IMLTOP,LIST(LLST1),X2(IX:),size(X2(IX:)),DPT2(IDTA:),size(DPT2(IDTA:)),Y(IY:),size(Y(IY:)))
    end if
    ! -----------------------------------------------
  case (19)
    ! FP&GP Two-el
    NI = NISH(ISYM12)
    NLST1 = NLIST(ISYM1,ISYM2,9)
    if (NI*NLST1 /= 0) then
      LLST1 = LLIST(ISYM1,ISYM2,9)
      VAL1(1) = -SQRI2
      VAL1(2) = -SQRI2
      IX = 1
      INCX1 = 1
      INCX2 = NAS1
      IDIT = 1+IOFDIT(ISYM12)
      INCF1 = NORB(ISYM12)
      INCF2 = 1
      IY = 1+NAS2*IOFCGP(ISYM2,ISYM12)
      INCY1 = 1
      INCY2 = NAS2*NI
      INCY3 = NAS2
      LEN1 = NIS1
      LEN2 = NI
      call MLTMV(IMLTOP,LIST(LLST1),NLST1,X2(IX),size(X2(IX:)),DPT2(IDIT),size(DPT2(IDIT:)),Y(IY),size(Y(IY:)))
    end if
    ! -----------------------------------------------
  case (20)
    ! FM&GM Two-el
    NI = NISH(ISYM12)
    NLST1 = NLIST(ISYM1,ISYM2,10)
    if (NI*NLST1 /= 0) then
      LLST1 = LLIST(ISYM1,ISYM2,10)
      VAL1(1) = -SQRI6
      VAL1(2) = SQRI6
      IX = 1
      INCX1 = 1
      INCX2 = NAS1
      IDIT = 1+IOFDIT(ISYM12)
      INCF1 = NORB(ISYM12)
      INCF2 = 1
      IY = 1+NAS2*IOFCGM(ISYM2,ISYM12)
      INCY1 = 1
      INCY2 = NAS2*NI
      INCY3 = NAS2
      LEN1 = NIS1
      LEN2 = NI
      call MLTMV(IMLTOP,LIST(LLST1),NLST1,X2(IX),size(X2(IX:)),DPT2(IDIT),size(DPT2(IDIT:)),Y(IY),size(Y(IY:)))
    end if
    ! -----------------------------------------------
  case (21)
    ! GP&HP Two-el
    LLST1 = LLIST(ISYM12,ISYM2,14)
    NLST1 = NLIST(ISYM12,ISYM2,14)
    if (NLST1 /= 0) then
      VAL1(1) = -SQRI2
      VAL1(2) = -One
      IX = 1+NAS1*IOFCGP(ISYM1,ISYM12)
      INCX1 = NAS1
      INCX2 = 1
      INCX3 = NAS1*NISH(ISYM12)
      IDIT = 1+IOFDIT(ISYM1)
      INCF1 = 1
      INCF2 = NORB(ISYM1)
      IY = 1
      INCY1 = NAS2
      INCY2 = 1
      LEN1 = NAS1
      LEN2 = NAS2
      call MLTR1(IMLTOP,LIST(LLST1),X2(IX:),size(X2(IX:)),DPT2(IDIT:),size(DPT2(IDIT:)),Y(IY:),size(Y(IY:)))
    end if
    ! -----------------------------------------------
  case (22)
    ! GM&HM Two-el
    LLST1 = LLIST(ISYM12,ISYM2,15)
    NLST1 = NLIST(ISYM12,ISYM2,15)
    if (NLST1 /= 0) then
      VAL1(1) = SQRI2
      VAL1(2) = -SQRI2
      IX = 1+NAS1*IOFCGM(ISYM1,ISYM12)
      INCX1 = NAS1
      INCX2 = 1
      INCX3 = NAS1*NISH(ISYM12)
      IDIT = 1+IOFDIT(ISYM1)
      INCF1 = 1
      INCF2 = NORB(ISYM1)
      IY = 1
      INCY1 = NAS2
      INCY2 = 1
      LEN1 = NAS1
      LEN2 = NAS2
      call MLTR1(IMLTOP,LIST(LLST1),X2(IX:),size(X2(IX:)),DPT2(IDIT:),size(DPT2(IDIT:)),Y(IY:),size(Y(IY:)))
    end if
    ! -----------------------------------------------
  case (23)
    ! D&HP One-el
    if (ISYM1 == 1) then
      IOXIA = 0
      do ISYMI=1,NSYM
        NI = NISH(ISYMI)
        ISYMA = ISYMI
        NA = NSSH(ISYMA)
        ISYMJ = Mul(ISYMI,ISYM2)
        NJ = NISH(ISYMJ)
        ISYMB = ISYMJ
        NB = NSSH(ISYMB)
        NLST1 = NLIST(ISYMI,ISYM2,14)
        NLST2 = NLIST(ISYMA,ISYM2,16)
        if (NI*NA*NJ*NB*NLST1*NLST2 /= 0) then
          LLST1 = LLIST(ISYMI,ISYM2,14)
          VAL1(1) = SQRI2
          VAL1(2) = One
          LLST2 = LLIST(ISYMA,ISYM2,16)
          VAL2(1) = SQRI2
          VAL2(2) = One
          IXIA = IOXIA+1
          INCX1 = 1
          INCX2 = NI
          IDJB = 1+IOFDIA(ISYMB)
          INCF1 = 1
          INCF2 = NORB(ISYMJ)
          IY = 1
          INCY1 = NAGEB(ISYM2)
          INCY2 = 1
          call MLTSCA(IMLTOP,LIST(LLST1),NLST1,LIST(LLST2),NLST2,X1(IXIA),size(X1(IXIA:)),DPT2(IDJB),size(DPT2(IDJB:)),Y(IY), &
                      size(Y(IY:)))
        end if
        IOXIA = IOXIA+NI*NA
      end do
    end if
    ! ---------------------------
  case (24)
    ! D&HM One-el
    if (ISYM1 == 1) then
      IOXIA = 0
      do ISYMI=1,NSYM
        NI = NISH(ISYMI)
        ISYMA = ISYMI
        NA = NSSH(ISYMA)
        ISYMJ = Mul(ISYMI,ISYM2)
        NJ = NISH(ISYMJ)
        ISYMB = ISYMJ
        NB = NSSH(ISYMB)
        NLST1 = NLIST(ISYMI,ISYM2,15)
        NLST2 = NLIST(ISYMA,ISYM2,17)
        if (NI*NA*NJ*NB*NLST1*NLST2 /= 0) then
          LLST1 = LLIST(ISYMI,ISYM2,15)
          VAL1(1) = SQR3*Half
          VAL1(2) = -VAL1(1)
          LLST2 = LLIST(ISYMA,ISYM2,17)
          VAL2(1) = One
          VAL2(2) = -One
          IXIA = IOXIA+1
          INCX1 = 1
          INCX2 = NI
          IDJB = 1+IOFDIA(ISYMB)
          INCF1 = 1
          INCF2 = NORB(ISYMJ)
          IY = 1
          INCY1 = NAGTB(ISYM2)
          INCY2 = 1
          call MLTSCA(IMLTOP,LIST(LLST1),NLST1,LIST(LLST2),NLST2,X1(IXIA),size(X1(IXIA:)),DPT2(IDJB),size(DPT2(IDJB:)),Y(IY), &
                      size(Y(IY:)))
        end if
        IOXIA = IOXIA+NI*NA
      end do
    end if
    ! ---------------------------
  case default
    write(u6,*) ' INTERNAL ERROR: OffDns reached invalid KOD=',KOD
    call Abend()
end select

end subroutine OFFDNS
