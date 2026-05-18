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
! SVC-20120305: This is a new parallel version of the SGM subroutine.
! Instead of the full array Y, distributed array indices lg_Y is passed,
! which refer to distributed arrays.

! Currently, only case H will be passed as a distributed array, since
! this is the largest array (about a factor of NI larger than case G)
! and it avoids communications of this array: if it is updated, we use
! only the local chunk. If is is used to update another case, that case
! is partially updated and then communicated.  The special routines are
! called MLTR1_GH, MLTR1_EH and MLTSCA_DH.

! The distributed arrays are (currently) layed out as 2-dimensional
! distributed arrays, with NAS rows and NIS columns (see subroutine
! RHS_ALLO in file par_rhs). The chunks are along the column indices,
! so each chunk has all the row indices (full columns).

subroutine SGM(IMLTOP,ISYM1,ICASE1,ISYM2,ICASE2,X1,nX1,X2,nX2,lg_Y,LIST,mList)

use Symmetry_Info, only: Mul
use Fockof, only: FAI, FAT, FIA, FIT, FTA, FTI, IOFFIA
use EQSOLV, only: IfCoup, LList, nList
use Sigma_data, only: INCF1, INCF2, INCX1, INCX2, INCX3, INCY1, INCY2, INCY3, LEN1, LEN2, nLst1, nLst2, Val1, Val2
use fake_GA, only: GA_Arrays
use caspt2_module, only: nAGEB, nAGTB, nAsh, nASUP, nIGEJ, nIGTJ, nIsh, nISUP, nSSh, nSym, nTGEU, nTGEU, nTGTU, nTUV
#ifdef _DEBUGPRINT_
use caspt2_module, only: Cases
#endif
use Constants, only: One, Two, Three, Six, Half, OneHalf
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IMLTOP, ISYM1, ICASE1, ISYM2, ICASE2, nX1, nX2, mList, LIST(mLIST)
real(kind=wp), intent(inout) :: X1(nX1), X2(nX2)
integer(kind=iwp) :: ICD, ICEM, ICEP, ICGM, ICGP, IJSYM, IOFCD(8,8), IOFCEM(8,8), IOFCEP(8,8), IOFCGM(8,8), IOFCGP(8,8), ISYM, &
                     ISYM12, ISYMA, ISYMAB, ISYMB, ISYMI, ISYMIJ, ISYMJ, IX, IXIA, IXTA, IXTI, IY, JSYM, JXOFF, KOD, lg_Y, LLST1, &
                     LLST2, NA, NAS1, NAS2, NB, NFA, NFI, NFT, NI, NIS1, NIS2, NJ, NT, NU
real(kind=wp), parameter :: SQR2 = sqrt(Two), SQR3 = sqrt(Three), SQR32 = sqrt(OneHalf), SQRI2 = sqrt(Half), SQRI6 = One/sqrt(Six)

! Compute a contribution from a single block (ISYM2,ICASE2)
! of expansion vector to a single block (ISYM1,ICASE1)
! of a sigma vector, if IMLTOP=0. If IMLTOP=1, the role of the
! two blocks are reversed. It is assumed that ICASE2 is
! greater than ICASE1, hence the reverse option.

!FUE if (((ICASE1 == 1) .and. (ICASE2 == 2)) .and. ((ISYM1 == 1) .and. (ISYM2 == 1))) IFTEST = 1
!PAM if ((ICASE1 == 5) .and. (ICASE2 > 11)) IFTEST = 1
#ifdef _DEBUGPRINT_
write(u6,'(A,10I5)') ' ENTERING SGM.'
write(u6,'(A,10I5)') '       IMLTOP:',IMLTOP
write(u6,'(A,10I5)') ' ISYM1,ICASE1:',ISYM1,ICASE1
write(u6,'(A,10I5)') ' ISYM2,ICASE2:',ISYM2,ICASE2
write(u6,'(A,A,A)') ' CASES: ',CASES(ICASE1),CASES(ICASE2)
#endif

if ((IMLTOP < 0) .or. (IMLTOP > 2)) then
  write(u6,*) 'Error in SGM: IMLTOP = ',IMLTOP
  call AbEnd()
end if

! SVC: IFCOUP is set in SIGMA_CASPT2
KOD = IFCOUP(ICASE2,ICASE1)
if (KOD == 0) return

ISYM12 = Mul(ISYM1,ISYM2)
NAS1 = NASUP(ISYM1,ICASE1)
NIS1 = NISUP(ISYM1,ICASE1)
NAS2 = NASUP(ISYM2,ICASE2)
NIS2 = NISUP(ISYM2,ICASE2)
do ISYM=1,NSYM
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

! SVC: this is an extra check, since coupling cases that are 0 should
! not have entered the sgm subroutine

select case (KOD)
  ! -----------------------------------------------
  case (1)
    !ICASE1 = 1
    !ICASE2 = 2

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
      INCF1 = NISH(ISYM12)
      INCF2 = 1
      IY = 1
      INCY1 = 1
      INCY2 = NTGEU(ISYM2)
      call MLTSCA(IMLTOP,LIST(LLST1),NLST1,LIST(LLST2),NLST2,X1(IXTI),size(X1(IXTI:)),FIT(ISYM12)%A,size(FIT(ISYM12)%A), &
                  GA_Arrays(lg_Y)%A(IY),size(GA_Arrays(lg_Y)%A(IY:)))
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
      INCF1 = NISH(ISYM12)
      INCF2 = 1
      IY = 1
      INCY1 = 1
      INCY2 = NTGEU(ISYM2)
      call MLTSCA(IMLTOP,LIST(LLST1),NLST1,LIST(LLST2),NLST2,X2(IX),size(X2(IX:)),FIT(ISYM12)%A,size(FIT(ISYM12)%A), &
                  GA_Arrays(lg_Y)%A(IY),size(GA_Arrays(lg_Y)%A(IY:)))
    end if
    ! -----------------------------------------------
  case (2)
    !ICASE1 = 1
    !ICASE2 = 3

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
      INCF1 = NISH(ISYM12)
      INCF2 = 1
      IY = 1
      INCY1 = 1
      INCY2 = NTGTU(ISYM2)
      call MLTSCA(IMLTOP,LIST(LLST1),NLST1,LIST(LLST2),NLST2,X1(IXTI),size(X1(IXTI:)),FIT(ISYM12)%A,size(FIT(ISYM12)%A), &
                  GA_Arrays(lg_Y)%A(IY),size(GA_Arrays(lg_Y)%A(IY:)))
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
      INCF1 = NISH(ISYM12)
      INCF2 = 1
      IY = 1
      INCY1 = 1
      INCY2 = NTGTU(ISYM2)
      call MLTSCA(IMLTOP,LIST(LLST1),NLST1,LIST(LLST2),NLST2,X2(IX),size(X2(IX:)),FIT(ISYM12)%A,size(FIT(ISYM12)%A), &
                  GA_Arrays(lg_Y)%A(IY),size(GA_Arrays(lg_Y)%A(IY:)))
    end if
    ! -----------------------------------------------
  case (3)
    !ICASE1 = 1
    !ICASE2 = 5

    ! A&D  Two-el
    NLST1 = NLIST(ISYM1,ISYM2,1)
    if (NLST1 /= 0) then
      LLST1 = LLIST(ISYM1,ISYM2,1)
      VAL1(1) = One
      VAL1(2) = One
      IX = 1
      INCX1 = 1
      INCX2 = NAS1
      INCF1 = NSSH(ISYM12)
      INCF2 = 1
      IY = 1+NAS2*IOFCD(ISYM2,ISYM12)
      INCY1 = 1
      INCY2 = NAS2
      INCY3 = NAS2*NISH(ISYM1)
      LEN1 = NISH(ISYM1)
      LEN2 = NSSH(ISYM12)
      call MLTMV(IMLTOP,LIST(LLST1),NLST1,X2(IX),size(X2(IX:)),FAT(ISYM12)%A,size(FAT(ISYM12)%A),GA_Arrays(lg_Y)%A(IY), &
                 size(GA_Arrays(lg_Y)%A(IY:)))
    end if
    ! -----------------------------------------------
  case (4)
    !ICASE1 = 1
    !ICASE2 = 6

    ! A&EP One-el
    if (ISYM2 == ISYM1) then
      NT = NASH(ISYM1)
      if (NT /= 0) then
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
            INCF1 = NSSH(ISYMA)
            INCF2 = 1
            IY = 1+NAS2*IOFCEP(ISYM2,ISYMA)
            INCY1 = NT*NA
            INCY2 = 1
            INCY3 = NT
            LEN1 = NT
            LEN2 = NA
            call MLTMV(IMLTOP,LIST(LLST1),NLST1,X1(IXTI),size(X1(IXTI:)),FAI(ISYMA)%A,size(FAI(ISYMA)%A),GA_Arrays(lg_Y)%A(IY), &
                       size(GA_Arrays(lg_Y)%A(IY:)))
          end if
        end do
      end if
    end if
    ! -----------------------------------------------
  case (5)
    !ICASE1 = 1
    !ICASE2 = 7

    ! A&EM One-el
    if (ISYM2 == ISYM1) then
      NT = NASH(ISYM1)
      if (NT /= 0) then
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
            INCF1 = NSSH(ISYMA)
            INCF2 = 1
            IY = 1+NAS2*IOFCEM(ISYM2,ISYMA)
            INCY1 = NT*NA
            INCY2 = 1
            INCY3 = NT
            LEN1 = NT
            LEN2 = NA
            call MLTMV(IMLTOP,LIST(LLST1),NLST1,X1(IXTI),size(X1(IXTI:)),FAI(ISYMA)%A,size(FAI(ISYMA)%A),GA_Arrays(lg_Y)%A(IY), &
                       size(GA_Arrays(lg_Y)%A(IY:)))
          end if
        end do
      end if
    end if
    ! -----------------------------------------------
  case (6)
    !ICASE1 = 2
    !ICASE2 = 6

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
      INCF1 = NSSH(ISYM12)
      INCF2 = 1
      IY = 1+NAS2*IOFCEP(ISYM2,ISYM12)
      INCY1 = 1
      INCY2 = NAS2*NA
      INCY3 = NAS2
      LEN1 = NIS1
      LEN2 = NA
      call MLTMV(IMLTOP,LIST(LLST1),NLST1,X2(IX),size(X2(IX:)),FAT(ISYM12)%A,size(FAT(ISYM12)%A),GA_Arrays(lg_Y)%A(IY), &
                 size(GA_Arrays(lg_Y)%A(IY:)))
    end if
    ! -----------------------------------------------
  case (7)
    !ICASE1 = 3
    !ICASE2 = 7

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
      INCF1 = NSSH(ISYM12)
      INCF2 = 1
      IY = 1+NAS2*IOFCEM(ISYM2,ISYM12)
      INCY1 = 1
      INCY2 = NAS2*NA
      INCY3 = NAS2
      LEN1 = NIS1
      LEN2 = NA
      call MLTMV(IMLTOP,LIST(LLST1),NLST1,X2(IX),size(X2(IX:)),FAT(ISYM12)%A,size(FAT(ISYM12)%A),GA_Arrays(lg_Y)%A(IY), &
                 size(GA_Arrays(lg_Y)%A(IY:)))
    end if
    ! -----------------------------------------------
  case (8)
    !ICASE1 = 4
    !ICASE2 = 5

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
      INCF1 = NI
      INCF2 = 1
      IY = 1+NAS2*IOFCD(ISYM2,ISYM1)
      INCY1 = 1
      INCY2 = NAS2*NI
      INCY3 = NAS2
      LEN1 = NSSH(ISYM1)
      LEN2 = NI
      call MLTMV(IMLTOP,LIST(LLST1),NLST1,X1(IXTA),size(X1(IXTA:)),FIT(ISYM12)%A,size(FIT(ISYM12)%A),GA_Arrays(lg_Y)%A(IY), &
                 size(GA_Arrays(lg_Y)%A(IY:)))
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
      INCF1 = NI
      INCF2 = 1
      IY = 1+NAS2*IOFCD(ISYM2,ISYM1)
      INCY1 = 1
      INCY2 = NAS2*NI
      INCY3 = NAS2
      LEN1 = NSSH(ISYM1)
      LEN2 = NI
      call MLTMV(IMLTOP,LIST(LLST1),NLST1,X2(IX),size(X2(IX:)),FIT(ISYM12)%A,size(FIT(ISYM12)%A),GA_Arrays(lg_Y)%A(IY), &
                 size(GA_Arrays(lg_Y)%A(IY:)))
    end if
    ! -----------------------------------------------
  case (9)
    !ICASE1 = 4
    !ICASE2 = 8

    ! C&FP One-el
    NLST1 = NLIST(ISYM1,ISYM2,12)
    NLST2 = NLIST(ISYM1,ISYM2,16)
    if (NLST1*NLST2 /= 0) then
      LLST1 = LLIST(ISYM1,ISYM2,12)
      VAL1(1) = -One
      VAL1(2) = -Two
      LLST2 = LLIST(ISYM1,ISYM2,16)
      VAL2(1) = One
      VAL2(2) = SQR2
      IXTA = 1
      INCX1 = 1
      INCX2 = NASH(ISYM1)
      INCF1 = 1
      INCF2 = NASH(ISYM12)
      IY = 1
      INCY1 = 1
      INCY2 = NTGEU(ISYM2)
      call MLTSCA(IMLTOP,LIST(LLST1),NLST1,LIST(LLST2),NLST2,X1(IXTA),size(X1(IXTA:)),FTA(ISYM12)%A,size(FTA(ISYM12)%A), &
                  GA_Arrays(lg_Y)%A(IY),size(GA_Arrays(lg_Y)%A(IY:)))
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
      INCF1 = 1
      INCF2 = NASH(ISYM12)
      IY = 1
      INCY1 = 1
      INCY2 = NTGEU(ISYM2)
      call MLTSCA(IMLTOP,LIST(LLST1),NLST1,LIST(LLST2),NLST2,X2(IX),size(X2(IX:)),FTA(ISYM12)%A,size(FTA(ISYM12)%A), &
                  GA_Arrays(lg_Y)%A(IY),size(GA_Arrays(lg_Y)%A(IY:)))
    end if
    ! -----------------------------------------------
  case (10)
    !ICASE1 = 4
    !ICASE2 = 9

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
      INCF1 = 1
      INCF2 = NASH(ISYM12)
      IY = 1
      INCY1 = 1
      INCY2 = NTGTU(ISYM2)
      call MLTSCA(IMLTOP,LIST(LLST1),NLST1,LIST(LLST2),NLST2,X1(IXTA),size(X1(IXTA:)),FTA(ISYM12)%A,size(FTA(ISYM12)%A), &
                  GA_Arrays(lg_Y)%A(IY),size(GA_Arrays(lg_Y)%A(IY:)))
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
      INCF1 = 1
      INCF2 = NASH(ISYM12)
      IY = 1
      INCY1 = 1
      INCY2 = NTGTU(ISYM2)
      call MLTSCA(IMLTOP,LIST(LLST1),NLST1,LIST(LLST2),NLST2,X2(IX),size(X2(IX:)),FTA(ISYM12)%A,size(FTA(ISYM12)%A), &
                  GA_Arrays(lg_Y)%A(IY),size(GA_Arrays(lg_Y)%A(IY:)))
    end if
    ! -----------------------------------------------
  case (11)
    !ICASE1 = 4
    !ICASE2 = 10

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
          INCF1 = NI
          INCF2 = 1
          IY = 1+NAS2*IOFCGP(ISYM2,ISYMI)
          INCY1 = NT*NI
          INCY2 = 1
          INCY3 = NT
          LEN1 = NT
          LEN2 = NI
          call MLTMV(IMLTOP,LIST(LLST1),NLST1,X1(IXTA),size(X1(IXTA:)),FIA(ISYMI)%A,size(FIA(ISYMI)%A),GA_Arrays(lg_Y)%A(IY), &
                     size(GA_Arrays(lg_Y)%A(IY:)))
        end if
      end do
    end if
    ! -----------------------------------------------
  case (12)
    !ICASE1 = 4
    !ICASE2 = 11

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
          INCF1 = NI
          INCF2 = 1
          IY = 1+NAS2*IOFCGM(ISYM2,ISYMI)
          INCY1 = NT*NI
          INCY2 = 1
          INCY3 = NT
          LEN1 = NT
          LEN2 = NI
          call MLTMV(IMLTOP,LIST(LLST1),NLST1,X1(IXTA),size(X1(IXTA:)),FIA(ISYMI)%A,size(FIA(ISYMI)%A),GA_Arrays(lg_Y)%A(IY), &
                     size(GA_Arrays(lg_Y)%A(IY:)))
        end if
      end do
    end if
    ! -----------------------------------------------
  case (13)
    !ICASE1 = 5
    !ICASE2 = 6

    ! D&EP One-el
    NT = NASH(ISYM2)
    if ((ISYM1 == 1) .and. (NT /= 0)) then
      do ISYMI=1,NSYM
        NI = NISH(ISYMI)
        NA = NSSH(ISYMI)
        ISYMIJ = Mul(ISYMI,ISYM2)
        NLST1 = NLIST(ISYMI,ISYMIJ,14)
        if (NI*NA*NLST1 /= 0) then
          LLST1 = LLIST(ISYMI,ISYMIJ,14)
          VAL1(1) = SQRI2
          VAL1(2) = One
          IXIA = 1+IOFFIA(ISYMI)
          INCX1 = 1
          INCX2 = NI
          INCF1 = NASH(ISYM2)
          INCF2 = 1
          IY = 1+NAS2*IOFCEP(ISYM2,ISYMI)
          INCY1 = NT*NA
          INCY2 = NT
          INCY3 = 1
          LEN1 = NA
          LEN2 = NT
          call MLTMV(IMLTOP,LIST(LLST1),NLST1,X1(IXIA),size(X1(IXIA:)),FTI(ISYM2)%A,size(FTI(ISYM2)%A),GA_Arrays(lg_Y)%A(IY), &
                     size(GA_Arrays(lg_Y)%A(IY:)))
        end if
      end do
    end if

    ! D&EP Two-el
    LLST1 = LLIST(ISYM1,ISYM2,7)
    NLST1 = NLIST(ISYM1,ISYM2,7)
    NU = NASH(ISYM2)
    if (NLST1*NU /= 0) then
      VAL1(1) = -One
      VAL1(2) = -One
      do ISYMA=1,NSYM
        NA = NSSH(ISYMA)
        ISYMI = Mul(ISYMA,ISYM1)
        NI = NISH(ISYMI)
        ISYMIJ = Mul(ISYMI,ISYM12)
        NLST2 = NLIST(ISYMI,ISYMIJ,14)
        if (NA*NI*NLST2 /= 0) then
          LLST2 = LLIST(ISYMI,ISYMIJ,14)
          VAL2(1) = SQRI2
          VAL2(2) = One
          IX = 1+NAS1*IOFCD(ISYM1,ISYMA)
          INCX1 = 1
          INCX2 = NAS1
          INCX3 = NAS1*NI
          INCF1 = NISH(ISYM12)
          INCF2 = 1
          IY = 1+NU*IOFCEP(ISYM2,ISYMA)
          INCY1 = 1
          INCY2 = NU*NA
          INCY3 = NU
          LEN1 = NA
          call MLTDXP(IMLTOP,LIST(LLST1),LIST(LLST2),X2(IX:),size(X2(IX:)),FIT(ISYM12)%A(:),size(FIT(ISYM12)%A(:)), &
                      GA_Arrays(lg_Y)%A(IY:),size(GA_Arrays(lg_Y)%A(IY:)))
        end if
      end do
    end if
    ! -----------------------------------------------
  case (14)
    !ICASE1 = 5
    !ICASE2 = 7

    ! D&EM One-el
    NT = NASH(ISYM2)
    if ((ISYM1 == 1) .and. (NT /= 0)) then
      do ISYMI=1,NSYM
        NI = NISH(ISYMI)
        NA = NSSH(ISYMI)
        ISYMIJ = Mul(ISYMI,ISYM2)
        NLST1 = NLIST(ISYMI,ISYMIJ,15)
        if (NI*NA*NLST1 /= 0) then
          LLST1 = LLIST(ISYMI,ISYMIJ,15)
          VAL1(1) = SQR32
          VAL1(2) = -SQR32
          IXIA = 1+IOFFIA(ISYMI)
          INCX1 = 1
          INCX2 = NI
          INCF1 = NASH(ISYM2)
          INCF2 = 1
          IY = 1+NAS2*IOFCEM(ISYM2,ISYMI)
          INCY1 = NT*NA
          INCY2 = NT
          INCY3 = 1
          LEN1 = NA
          LEN2 = NT
          call MLTMV(IMLTOP,LIST(LLST1),NLST1,X1(IXIA),size(X1(IXIA:)),FTI(ISYM2)%A,size(FTI(ISYM12)%A),GA_Arrays(lg_Y)%A(IY), &
                     size(GA_Arrays(lg_Y)%A(IY:)))
        end if
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
          INCF1 = NISH(ISYM12)
          INCF2 = 1
          IY = 1+NU*IOFCEM(ISYM2,ISYMA)
          INCY1 = 1
          INCY2 = NU*NA
          INCY3 = NU
          LEN1 = NA
          call MLTDXP(IMLTOP,LIST(LLST1),LIST(LLST2),X2(IX:),size(X2(IX:)),FIT(ISYM12)%A(:),size(FIT(ISYM12)%A(:)), &
                      GA_Arrays(lg_Y)%A(IY:),size(GA_Arrays(lg_Y)%A(IY:)))
        end if
      end do
    end if
    ! -----------------------------------------------
  case (15)
    !ICASE1 = 5
    !ICASE2 = 10

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
          INCF1 = 1
          INCF2 = NASH(ISYM12)
          IY = 1+NU*IOFCGP(ISYM2,ISYMI)
          INCY1 = 1
          INCY2 = NU*NI
          INCY3 = NU
          LEN1 = NI
          call MLTDXP(IMLTOP,LIST(LLST1),LIST(LLST2),X2(IX:),size(X2(IX:)),FTA(ISYM12)%A(:),size(FTA(ISYM12)%A(:)), &
                      GA_Arrays(lg_Y)%A(IY:),size(GA_Arrays(lg_Y)%A(IY:)))
        end if
      end do
    end if
    ! -----------------------------------------------
  case (16)
    !ICASE1 = 5
    !ICASE2 = 11

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
          INCF1 = 1
          INCF2 = NASH(ISYM12)
          IY = 1+NU*IOFCGM(ISYM2,ISYMI)
          INCY1 = 1
          INCY2 = NU*NI
          INCY3 = NU
          LEN1 = NI
          call MLTDXP(IMLTOP,LIST(LLST1),LIST(LLST2),X2(IX:),size(X2(IX:)),FTA(ISYM12)%A(:),size(FTA(ISYM12)%A(:)), &
                      GA_Arrays(lg_Y)%A(IY:),size(GA_Arrays(lg_Y)%A(IY:)))
        end if
      end do
    end if
    ! -----------------------------------------------
  case (17)
    !ICASE1 = 6
    !ICASE2 = 12

    ! EP&HP Two-el
    NA = NSSH(ISYM12)
    NLST1 = NLIST(ISYM12,ISYM2,16)
    if (NA*NLST1 /= 0) then
      LLST1 = LLIST(ISYM12,ISYM2,16)
      VAL1(1) = SQRI2
      VAL1(2) = One
      JXOFF = IOFCEP(ISYM1,ISYM12)
      INCX3 = NAS1*NA
      NFT = NASH(ISYM1)
      NFA = NSSH(ISYM1)
      call PMLTR1(KOD,IMLTOP,LIST(LLST1),X2(:),size(X2(:)),NAS1,NIS1,JXOFF,FTA(ISYM1)%A(:),NFT,NFA,lg_Y,NAS2,NIS2)
    end if
    ! -----------------------------------------------
  case (18)
    !ICASE1 = 7
    !ICASE2 = 13

    ! EM&HM Two-el
    NA = NSSH(ISYM12)
    NLST1 = NLIST(ISYM12,ISYM2,17)
    if (NA*NLST1 /= 0) then
      LLST1 = LLIST(ISYM12,ISYM2,17)
      VAL1(1) = SQRI2
      VAL1(2) = -SQRI2
      JXOFF = IOFCEM(ISYM1,ISYM12)
      INCX3 = NAS1*NA
      NFT = NASH(ISYM1)
      NFA = NSSH(ISYM1)
      call PMLTR1(KOD,IMLTOP,LIST(LLST1),X2(:),size(X2(:)),NAS1,NIS1,JXOFF,FTA(ISYM1)%A(:),NFT,NFA,lg_Y,NAS2,NIS2)
    end if
    ! -----------------------------------------------
  case (19)
    !ICASE1 = 8
    !ICASE2 = 10

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
      INCF1 = NI
      INCF2 = 1
      IY = 1+NAS2*IOFCGP(ISYM2,ISYM12)
      INCY1 = 1
      INCY2 = NAS2*NI
      INCY3 = NAS2
      LEN1 = NIS1
      LEN2 = NI
      call MLTMV(IMLTOP,LIST(LLST1),NLST1,X2(IX),size(X2(IX:)),FIT(ISYM12)%A,size(FIT(ISYM12)%A),GA_Arrays(lg_Y)%A(IY), &
                 size(GA_Arrays(lg_Y)%A(IY:)))
    end if
    ! -----------------------------------------------
  case (20)
    !ICASE1 = 9
    !ICASE2 = 11

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
      INCF1 = NI
      INCF2 = 1
      IY = 1+NAS2*IOFCGM(ISYM2,ISYM12)
      INCY1 = 1
      INCY2 = NAS2*NI
      INCY3 = NAS2
      LEN1 = NIS1
      LEN2 = NI
      call MLTMV(IMLTOP,LIST(LLST1),NLST1,X2(IX),size(X2(IX:)),FIT(ISYM12)%A,size(FIT(ISYM12)%A),GA_Arrays(lg_Y)%A(IY), &
                 size(GA_Arrays(lg_Y)%A(IY:)))
    end if
    ! -----------------------------------------------
  case (21)
    !ICASE1 = 10
    !ICASE2 = 12

    ! GP&HP Two-el
    NLST1 = NLIST(ISYM12,ISYM2,14)
    if (NLST1 /= 0) then
      LLST1 = LLIST(ISYM12,ISYM2,14)
      VAL1(1) = -SQRI2
      VAL1(2) = -One
      JXOFF = IOFCGP(ISYM1,ISYM12)
      INCX3 = NAS1*NISH(ISYM12)
      NFT = NASH(ISYM1)
      NFI = NISH(ISYM1)
      call PMLTR1(KOD,IMLTOP,LIST(LLST1),X2(:),size(X2(:)),NAS1,NIS1,JXOFF,FTI(ISYM1)%A(:),NFT,NFI,lg_Y,NAS2,NIS2)
    end if
    ! -----------------------------------------------
  case (22)
    !ICASE1 = 11
    !ICASE2 = 13

    ! GM&HM Two-el
    NLST1 = NLIST(ISYM12,ISYM2,15)
    if (NLST1 /= 0) then
      LLST1 = LLIST(ISYM12,ISYM2,15)
      VAL1(1) = SQRI2
      VAL1(2) = -SQRI2
      JXOFF = IOFCGM(ISYM1,ISYM12)
      INCX3 = NAS1*NISH(ISYM12)
      NFT = NASH(ISYM1)
      NFI = NISH(ISYM1)
      call PMLTR1(KOD,IMLTOP,LIST(LLST1),X2(:),size(X2(:)),NAS1,NIS1,JXOFF,FTI(ISYM1)%A(:),NFT,NFI,lg_Y,NAS2,NIS2)
    end if
    ! -----------------------------------------------
  case (23)
    !ICASE1 = 5
    !ICASE2 = 12

    ! D&HP One-el
    if (ISYM1 == 1) then
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
          IXIA = 1+IOFFIA(ISYMI)
          INCX1 = 1
          INCX2 = NI
          call PMLTSCA(KOD,IMLTOP,LIST(LLST1),LIST(LLST2),X1(IXIA),NI,NA,FIA(ISYMB)%A,NJ,NB,lg_Y,NAS2,NIS2)
        end if
      end do
    end if
    !---------------------------
  case (24)
    !ICASE1 = 5
    !ICASE2 = 13

    ! D&HM One-el
    if (ISYM1 == 1) then
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
          IXIA = 1+IOFFIA(ISYMI)
          INCX1 = 1
          INCX2 = NI
          call PMLTSCA(KOD,IMLTOP,LIST(LLST1),LIST(LLST2),X1(IXIA),NI,NA,FIA(ISYMB)%A,NJ,NB,lg_Y,NAS2,NIS2)
        end if
      end do
    end if
    !---------------------------
  case default
    write(u6,*) ' INTERNAL ERROR: SGM reached invalid KOD=',KOD
    call Abend()
end select

end subroutine SGM
