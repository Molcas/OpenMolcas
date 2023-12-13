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

subroutine CHO_DBGINT_CHO(XINT,NCD,NAB,WRK,LWRK,ERRMAX,ERRMIN,ERRRMS,NCMP,ISHLCD,ISHLAB)
!
! Purpose: calculate integrals in shell quadruple (CD|AB) from
!          Cholesky vectors on disk and compare to those in
!          XINT (for debugging).
!
! NOTE: this is *only* for debugging.

use Index_Functions, only: nTri_Elem
use Cholesky, only: iiBstR, iiBstRSh, IndRed, iSP2F, nBstSh, nnBstR, nnBstRSh, nSym, nSys_call, NumCho
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NCD, NAB, LWRK, ISHLCD, ISHLAB
real(kind=wp), intent(in) :: XINT(NCD,NAB)
real(kind=wp), intent(out) :: WRK(LWRK), ERRMAX, ERRMIN, ERRRMS
integer(kind=iwp), intent(out) :: NCMP
integer(kind=iwp) :: IAB, IBATCH, ICD, ICDAB, ISHLA, ISHLB, ISHLC, ISHLD, ISYM, IVEC, JAB, JCD, JVEC, JVEC1, KAB, KCD, KCDAB, &
                     KCHOAB, KCHOCD, KEND0, KEND1, KEND2, KINT, KOFF1, KOFF2, KREAD, KVEC1, KXINT, LCDABT, LENint, LREAD, LVEC1, &
                     LWRK0, LWRK1, LWRK2, MINM, NABL, NBATCH, NCDL, NSCALL, NUMAB, NUMCD, NUMV, NVEC
real(kind=wp) :: DIFF
character(len=*), parameter :: SECNAM = 'CHO_DBGINT_CHO'
integer(kind=iwp), external :: CHO_LREAD

! Initializations.
! ----------------

call CHO_INVPCK(ISP2F(ISHLCD),ISHLC,ISHLD,.true.)
call CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.true.)
ERRMAX = -1.0e12_wp
ERRMIN = 1.0e12_wp
ERRRMS = Zero
NCMP = 0
LCDABT = NCD*NAB

if (ISHLC == ISHLD) then
  NCDL = nTri_Elem(NBSTSH(ISHLC))
else
  NCDL = NBSTSH(ISHLC)*NBSTSH(ISHLD)
end if
if (ISHLA == ISHLB) then
  NABL = nTri_Elem(NBSTSH(ISHLA))
else
  NABL = NBSTSH(ISHLA)*NBSTSH(ISHLB)
end if
if (NCDL > NCD) call CHO_QUIT('NCD error in '//SECNAM,104)
if (NABL > NAB) call CHO_QUIT('NAB error in '//SECNAM,104)
if ((NAB < 1) .or. (NCD < 1)) return

! Save read-call counter.
! -----------------------

NSCALL = NSYS_CALL

! Get a copy of XINT.
! -------------------

KXINT = 1
KEND0 = KXINT+LCDABT
LWRK0 = LWRK-KEND0
if (LWRK0 <= 0) call CHO_QUIT('Insufficient memory in '//SECNAM//' [0]',101)

WRK(1:LCDABT) = reshape(XINT(:,:),[LCDABT])

! Start symmetry loop.
! --------------------

do ISYM=1,NSYM

  NUMCD = NNBSTRSH(ISYM,ISHLCD,2)
  NUMAB = NNBSTRSH(ISYM,ISHLAB,2)

  if ((NUMCD > 0) .and. (NUMAB > 0) .and. (NUMCHO(ISYM) > 0)) then

    ! Allocate space for integrals and for Cholesky reading.
    ! ------------------------------------------------------

    LENint = NUMCD*NUMAB
    LREAD = CHO_LREAD(ISYM,LWRK)
    LVEC1 = NNBSTR(ISYM,2)

    KINT = KEND0
    KREAD = KINT+LENint
    KVEC1 = KREAD+LREAD
    KEND1 = KVEC1+LVEC1
    LWRK1 = LWRK-KEND1+1

    if (LWRK1 <= 0) call CHO_QUIT('Insufficient memory in '//SECNAM,104)

    ! Initialize integral array.
    ! --------------------------

    WRK(KINT:KINT+LENint-1) = Zero

    ! Set up batch over Cholesky vectors.
    ! -----------------------------------

    MINM = NUMCD+NUMAB
    NVEC = min(LWRK1/MINM,NUMCHO(ISYM))
    if (NVEC < 1) call CHO_QUIT('Batch problem in '//SECNAM,104)
    NBATCH = (NUMCHO(ISYM)-1)/NVEC+1

    ! Start batch loop.
    ! -----------------

    do IBATCH=1,NBATCH

      if (IBATCH == NBATCH) then
        NUMV = NUMCHO(ISYM)-NVEC*(NBATCH-1)
      else
        NUMV = NVEC
      end if
      JVEC1 = NVEC*(IBATCH-1)+1

      KCHOCD = KEND1
      KCHOAB = KCHOCD+NUMCD*NUMV
      KEND2 = KCHOAB+NUMAB*NUMV
      LWRK2 = LWRK-KEND2+1

      if (LWRK2 < 0) call CHO_QUIT('Batch error in '//SECNAM,104)

      ! Read vectors.
      ! -------------

      do IVEC=1,NUMV
        JVEC = JVEC1+IVEC-1
        call CHO_GETVEC(WRK(KVEC1),LVEC1,1,JVEC,ISYM,WRK(KREAD),LREAD)
        KOFF1 = KVEC1+IIBSTRSH(ISYM,ISHLCD,2)
        KOFF2 = KCHOCD+NUMCD*(IVEC-1)
        WRK(KOFF2:KOFF2+NUMCD-1) = WRK(KOFF1:KOFF1+NUMCD-1)
        KOFF1 = KVEC1+IIBSTRSH(ISYM,ISHLAB,2)
        KOFF2 = KCHOAB+NUMAB*(IVEC-1)
        WRK(KOFF2:KOFF2+NUMAB-1) = WRK(KOFF1:KOFF1+NUMAB-1)
      end do

      ! Calculate contribution.
      ! -----------------------

      call DGEMM_('N','T',NUMCD,NUMAB,NUMV,One,WRK(KCHOCD),NUMCD,WRK(KCHOAB),NUMAB,One,WRK(KINT),NUMCD)

    end do

    ! Subtract contribution from full shell pair.
    ! -------------------------------------------

    do IAB=1,NUMAB
      JAB = IIBSTR(ISYM,2)+IIBSTRSH(ISYM,ISHLAB,2)+IAB
      KAB = INDRED(INDRED(JAB,2),1)
      do ICD=1,NUMCD
        JCD = IIBSTR(ISYM,2)+IIBSTRSH(ISYM,ISHLCD,2)+ICD
        KCD = INDRED(INDRED(JCD,2),1)
        ICDAB = KINT+NUMCD*(IAB-1)+ICD-1
        KCDAB = KXINT+NCD*(KAB-1)+KCD-1
        WRK(KCDAB) = WRK(KCDAB)-WRK(ICDAB)
      end do
    end do

  end if

end do

! Compare full shell pair.
! ------------------------

do KAB=1,NAB
  do KCD=1,NCD
    KCDAB = KXINT+NCD*(KAB-1)+KCD-1
    DIFF = WRK(KCDAB)
    NCMP = NCMP+1
    if (NCMP == 1) then
      ERRMAX = DIFF
      ERRMIN = DIFF
    else
      if (abs(DIFF) > abs(ERRMAX)) ERRMAX = DIFF
      if (abs(DIFF) < abs(ERRMIN)) ERRMIN = DIFF
    end if
    ERRRMS = ERRRMS+DIFF*DIFF
  end do
end do

! Restore read-call counter.
! --------------------------

NSYS_CALL = NSCALL

end subroutine CHO_DBGINT_CHO
