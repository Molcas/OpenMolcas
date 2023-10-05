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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************

subroutine CHO_REORDR(irc,scr,lscr,jVref,JVEC1,JNUM,NUMV,JSYM,IREDC,iSwap,ipChoV,Arr,iSkip)
!***********************************************************
! Author: F. Aquilante
!
! Purpose:  SCR(lscr) contains JNUM cholesky vectors
!           starting from JVEC1 and stored in reduced
!           sets. The routine performs a reallocation
!           of these elements in a set of target
!           arrays identified by the pointers ipChoV.
!           In the target arrays, the vectors are
!           stored in full dimension and as a
!           subset of a a given NUMV number of vectors.
!           Each pointer should thereby point to a
!           location where the corresponding Cholesky
!           vector of a given unique symmetry pair
!           of indices has to be stored
!
! Input:
!     Ivec1 =  first vector to be copied
!     JNum  =  # of vectors to be copied
!
!     NumV  =  total # of vectors in the target arrays
!
!     iSwap :   = 0   L(a,b,J) is returned
!                     (in LT-storage if sym(a)=sym(b))
!               = 1   L(a,J,b) is returned
!               = 2   L(a,J,b) is returned
!                     (in SQ-storage if sym(a)=sym(b))
!
!     iSkip(syma)=0 : skip the symmetry block a.
!                  Any vector L(ab) or L(ba) with syma x symb=JSYM
!                  won't be returned in the target array
!
!     IREDC :  reduced set in core at the moment of
!              the call to the routine.
!              Can be set to -1 by the calling routine
!
!********************************************************

use Symmetry_Info, only: Mul
use Index_Functions, only: nTri_Elem
use Index_Functions, only: iTri
use Cholesky, only: iBas, iiBstR, IndRed, InfVec, iRS2F, nBas, nDimRS, nnBstR
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: lscr, jVref, JVEC1, JNUM, NUMV, JSYM, iSwap, ipChoV(*), iSkip(*)
real(kind=wp), intent(in) :: Scr(lscr)
integer(kind=iwp), intent(inout) :: IREDC
real(kind=wp), intent(inout) :: Arr(*)
integer(kind=iwp) :: iabf, iag, ias, ibg, ibs, iLoc, iRab, iSyma, iSymb, jRab, JRED, JVEC, kchov, kchov1, kchov2, kRab, kscr, NREAD
integer(kind=iwp), external :: cho_isao

!*********************************************************
!
! From Reduced sets to full storage
! ---------------------------------
!
! iSwap = 0
!
!  L{a,b,J} ---> L(a,b,J)
!
! iSwap = 1
!
!  L{a,b,J} ---> L(a,J,b)
!
! iSwap = 2
!
!  L{a,b,J} ---> L(ab,J)  ! with squaring of the
!                         ! "diagonal" symmetry blocks
!
!*********************************************************

iLoc = 3 ! use scratch location in reduced index arrays

if ((jSym == 1) .and. (iSwap == 0)) then  ! L(ab),J

  NREAD = 0
  kchov = 0
  kchov2 = 0
  do JVEC=1,JNUM

    JRED = InfVec(JVEC1-1+JVEC,2,jSym)

    if (JRED /= IREDC) then ! JRED is not the rs in core
      call Cho_X_SetRed(irc,iLoc,JRED) !set indx arrays at iLoc
      IREDC = JRED
    end if

    kscr = NREAD
    NREAD = NREAD+nDimRS(jSym,JRED)

    do jRab=1,nnBstR(jSym,iLoc)

      kRab = iiBstr(jSym,iLoc)+jRab
      iRab = IndRed(kRab,iLoc)

      iag = iRS2F(1,iRab)  !global address
      ibg = iRS2F(2,iRab)

      iSyma = cho_isao(iag)  !symmetry block

      kscr = kscr+1

      if (iSkip(iSyma) /= 0) then

        ias = iag-ibas(iSyma)  !addr within that symm block
        ibs = ibg-ibas(iSyma)

        iabf = iTri(ias,ibs)
        kchov = (JVEC-1)*nTri_Elem(nBas(iSyma))+iabf+ipChoV(iSyma)-1

        Arr(kchov) = Scr(kscr)

      end if

    end do

  end do

else if ((jSym == 1) .and. (iSwap == 1)) then  ! LaJ,b

  NREAD = 0
  kchov = 0
  kchov2 = 0
  do JVEC=1,JNUM

    JRED = InfVec(JVEC1-1+JVEC,2,jSym)

    if (JRED /= IREDC) then ! JRED is not the rs in core
      call Cho_X_SetRed(irc,iLoc,JRED) !set indx arrays at iLoc
      IREDC = JRED
    end if

    kscr = NREAD
    NREAD = NREAD+nDimRS(jSym,JRED)

    do jRab=1,nnBstR(jSym,iLoc)

      kRab = iiBstr(jSym,iLoc)+jRab
      iRab = IndRed(kRab,iLoc)

      iag = iRS2F(1,iRab)  !global address
      ibg = iRS2F(2,iRab)

      iSyma = cho_isao(iag)  !symmetry block

      kscr = kscr+1

      if (iSkip(iSyma) /= 0) then

        ias = iag-ibas(iSyma)  !addr within that symm block
        ibs = ibg-ibas(iSyma)

        kchov1 = nBas(iSyma)*NUMV*(ibs-1)+nBas(iSyma)*(jVref+JVEC-2)+ias+ipChoV(iSyma)-1
        kchov2 = nBas(iSyma)*NUMV*(ias-1)+nBas(iSyma)*(jVref+JVEC-2)+ibs+ipChoV(iSyma)-1

        Arr(kchov1) = Scr(kscr)
        Arr(kchov2) = Scr(kscr)

      end if

    end do

  end do

else if ((jSym == 1) .and. (iSwap == 2)) then  ! L[ab],J

  NREAD = 0
  kchov = 0
  kchov2 = 0
  do JVEC=1,JNUM

    JRED = InfVec(JVEC1-1+JVEC,2,jSym)

    if (JRED /= IREDC) then ! JRED is not the rs in core
      call Cho_X_SetRed(irc,iLoc,JRED) !set indx arrays at iLoc
      IREDC = JRED
    end if

    kscr = NREAD
    NREAD = NREAD+nDimRS(jSym,JRED)

    do jRab=1,nnBstR(jSym,iLoc)

      kRab = iiBstr(jSym,iLoc)+jRab
      iRab = IndRed(kRab,iLoc)

      iag = iRS2F(1,iRab)  !global address
      ibg = iRS2F(2,iRab)

      iSyma = cho_isao(iag)  !symmetry block

      kscr = kscr+1

      if (iSkip(iSyma) /= 0) then

        ias = iag-ibas(iSyma)  !addr within that symm block
        ibs = ibg-ibas(iSyma)

        kchov1 = nBas(iSyma)*nBas(iSyma)*(JVEC-1)+nBas(iSyma)*(ibs-1)+ias+ipChoV(iSyma)-1
        kchov2 = nBas(iSyma)*nBas(iSyma)*(JVEC-1)+nBas(iSyma)*(ias-1)+ibs+ipChoV(iSyma)-1

        Arr(kchov1) = Scr(kscr)
        Arr(kchov2) = Scr(kscr)

      end if

    end do

  end do

else if ((jSym > 1) .and. (iSwap == 0)) then

  NREAD = 0
  kchov = 0
  kchov2 = 0
  do JVEC=1,JNUM

    JRED = InfVec(JVEC1-1+JVEC,2,jSym)

    if (JRED /= IREDC) then ! JRED is not the rs in core
      call Cho_X_SetRed(irc,iLoc,JRED) !set indx arrays at iLoc
      IREDC = JRED
    end if

    kscr = NREAD
    NREAD = NREAD+nDimRS(jSym,JRED)

    do jRab=1,nnBstR(jSym,iLoc)

      kRab = iiBstr(jSym,iLoc)+jRab
      iRab = IndRed(kRab,iLoc)

      iag = iRS2F(1,iRab)  !global address
      ibg = iRS2F(2,iRab)

      iSyma = cho_isao(iag)  !symmetry block
      iSymb = Mul(jSym,iSyma) ! sym(a) > sym(b)

      kscr = kscr+1

      if (iSkip(iSyma) /= 0) then

        ias = iag-ibas(iSyma)  !addr within that symm block
        ibs = ibg-ibas(iSymb)

        kchov = nBas(iSyma)*nBas(iSymb)*(JVEC-1)+nBas(iSyma)*(ibs-1)+ias+ipChoV(iSyma)-1

        Arr(kchov) = Scr(kscr)

      end if

    end do

  end do

else if ((jSym > 1) .and. (iSwap == 1)) then

  NREAD = 0
  kchov = 0
  kchov2 = 0
  do JVEC=1,JNUM

    JRED = InfVec(JVEC1-1+JVEC,2,jSym)

    if (JRED /= IREDC) then ! JRED is not the rs in core
      call Cho_X_SetRed(irc,iLoc,JRED) !set indx arrays at iLoc
      IREDC = JRED
    end if

    kscr = NREAD
    NREAD = NREAD+nDimRS(jSym,JRED)

    do jRab=1,nnBstR(jSym,iLoc)

      kRab = iiBstr(jSym,iLoc)+jRab
      iRab = IndRed(kRab,iLoc)

      iag = iRS2F(1,iRab)  !global address
      ibg = iRS2F(2,iRab)

      iSyma = cho_isao(iag)  !symmetry block
      iSymb = Mul(jSym,iSyma) ! sym(a) > sym(b)

      kscr = kscr+1

      if (iSkip(iSyma) /= 0) then

        ias = iag-ibas(iSyma)  !addr within that symm block
        ibs = ibg-ibas(iSymb)

        kchov = nBas(iSyma)*NUMV*(ibs-1)+nBas(iSyma)*(jVref+JVEC-2)+ias+ipChoV(iSyma)-1

        Arr(kchov) = Scr(kscr)

      end if

    end do

  end do

else

  write(u6,*) 'Wrong parameters combination. JSYM,iSwap= ',JSYM,iSwap
  irc = 66
  return

end if

irc = 0

return

end subroutine CHO_REORDR
