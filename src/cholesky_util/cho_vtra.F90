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

subroutine CHO_VTRA(irc,scr,lscr,jVref,JVEC1,JNUM,NUMV,JSYM,IREDC,iSwap,nDen,kDen,MOs,ChoT)
!********************************************************
! Author: F. Aquilante
!
! Purpose:  SCR(lscr) contains JNUM cholesky vectors
!           starting from JVEC1 and stored in reduced
!           sets. The routine performs an MOs half transformation
!           of these elements in a set of target
!           arrays, ChoT, identified the data type SBA_Type.
!           In the target arrays, the vectors are
!           stored in full dimension and as a
!           subset of a a given NUMV number of vectors.
!
! Input:
!     jVref =  index of the first vector to be transformed
!              computed wrt to the first vector in the
!              target arrays (see calling routine)
!
!     Jvec1 =  first vector to be transformed
!     JNum  =  # of vectors to be transformed
!
!     NumV  =  total # of vectors in the target arrays
!
!     nDen  =  total # of densities to which MOs refer
!     kDen  =  first density to be treated
!
!     iSwap :   = 0   L(k,b,J) is returned
!               = 1   L(a,k,J) is returned
!               = 2   L(k,J,b) is returned
!               = 3   L(a,J,k) is returned
!
!     IREDC :  reduced set in core at the moment of
!              the first call to the routine.
!              Can be set to -1 by the calling routine
!
!********************************************************

use Symmetry_Info, only: Mul
use Cholesky, only: iBas, iiBstR, IndRed, InfVec, iRS2F, nBas, nDimRS, nnBstR, nSym
use Data_Structures, only: DSBA_Type, SBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: lScr, jVref, JVEC1, JNUM, NUMV, JSYM, iSWap, nDen, kDen
real(kind=wp), intent(in) :: Scr(lscr)
integer(kind=iwp), intent(inout) :: IREDC
type(DSBA_Type), intent(in) :: MOs(nDen)
type(SBA_Type), intent(inout) :: ChoT(nDen)
integer(kind=iwp) :: iag, ias, ibg, ibs, iDen, iE, ij, iLoc, iRab, iSym, iSyma, iSymb, jDen, jRab, JRED, JVEC, kRab, kscr, kVEC, &
                     LVEC, n1, NREAD
real(kind=wp) :: xfd
integer(kind=iwp), allocatable :: nPorb(:,:)
real(kind=wp), parameter :: Fac(0:1) = [Half,One]
character(len=*), parameter :: SECNAM = 'CHO_VTRA'
integer(kind=iwp), external :: cho_isao

call mma_allocate(nPorb,8,nDen,Label='nPorb')
do iDen=1,nDen
  do iSym=1,nSym
    if (associated(MOs(iDen)%SB(iSym)%A2)) then
      nPorb(iSym,iDen) = size(MOs(iDen)%SB(iSym)%A2,1)
    else
      nPorb(iSym,iDen) = 0
    end if
  end do
end do

!*********************************************************
!
! From Reduced sets to half-MOs full storage
! ------------------------------------------
!
! iSwap = 0
!
!  L{a,b,J} ---> L(p,b,J)   ! stride-1 transformation
!
! iSwap = 1
!
!  L{a,b,J} ---> L(a,q,J)
!
! iSwap = 2
!
!  L{a,b,J} ---> L(p,J,b)   ! stride-1 transformation
!
! iSwap = 3
!
!  L{a,b,J} ---> L(a,J,q)
!*********************************************************

iLoc = 3 ! use scratch location in reduced index arrays

!                                                                      *
!***********************************************************************
!                                                                      *
select case (iSwap)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  case (0)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    NREAD = 0
    do JVEC=1,JNUM   ! Relative index in the JNUM batch

      LVEC = JVEC-1+JVREF  ! Relative index in the NUMV batch
      kVEC = JVEC-1+JVEC1  ! Absolute index
      JRED = InfVec(KVEC,2,jSym)

      if (JRED /= IREDC) then ! JRED is not the rs in core
        call Cho_X_SetRed(irc,iLoc,JRED) !set indx arrays at iLoc
        IREDC = JRED
      end if

      kscr = NREAD
      NREAD = NREAD+nDimRS(jSym,JRED)

      if (JSYM == 1) then  ! L(a,b,J)=L(b,a,J); only a >= b presnt

        do jRab=1,nnBstR(jSym,iLoc)

          kRab = iiBstr(jSym,iLoc)+jRab
          iRab = IndRed(kRab,iLoc)

          iag = iRS2F(1,iRab)  !global address
          ibg = iRS2F(2,iRab)

          iSyma = cho_isao(iag) !symmetry block; S(b)=S(a)=S(p)

          kscr = kscr+1

          ias = iag-ibas(iSyma) !address within that sym
          ibs = ibg-ibas(iSyma)
          xfd = Fac(min(abs(ias-ibs),1)) !fac for diag

          ! L(p,b,J) = sum_a  xfd* L(a,b,J) * C(p,a)
          ! ----------------------------------------
          do jDen=kDen,nDen

            if (.not. associated(ChoT(jDen)%SB(iSyma)%A3)) cycle

            ! C(1,b)
            ChoT(jDen)%SB(iSyma)%A3(:,ias,LVEC) = ChoT(jDen)%SB(iSyma)%A3(:,ias,LVEC)+xfd*Scr(kscr)*MOs(JDen)%SB(iSyma)%A2(:,ibs)

            ! C(1,a)
            ChoT(jDen)%SB(iSyma)%A3(:,ibs,LVEC) = ChoT(jDen)%SB(iSyma)%A3(:,ibs,LVEC)+xfd*Scr(kscr)*MOs(JDen)%SB(iSyma)%A2(:,ias)

          end do  ! loop over densities

        end do  ! jRab loop

      else  ! jSym /= 1

        do jRab=1,nnBstR(jSym,iLoc)

          kRab = iiBstr(jSym,iLoc)+jRab
          iRab = IndRed(kRab,iLoc)

          iag = iRS2F(1,iRab)  !global address
          ibg = iRS2F(2,iRab)

          iSyma = cho_isao(iag)  !symmetry block
          iSymb = mul(jSym,iSyma)

          kscr = kscr+1

          ias = iag-ibas(iSyma)  !address within that sym
          ibs = ibg-ibas(iSymb)

          ! L(p,b,J) = sum_a  L(a,b,J) * C(p,a)
          ! -----------------------------------
          do jDen=kDen,nDen

            if (.not. associated(ChoT(jDen)%SB(iSyma)%A3)) cycle

            ChoT(jDen)%SB(iSyma)%A3(:,ibs,LVEC) = ChoT(jDen)%SB(iSyma)%A3(:,ibs,LVEC)+Scr(kscr)*MOs(jDen)%SB(iSyma)%A2(:,ias)

          end do

          ! L(p,a,J) = sum_b  L(a,b,J) * C(p,b)
          ! -----------------------------------
          do jDen=kDen,nDen

            if (.not. associated(ChoT(jDen)%SB(iSymb)%A3)) cycle

            ChoT(jDen)%SB(iSymb)%A3(:,ias,LVEC) = ChoT(jDen)%SB(iSymb)%A3(:,ias,LVEC)+Scr(kscr)*MOs(jDen)%SB(iSymb)%A2(:,ibs)

          end do

        end do  ! jRab loop

      end if ! total symmetric vectors check

    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  case (1)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    NREAD = 0
    do JVEC=1,JNUM

      LVEC = JVEC-1+JVREF
      kVEC = JVEC-1+JVEC1  ! Absolute index
      JRED = InfVec(KVEC,2,jSym)

      if (JRED /= IREDC) then ! JRED is not the reduced set in core
        call Cho_X_SetRed(irc,iLoc,JRED) !set index arrays at iLoc
        IREDC = JRED
      end if

      kscr = NREAD
      NREAD = NREAD+nDimRS(jSym,JRED)

      if (JSYM == 1) then  ! L(a,b,J)=L(b,a,J); only a >= b present

        do jRab=1,nnBstR(jSym,iLoc)

          kRab = iiBstr(jSym,iLoc)+jRab
          iRab = IndRed(kRab,iLoc)

          iag = iRS2F(1,iRab)  !global address
          ibg = iRS2F(2,iRab)

          iSyma = cho_isao(iag)  !symmetry block; Sym(b)=Sym(a)=Sym(p)

          kscr = kscr+1

          ias = iag-ibas(iSyma)  !address within that symm block
          ibs = ibg-ibas(iSyma)
          xfd = Fac(min(abs(ias-ibs),1)) !fac for diagonal elements

          ! L(a,p,J) = sum_b  xfd* L(a,b,J) * C(p,b)
          ! ----------------------------------------
          do jDen=kDen,nDen

            if (.not. associated(ChoT(jDen)%SB(iSyma)%A2)) cycle

            iE = size(ChoT(jDen)%SB(iSyma)%A2,1)
            call DAXPY_(nPorb(iSyma,jDen),xfd*Scr(kscr),MOs(jDen)%SB(iSyma)%A2(:,ibs),1,ChoT(jDen)%SB(iSyma)%A2(ias:iE,LVEC), &
                        nBas(iSyma))

            call DAXPY_(nPorb(iSyma,jDen),xfd*Scr(kscr),MOs(jDen)%SB(iSyma)%A2(:,ias),1,ChoT(jDen)%SB(iSyma)%A2(ibs:iE,LVEC), &
                        nBas(iSyma))

          end do  ! loop over densities

        end do  ! jRab loop

      else  ! jSym /= 1

        do jRab=1,nnBstR(jSym,iLoc)

          kRab = iiBstr(jSym,iLoc)+jRab
          iRab = IndRed(kRab,iLoc)

          iag = iRS2F(1,iRab)  !global address
          ibg = iRS2F(2,iRab)

          iSyma = cho_isao(iag)  !symmetry block
          iSymb = mul(jSym,iSyma)

          kscr = kscr+1

          ias = iag-ibas(iSyma)  !address within that symm block
          ibs = ibg-ibas(iSymb)

          ! L(a,q,J) = sum_b  L(a,b,J) * C(q,b)
          ! -----------------------------------
          do jDen=kDen,nDen

            if (.not. associated(ChoT(jDen)%SB(iSyma)%A2)) cycle

            iE = size(ChoT(jDen)%SB(iSyma)%A2,1)
            call DAXPY_(nPorb(iSymb,jDen),Scr(kscr),MOs(jDen)%SB(iSymb)%A2(:,ibs),1,ChoT(jDen)%SB(iSyma)%A2(ias:iE,LVEC), &
                        nBas(iSyma))

          end do

          ! L(b,q,J) = sum_a  L(a,b,J) * C(q,a)
          ! -----------------------------------
          do jDen=kDen,nDen

            if (.not. associated(ChoT(jDen)%SB(iSymb)%A2)) cycle

            iE = size(ChoT(jDen)%SB(iSymb)%A2,1)
            call DAXPY_(nPorb(iSyma,jDen),Scr(kscr),MOs(jDen)%SB(iSyma)%A2(:,ias),1,ChoT(jDen)%SB(iSymb)%A2(ibs:iE,LVEC), &
                        nBas(iSyma))

          end do

        end do  ! jRab loop

      end if ! total symmetric vectors check

    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  case (2)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    NREAD = 0
    do JVEC=1,JNUM

      LVEC = JVEC-1+JVREF
      kVEC = JVEC-1+JVEC1  ! Absolute index
      JRED = InfVec(KVEC,2,jSym)

      if (JRED /= IREDC) then
        call Cho_X_SetRed(irc,iLoc,JRED)
        IREDC = JRED
      end if

      kscr = NREAD
      NREAD = NREAD+nDimRS(jSym,JRED)

      if (JSYM == 1) then  ! L(a,b,J)=L(b,a,J); only a >= b

        do jRab=1,nnBstR(jSym,iLoc)

          kRab = iiBstr(jSym,iLoc)+jRab
          iRab = IndRed(kRab,iLoc)

          iag = iRS2F(1,iRab)  !global address
          ibg = iRS2F(2,iRab)

          iSyma = cho_isao(iag)

          kscr = kscr+1

          ias = iag-ibas(iSyma)
          ibs = ibg-ibas(iSyma)
          xfd = Fac(min(abs(ias-ibs),1))

          ! L(p,J,b) = sum_a  xfd* L(a,b,J) * C(p,a)
          ! ----------------------------------------
          do jDen=kDen,nDen

            if (.not. associated(ChoT(jDen)%SB(iSyma)%A3)) cycle

            ChoT(jDen)%SB(iSyma)%A3(:,LVEC,ias) = ChoT(jDen)%SB(iSyma)%A3(:,LVEC,ias)+xfd*Scr(kscr)*MOs(jDen)%SB(iSyma)%A2(:,ibs)

            ChoT(jDen)%SB(iSyma)%A3(:,LVEC,ibs) = ChoT(jDen)%SB(iSyma)%A3(:,LVEC,ibs)+xfd*Scr(kscr)*MOs(jDen)%SB(iSyma)%A2(:,ias)

          end do  ! loop over densities

        end do  ! jRab loop

      else  ! jSym /= 1

        do jRab=1,nnBstR(jSym,iLoc)

          kRab = iiBstr(jSym,iLoc)+jRab
          iRab = IndRed(kRab,iLoc)

          iag = iRS2F(1,iRab)  !global address
          ibg = iRS2F(2,iRab)

          iSyma = cho_isao(iag)  !symmetry block
          iSymb = mul(jSym,iSyma)  !(syma>symb)

          kscr = kscr+1

          ias = iag-ibas(iSyma)
          ibs = ibg-ibas(iSymb)

          ! L(p,J,b) = sum_a  L(a,b,J) * C(p,a)
          ! -----------------------------------
          do jDen=kDen,nDen

            if (.not. associated(ChoT(jDen)%SB(iSyma)%A3)) cycle

            ChoT(jDen)%SB(iSyma)%A3(:,LVEC,ibs) = ChoT(jDen)%SB(iSyma)%A3(:,LVEC,ibs)+Scr(kscr)*MOs(jDen)%SB(iSyma)%A2(:,ias)

          end do

          ! L(p,J,a) = sum_b  L(a,b,J) * C(p,b)
          ! -----------------------------------
          do jDen=kDen,nDen

            if (.not. associated(ChoT(jDen)%SB(iSymb)%A3)) cycle

            ChoT(jDen)%SB(iSymb)%A3(:,LVEC,ias) = ChoT(jDen)%SB(iSymb)%A3(:,LVEC,ias)+Scr(kscr)*MOs(jDen)%SB(iSymb)%A2(:,ibs)

          end do

        end do  ! jRab loop

      end if ! total symmetric vectors check

    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  case (3)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    NREAD = 0
    do JVEC=1,JNUM

      LVEC = JVEC-1+JVREF
      KVEC = JVEC-1+JVEC1  ! Absolute index
      JRED = InfVec(KVEC,2,jSym)

      if (JRED /= IREDC) then ! JRED is not the reduced set in core
        call Cho_X_SetRed(irc,iLoc,JRED) !set index arrays at iLoc
        IREDC = JRED
      end if

      kscr = NREAD
      NREAD = NREAD+nDimRS(jSym,JRED)

      if (JSYM == 1) then  ! L(a,b,J)=L(b,a,J); only a >= b present

        do jRab=1,nnBstR(jSym,iLoc)

          kRab = iiBstr(jSym,iLoc)+jRab
          iRab = IndRed(kRab,iLoc)

          iag = iRS2F(1,iRab)  !global address
          ibg = iRS2F(2,iRab)

          iSyma = cho_isao(iag)  !symmetry block; Sym(b)=Sym(a)=Sym(p)

          kscr = kscr+1

          ias = iag-ibas(iSyma)  !address within that sym block
          ibs = ibg-ibas(iSyma)
          xfd = Fac(min(abs(ias-ibs),1)) !scale fac for diag

          ! L(a,J,p) = sum_b  xfd* L(a,b,J) * C(p,b)
          ! ----------------------------------------
          do jDen=kDen,nDen

            if (.not. associated(ChoT(jDen)%SB(iSyma)%A1)) cycle

            n1 = size(ChoT(jDen)%SB(iSyma)%A3,1)
            ij = ias+n1*(LVEC-1)
            call DAXPY_(nPorb(iSyma,jDen),xfd*Scr(kscr),MOs(jDen)%SB(iSyma)%A2(:,ibs),1,ChoT(jDen)%SB(iSyma)%A1(ij:), &
                        nBas(iSyma)*NUMV)

            ij = ibs+n1*(LVEC-1)
            call DAXPY_(nPorb(iSyma,jDen),xfd*Scr(kscr),MOs(jDen)%SB(iSyma)%A2(:,ias),1,ChoT(jDen)%SB(iSyma)%A1(ij:), &
                        nBas(iSyma)*NUMV)

          end do  ! loop over densities

        end do  ! jRab loop

      else  ! jSym /= 1

        do jRab=1,nnBstR(jSym,iLoc)

          kRab = iiBstr(jSym,iLoc)+jRab
          iRab = IndRed(kRab,iLoc)

          iag = iRS2F(1,iRab)  !global address
          ibg = iRS2F(2,iRab)

          iSyma = cho_isao(iag)  !symmetry block
          iSymb = mul(jSym,iSyma)

          kscr = kscr+1

          ias = iag-ibas(iSyma)  !address within that symm block
          ibs = ibg-ibas(iSymb)

          ! L(a,J,q) = sum_b  L(a,b,J) * C(q,b)
          ! -----------------------------------
          do jDen=kDen,nDen

            if (.not. associated(ChoT(jDen)%SB(iSyma)%A1)) cycle

            n1 = size(ChoT(jDen)%SB(iSyma)%A3,1)
            ij = ias+n1*(LVEC-1)
            call DAXPY_(nPorb(iSymb,jDen),Scr(kscr),MOs(jDen)%SB(iSymb)%A2(:,ibs),1,ChoT(jDen)%SB(iSyma)%A1(ij:),nBas(iSyma)*NUMV)

          end do

          ! L(b,J,q) = sum_a  L(a,b,J) * C(q,a)
          ! -----------------------------------
          do jDen=kDen,nDen

            if (.not. associated(ChoT(jDen)%SB(iSymb)%A1)) cycle

            n1 = size(ChoT(jDen)%SB(iSymb)%A3,1)
            ij = ibs+n1*(LVEC-1)
            call DAXPY_(nPorb(iSyma,jDen),Scr(kscr),MOs(jDen)%SB(iSyma)%A2(:,ias),1,ChoT(jDen)%SB(iSymb)%A1(ij:),nBas(iSyma)*NUMV)

          end do

        end do  ! jRab loop

      end if ! total symmetric vectors check

    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  case Default
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    write(u6,*) SECNAM//': invalid argument. Iswap= ',Iswap
    irc = 66
    return

end select

call mma_deallocate(nPorb)
irc = 0

return

end subroutine CHO_VTRA
