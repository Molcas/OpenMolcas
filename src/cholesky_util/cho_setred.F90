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

subroutine CHO_SETRED(DIAG)
!
! Purpose: set next reduced set. A copy of the previous set
!          is stored in location 3.

use Cholesky, only: Cho_TrcNeg, Cho_UseAbs, Damp, iAtomShl, iiBstR, iiBstRSh, IndRed, iSP2F, LuPri, Mode_Screen, nnBstR, nnBstRSh, &
                    nnBstRT, nnShl, nSym, ScDiag, ThrCom
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: DIAG(*)
integer(kind=iwp) :: IAB, INEG, ISHLA, ISHLAB, ISHLB, ISYM, JAB, JAB1, JAB2, KAB, MSYM, NNEG
real(kind=wp) :: TST, XM
character(len=*), parameter :: SECNAM = 'CHO_SETRED'

MSYM = size(iiBstRSh,1)

! Debug print.
! ------------

if (CHO_TRCNEG) then
  write(LUPRI,*)
  write(LUPRI,*) SECNAM,': tracing of negative diagonals activated.'
  write(LUPRI,*) SECNAM,': flag SCDIAG     is ',SCDIAG
  write(LUPRI,*) SECNAM,': flag CHO_USEABS is ',CHO_USEABS
  if (SCDIAG) write(LUPRI,*) SECNAM,': MODE_SCREEN     is ',MODE_SCREEN
  write(LUPRI,*) SECNAM,': checking for negative diagonals in first reduced set:'
  NNEG = 0
  do ISYM=1,NSYM
    JAB1 = IIBSTR(ISYM,1)+1
    JAB2 = JAB1+NNBSTR(ISYM,1)-1
    INEG = 0
    do JAB=JAB1,JAB2
      if (DIAG(JAB) < Zero) INEG = INEG+1
    end do
    NNEG = NNEG+INEG
    write(LUPRI,*) SECNAM,': #negative in symmetry ',ISYM,': ',INEG
  end do
  write(LUPRI,*) SECNAM,': total #negative: ',NNEG
end if

! Copy index arrays from location 2 to location 3.
! ------------------------------------------------

call CHO_RSCOPY(2,3)

! Re-initialize index arrays at location 2.
! -----------------------------------------

IndRed(:,2) = 0
iiBstRSh(:,:,2) = 0
nnBstRSh(:,:,2) = 0
iiBstr(1:MSYM,2) = 0
nnBstr(1:MSYM,2) = 0
NNBSTRT(2) = 0

! Set new reduced set: mapping and SP counter.
! --------------------------------------------

if (SCDIAG) then  ! do screening

  if (MODE_SCREEN == 1) then ! both conv. and unconv. included

    if (CHO_USEABS) then ! neg. diag. might be included

      KAB = 0
      do ISYM=1,NSYM
        if (NNBSTR(ISYM,3) > 0) then

          JAB1 = IIBSTR(ISYM,3)+1
          JAB2 = JAB1+NNBSTR(ISYM,3)-1

          IAB = INDRED(JAB1,3)
          XM = abs(DIAG(IAB))
          do JAB=JAB1+1,JAB2
            IAB = INDRED(JAB,3)
            XM = max(XM,abs(DIAG(IAB)))
          end do

          if (XM > THRCOM) then  ! only if not converged
            do ISHLAB=1,NNSHL
              JAB1 = IIBSTR(ISYM,3)+IIBSTRSH(ISYM,ISHLAB,3)+1
              JAB2 = JAB1+NNBSTRSH(ISYM,ISHLAB,3)-1
              do JAB=JAB1,JAB2
                IAB = INDRED(JAB,3)
                TST = sqrt(abs(DIAG(IAB))*XM)*DAMP(2)
                if (TST > THRCOM) then
                  KAB = KAB+1
                  INDRED(KAB,2) = IAB
                  NNBSTRSH(ISYM,ISHLAB,2) = NNBSTRSH(ISYM,ISHLAB,2)+1
                end if
              end do
            end do
          end if

        end if
      end do

    else ! neg. diag. excluded

      KAB = 0
      do ISYM=1,NSYM
        if (NNBSTR(ISYM,3) > 0) then

          JAB1 = IIBSTR(ISYM,3)+1
          JAB2 = JAB1+NNBSTR(ISYM,3)-1

          IAB = INDRED(JAB1,3)
          XM = abs(DIAG(IAB))
          do JAB=JAB1+1,JAB2
            IAB = INDRED(JAB,3)
            XM = max(XM,abs(DIAG(IAB)))
          end do

          if (XM > THRCOM) then  ! only if not converged
            do ISHLAB=1,NNSHL
              JAB1 = IIBSTR(ISYM,3)+IIBSTRSH(ISYM,ISHLAB,3)+1
              JAB2 = JAB1+NNBSTRSH(ISYM,ISHLAB,3)-1
              do JAB=JAB1,JAB2
                IAB = INDRED(JAB,3)
                if (DIAG(IAB) > Zero) then ! neg=>conv
                  TST = sqrt(DIAG(IAB)*XM)*DAMP(2)
                  if (TST > THRCOM) then
                    KAB = KAB+1
                    INDRED(KAB,2) = IAB
                    NNBSTRSH(ISYM,ISHLAB,2) = NNBSTRSH(ISYM,ISHLAB,2)+1
                  end if
                end if
              end do
            end do
          end if

        end if
      end do

    end if

  else if (MODE_SCREEN == 2) then ! only unconv. included

    if (CHO_USEABS) then ! neg. diag. might be included

      KAB = 0
      do ISYM=1,NSYM
        if (NNBSTR(ISYM,3) > 0) then

          do ISHLAB=1,NNSHL
            JAB1 = IIBSTR(ISYM,3)+IIBSTRSH(ISYM,ISHLAB,3)+1
            JAB2 = JAB1+NNBSTRSH(ISYM,ISHLAB,3)-1
            do JAB=JAB1,JAB2
              IAB = INDRED(JAB,3)
              if (abs(DIAG(IAB)) > THRCOM) then
                KAB = KAB+1
                INDRED(KAB,2) = IAB
                NNBSTRSH(ISYM,ISHLAB,2) = NNBSTRSH(ISYM,ISHLAB,2)+1
              end if
            end do
          end do

        end if
      end do

    else ! neg. diag. excluded

      KAB = 0
      do ISYM=1,NSYM
        if (NNBSTR(ISYM,3) > 0) then

          do ISHLAB=1,NNSHL
            JAB1 = IIBSTR(ISYM,3)+IIBSTRSH(ISYM,ISHLAB,3)+1
            JAB2 = JAB1+NNBSTRSH(ISYM,ISHLAB,3)-1
            do JAB=JAB1,JAB2
              IAB = INDRED(JAB,3)
              if (DIAG(IAB) > THRCOM) then
                KAB = KAB+1
                INDRED(KAB,2) = IAB
                NNBSTRSH(ISYM,ISHLAB,2) = NNBSTRSH(ISYM,ISHLAB,2)+1
              end if
            end do
          end do

        end if
      end do

    end if

  else if (MODE_SCREEN == 3) then ! only 1-center unconv. incl.

    if (CHO_USEABS) then ! neg. diag. might be included

      KAB = 0
      do ISYM=1,NSYM
        if (NNBSTR(ISYM,3) > 0) then

          do ISHLAB=1,NNSHL
            call CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.true.)
            if (IATOMSHL(ISHLA) == IATOMSHL(ISHLB)) then
              JAB1 = IIBSTR(ISYM,3)+IIBSTRSH(ISYM,ISHLAB,3)+1
              JAB2 = JAB1+NNBSTRSH(ISYM,ISHLAB,3)-1
              do JAB=JAB1,JAB2
                IAB = INDRED(JAB,3)
                if (abs(DIAG(IAB)) > THRCOM) then
                  KAB = KAB+1
                  INDRED(KAB,2) = IAB
                  NNBSTRSH(ISYM,ISHLAB,2) = NNBSTRSH(ISYM,ISHLAB,2)+1
                end if
              end do
            end if
          end do

        end if
      end do

    else ! neg. diag. excluded

      KAB = 0
      do ISYM=1,NSYM
        if (NNBSTR(ISYM,3) > 0) then

          do ISHLAB=1,NNSHL
            call CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.true.)
            if (IATOMSHL(ISHLA) == IATOMSHL(ISHLB)) then
              JAB1 = IIBSTR(ISYM,3)+IIBSTRSH(ISYM,ISHLAB,3)+1
              JAB2 = JAB1+NNBSTRSH(ISYM,ISHLAB,3)-1
              do JAB=JAB1,JAB2
                IAB = INDRED(JAB,3)
                if (DIAG(IAB) > THRCOM) then
                  KAB = KAB+1
                  INDRED(KAB,2) = IAB
                  NNBSTRSH(ISYM,ISHLAB,2) = NNBSTRSH(ISYM,ISHLAB,2)+1
                end if
              end do
            end if
          end do

        end if
      end do

    end if

  else ! MODE_SCREEN out of bounds

    write(LUPRI,*) SECNAM,': MODE_SCREEN = ',MODE_SCREEN
    call CHO_QUIT('MODE_SCREEN out of bounds in '//SECNAM,103)

  end if

else ! no screening; remove zero diagonals and check convergence

  if (CHO_USEABS) then ! neg diag might be incl.

    KAB = 0
    do ISYM=1,NSYM
      if (NNBSTR(ISYM,3) > 0) then

        JAB1 = IIBSTR(ISYM,3)+1
        JAB2 = JAB1+NNBSTR(ISYM,3)-1

        IAB = INDRED(JAB1,3)
        XM = abs(DIAG(IAB))
        do JAB=JAB1+1,JAB2
          IAB = INDRED(JAB,3)
          XM = max(XM,abs(DIAG(IAB)))
        end do

        if (XM > THRCOM) then  ! only if not converged
          do ISHLAB=1,NNSHL
            JAB1 = IIBSTR(ISYM,3)+IIBSTRSH(ISYM,ISHLAB,3)+1
            JAB2 = JAB1+NNBSTRSH(ISYM,ISHLAB,3)-1
            do JAB=JAB1,JAB2
              IAB = INDRED(JAB,3)
              if (abs(DIAG(IAB)) > Zero) then
                KAB = KAB+1
                INDRED(KAB,2) = IAB
                NNBSTRSH(ISYM,ISHLAB,2) = NNBSTRSH(ISYM,ISHLAB,2)+1
              end if
            end do
          end do
        end if

      end if
    end do

  else ! neg diag excl.

    KAB = 0
    do ISYM=1,NSYM
      if (NNBSTR(ISYM,3) > 0) then

        JAB1 = IIBSTR(ISYM,3)+1
        JAB2 = JAB1+NNBSTR(ISYM,3)-1

        IAB = INDRED(JAB1,3)
        XM = abs(DIAG(IAB))
        do JAB=JAB1+1,JAB2
          IAB = INDRED(JAB,3)
          XM = max(XM,abs(DIAG(IAB)))
        end do

        if (XM > THRCOM) then  ! only if not converged
          do ISHLAB=1,NNSHL
            JAB1 = IIBSTR(ISYM,3)+IIBSTRSH(ISYM,ISHLAB,3)+1
            JAB2 = JAB1+NNBSTRSH(ISYM,ISHLAB,3)-1
            do JAB=JAB1,JAB2
              IAB = INDRED(JAB,3)
              if (DIAG(IAB) > Zero) then
                KAB = KAB+1
                INDRED(KAB,2) = IAB
                NNBSTRSH(ISYM,ISHLAB,2) = NNBSTRSH(ISYM,ISHLAB,2)+1
              end if
            end do
          end do
        end if

      end if
    end do

  end if

end if

! Set remaining index arrays.
! ---------------------------

call CHO_SETREDIND(2)

! Debug print.
! ------------

if (CHO_TRCNEG) then
  write(LUPRI,*) SECNAM,': checking for negative diagonals in next reduced set:'
  NNEG = 0
  do ISYM=1,NSYM
    JAB1 = IIBSTR(ISYM,2)+1
    JAB2 = JAB1+NNBSTR(ISYM,2)-1
    INEG = 0
    do JAB=JAB1,JAB2
      IAB = INDRED(JAB,2)
      if (DIAG(IAB) < Zero) INEG = INEG+1
    end do
    NNEG = NNEG+INEG
    write(LUPRI,*) SECNAM,': #negative in symmetry ',ISYM,': ',INEG
  end do
  write(LUPRI,*) SECNAM,': total #negative: ',NNEG
end if

end subroutine CHO_SETRED
