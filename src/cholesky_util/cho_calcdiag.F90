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

subroutine CHO_CALCDIAG(BUF,IBUF,LENBUF,SCR,LENSCR,NDUMP)
!
! Purpose: shell-driven calculation of the integral diagonal and
!          setup of the first reduced set.

use Index_Functions, only: iTri, nTri_Elem
use Symmetry_Info, only: Mul
use Cholesky, only: CHO_NO2CENTER, Cho_PreScreen, CHO_USEABS, Damp, DiaMax, iAtomShl, iBasSh, IPRINT, iSP2F, lBuf, LuPri, LuScr, &
                    Mx2Sh, MySP, n_MySP, NBAST, nBasSh, nBstSh, nnBstRSh, nnShl, nShell, nSym, SCDIAG, ShA, ShB, ThrCom, ThrDiag, &
                    XlDiag
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LENBUF, LENSCR
real(kind=wp), intent(out) :: BUF(LENBUF), SCR(LENSCR)
integer(kind=iwp), intent(out) :: IBUF(4,LENBUF), NDUMP
integer(kind=iwp) :: I, I_MYSP, IA, IAA, IAB, IB, IBB, ICOUNT, IDUMP, IOPT, ISAB, ISHLA, ISHLAB, ISHLB, ISYM, ISYMA, ISYMAB, &
                     ISYMB, IUNIT, JUNIT, L, LENGTH, LINTD, ll, LSCR, n_NegCalcDiag, n_NegCalcDiag_local, NIATOMSHL, NUMA, NUMAB, &
                     NUMB
real(kind=wp) :: DEL1, DIAAB, DIAGAB, SAVD, SCRMAX(8), XLDIA, XMDIA, XNCD(1), XXX
real(kind=wp), allocatable :: NegCalcDiag(:)
integer(kind=iwp), parameter :: INFO_DEBUG = 4, INFO_INSANE = 10
character(len=*), parameter :: SECNAM = 'CHO_CALCDIAG'
integer(kind=iwp), external :: CHO_ISAOSH

! Check dimensions.
! -----------------

if (LENBUF < LBUF) then
  write(LUPRI,'(//,1X,A,A)') SECNAM,': LENBUF >= LBUF required!'
  write(LUPRI,'(1X,A,I10)') 'LENBUF = ',LENBUF
  write(LUPRI,'(1X,A,I10,/)') 'LBUF   = ',LBUF
  call CHO_QUIT('Buffer error in '//SECNAM,102)
end if
LSCR = MX2SH
if (LENSCR < LSCR) then
  write(LUPRI,'(//,1X,A,A)') SECNAM,': LENSCR >= MX2SH required!'
  write(LUPRI,'(1X,A,I10)') 'LENSCR = ',LENSCR
  write(LUPRI,'(1X,A,I10,/)') 'MX2SH  = ',LSCR
  call CHO_QUIT('Scratch space error in '//SECNAM,102)
end if

! Open scratch files.
! -------------------

IUNIT = -1
call CHO_OPEN(IUNIT,'_CHO_DIASCR2')
JUNIT = -1
call CHO_OPEN(JUNIT,'_CHO_DIASCR1')
rewind(IUNIT)
rewind(JUNIT)

! Make JUNIT available outside this routine.
! ------------------------------------------

LUSCR = JUNIT

! Initialize abs. max. diag. array.
! ---------------------------------

DIAMAX(1:NSYM) = Zero

! Allocate array for storing 10 most negative diagonals
! (there should be none, of course, but they do show up)
! ------------------------------------------------------

call mma_allocate(NegCalcDiag,10,Label='NegCalcDiag')
NegCalcDiag(:) = Zero
n_NegCalcDiag = 0

! Calculate diagonal in loop over shell-pairs.
! CHO_NO2CENTER on: skip all 2-center diagonals.
! ----------------------------------------------

NIATOMSHL = 0
if (allocated(IATOMSHL)) NIATOMSHL = size(IATOMSHL)
if (CHO_NO2CENTER .and. (NIATOMSHL < NSHELL)) call CHO_QUIT(SECNAM//': iAtomShl not allocated correctly!',103)

XLDIAG = Zero
ICOUNT = 0
NDUMP = 0
N_MYSP = 0
IOPT = 2
call CHO_P_DISTRIB_SP(IOPT,MYSP,N_MYSP)
call mma_maxDBLE(LINTD)
call XSETMEM_INTS(LINTD) ! set memory for seward
do I_MYSP=1,N_MYSP

  ISAB = MYSP(I_MYSP)

  ISHLAB = ISP2F(ISAB)
  call CHO_INVPCK(ISHLAB,ISHLA,ISHLB,.true.)

  if (CHO_NO2CENTER) then
    if (IATOMSHL(ISHLA) /= IATOMSHL(ISHLB)) cycle
  end if

  NUMA = NBSTSH(ISHLA)
  NUMB = NBSTSH(ISHLB)
  if (ISHLA == ISHLB) then
    NUMAB = nTri_Elem(NUMA)
  else
    NUMAB = NUMA*NUMB
  end if

  SHA = ISHLA
  SHB = ISHLB

  call CHO_MCA_DIAGINT(ISHLA,ISHLB,SCR,NUMAB)

  if (IPRINT >= INFO_INSANE) then
    if ((ISHLA == 1) .and. (ISHLB == 1)) then
      if (CHO_PRESCREEN) then
        call CHO_HEAD(SECNAM//': Prescreened Diagonal','=',80,LUPRI)
      else
        call CHO_HEAD(SECNAM//': Unscreened Diagonal','=',80,LUPRI)
      end if
    end if
    write(LUPRI,'(/,2X,A,I10,1X,I10,1X,I10)') 'Diagonal shell block A,B,AB = ',ISHLA,ISHLB,ITRI(ISHLA,ISHLB)
    if (ISHLA == ISHLB) then
      call CHO_OUTPAK(SCR,NUMA,1,LUPRI)
    else
      call CHO_OUTPUT(SCR,1,NUMA,1,NUMB,NUMA,NUMB,1,LUPRI)
    end if
  end if

  if (ISHLA == ISHLB) then
    do IA=1,NUMA
      ISYMA = CHO_ISAOSH(IA,ISHLA)
      do IB=1,IA
        ISYMB = CHO_ISAOSH(IB,ISHLB)
        ISYMAB = MUL(ISYMB,ISYMA)
        IAB = ITRI(IA,IB)
        DIAAB = SCR(IAB)
        if (DIAAB < Zero) then
          n_NegCalcDiag = n_NegCalcDiag+1
          call UpdateMostNegative(size(NegCalcDiag),NegCalcDiag,DIAAB)
        end if
        if (DIAAB > THRDIAG) then
          DIAMAX(ISYMAB) = max(DIAMAX(ISYMAB),DIAAB)
          ICOUNT = ICOUNT+1
          BUF(ICOUNT) = SCR(IAB)
          IBUF(1,ICOUNT) = ISAB
          IBUF(2,ICOUNT) = IAB
          IBUF(3,ICOUNT) = ISYMAB
          IBUF(4,ICOUNT) = IAB
          if (ICOUNT == LBUF) then
            call CHO_WRBUF(LBUF,BUF,IBUF,LBUF,IUNIT)
            XLDIAG = XLDIAG+real(LBUF,kind=wp)
            ICOUNT = 0
            NDUMP = NDUMP+1
          end if
        end if
      end do
    end do
  else
    do ISYMB=1,NSYM
      do IBB=1,NBASSH(ISYMB,ISHLB)
        IB = IBASSH(ISYMB,ISHLB)+IBB
        do ISYMA=1,NSYM
          do IAA=1,NBASSH(ISYMA,ISHLA)
            IA = IBASSH(ISYMA,ISHLA)+IAA
            ISYMAB = MUL(ISYMA,ISYMB)
            IAB = NUMA*(IB-1)+IA
            DIAAB = SCR(IAB)
            if (DIAAB < Zero) then
              n_NegCalcDiag = n_NegCalcDiag+1
              call UpdateMostNegative(size(NegCalcDiag),NegCalcDiag,DIAAB)
            end if
            if (DIAAB > THRDIAG) then
              DIAMAX(ISYMAB) = max(DIAMAX(ISYMAB),DIAAB)
              ICOUNT = ICOUNT+1
              BUF(ICOUNT) = SCR(IAB)
              IBUF(1,ICOUNT) = ISAB
              IBUF(2,ICOUNT) = IAB
              IBUF(3,ICOUNT) = ISYMAB
              IBUF(4,ICOUNT) = IAB
              if (ICOUNT == LBUF) then
                call CHO_WRBUF(LBUF,BUF,IBUF,LBUF,IUNIT)
                XLDIAG = XLDIAG+real(LBUF,kind=wp)
                ICOUNT = 0
                NDUMP = NDUMP+1
              end if
            end if
          end do
        end do
      end do
    end do
  end if

end do
if (ICOUNT > 0) then ! flush buffer
  if (ICOUNT > LBUF) call CHO_QUIT('Logical error in '//SECNAM,103)
  call CHO_WRBUF(ICOUNT,BUF,IBUF,LBUF,IUNIT)
  XLDIAG = XLDIAG+real(ICOUNT,kind=wp)
  ICOUNT = 0
  NDUMP = NDUMP+1
end if
call XRLSMEM_INTS ! release memory (seward)
call CHO_GADGOP(DIAMAX,NSYM,'max') ! sync abs. max. diag.
n_NegCalcDiag_local = n_NegCalcDiag
XNCD(1) = real(n_NegCalcDiag,kind=wp)
call CHO_GADGOP(XNCD,1,'+')
n_NegCalcDiag = int(XNCD(1))
if (n_NegCalcDiag > 0) then
  call WarningMessage(1,'WARNING: negative integral diagonal elements computed')
  write(LuPri,'(3X,A)') 'All negative integral diagonal elements have been removed (zeroed) - they are considered irrelevant!'
  write(LuPri,'(3X,A,I10)') 'Number of negative elements computed:   ',n_NegCalcDiag
  write(LuPri,'(3X,A,I10)') 'Number of negative elements (this node):',n_NegCalcDiag_local
  if (n_NegCalcDiag_local > 0) then
    ll = min(n_NegCalcDiag_local,size(NegCalcDiag))
    write(LuPri,'(I5,A)') ll,' most negative elements (this node):'
    write(LuPri,'(10ES12.4)') (NegCalcDiag(i),i=1,ll)
  end if
  call CHO_GADGOP(NegCalcDiag,1,'min')
  write(LuPri,'(3X,A,ES12.4)') 'Most negative element overall: ',NegCalcDiag(1)
end if
call mma_deallocate(NegCalcDiag)

if (IPRINT >= INFO_DEBUG) then
  call CHO_HEAD(SECNAM//': Diagonal Info','=',80,LUPRI)
  XXX = real(NBAST,kind=wp)
  XMDIA = XXX*(XXX+One)*Half
  XLDIA = XLDIAG
  SAVD = 1.0e2_wp*(XMDIA-XLDIA)/XMDIA
  write(LUPRI,'(/,2X,A,ES15.6)') 'Screening threshold for initial diagonal: ',THRDIAG
  write(LUPRI,'(2X,A,F15.1,/,2X,A,F15.1)') 'Dimension of unscreened initial diagonal: ',XMDIA, &
                                           'Dimension of   screened initial diagonal: ',XLDIA
  write(LUPRI,'(2X,A,7X,F8.3,A)') 'Saving from screening                   : ',SAVD,'%'
  do ISYM=1,NSYM
    write(LUPRI,'(2X,A,I2,12X,A,ES15.6)') 'Maximum diagonal, symmetry',ISYM,': ',DIAMAX(ISYM)
  end do
  write(LUPRI,'(2X,A,5X,I10)') 'Number of negative diagonals computed   : ',n_NegCalcDiag
end if

! Read through the file to get first reduced set.
! -----------------------------------------------

nnBstRSh(:,:,1) = 0

rewind(IUNIT)
rewind(JUNIT)
if (SCDIAG) then ! screen diagonal
  DEL1 = THRCOM*THRCOM/DAMP(1)
  do ISYM=1,NSYM
    if (abs(DIAMAX(ISYM)) > Zero) then
      SCRMAX(ISYM) = DEL1/DIAMAX(ISYM)
    else
      SCRMAX(ISYM) = 1.0e15_wp
    end if
  end do
  if (CHO_USEABS) then
    do IDUMP=1,NDUMP
      call CHO_RDBUF(LENGTH,BUF,IBUF,LBUF,IUNIT)
      if (IDUMP == NDUMP) call CHO_CLOSE(IUNIT,'DELETE')
      do L=1,LENGTH
        DIAGAB = BUF(L)
        ISYMAB = IBUF(3,L)
        if (abs(DIAGAB) < SCRMAX(ISYMAB)) then
          BUF(L) = Zero
          IBUF(2,L) = -1
        else
          ISHLAB = IBUF(1,L)
          NNBSTRSH(ISYMAB,ISHLAB,1) = NNBSTRSH(ISYMAB,ISHLAB,1)+1
          IBUF(2,L) = NNBSTRSH(ISYMAB,ISHLAB,1)
        end if
      end do
      call CHO_WRBUF(LENGTH,BUF,IBUF,LBUF,JUNIT)
    end do
  else
    do IDUMP=1,NDUMP
      call CHO_RDBUF(LENGTH,BUF,IBUF,LBUF,IUNIT)
      if (IDUMP == NDUMP) call CHO_CLOSE(IUNIT,'DELETE')
      do L=1,LENGTH
        DIAGAB = BUF(L)
        ISYMAB = IBUF(3,L)
        if (DIAGAB < SCRMAX(ISYMAB)) then
          BUF(L) = Zero
          IBUF(2,L) = -1
        else
          ISHLAB = IBUF(1,L)
          NNBSTRSH(ISYMAB,ISHLAB,1) = NNBSTRSH(ISYMAB,ISHLAB,1)+1
          IBUF(2,L) = NNBSTRSH(ISYMAB,ISHLAB,1)
        end if
      end do
      call CHO_WRBUF(LENGTH,BUF,IBUF,LBUF,JUNIT)
    end do
  end if
else ! no screening at all
  do IDUMP=1,NDUMP
    call CHO_RDBUF(LENGTH,BUF,IBUF,LBUF,IUNIT)
    if (IDUMP == NDUMP) call CHO_CLOSE(IUNIT,'DELETE')
    do L=1,LENGTH
      ISHLAB = IBUF(1,L)
      ISYMAB = IBUF(3,L)
      NNBSTRSH(ISYMAB,ISHLAB,1) = NNBSTRSH(ISYMAB,ISHLAB,1)+1
      IBUF(2,L) = NNBSTRSH(ISYMAB,ISHLAB,1)
    end do
    call CHO_WRBUF(LENGTH,BUF,IBUF,LBUF,JUNIT)
  end do
end if

call CHO_GAIGOP(NNBSTRSH(:,:,1),NSYM*NNSHL,'+') ! sync
call CHO_SETREDIND(1)

end subroutine CHO_CALCDIAG
