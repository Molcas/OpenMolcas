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

subroutine CHO_CC_drv(rc,CMO)
!***********************************************************************
!
! a,b,g,d:  AO-index
! p,q,r,s:  MO-indices belonging to (probably frozen excluded ?)
!
!***********************************************************************

use Cholesky, only: InfVec, nBas, nDimRS, nSym, NumCho, timings
use Symmetry_Info, only: Mul
use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type, SBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: rc
type(DSBA_Type), intent(in) :: CMO
integer(kind=iwp) :: i, iBatch, idisk, iE, iLoc, irc, IREDC, iS, iSwap, iSymb, iSymp, IVEC2, iVrs, JNUM, JRED, JRED1, JRED2, jSym, &
                     JVC, JVEC, k, l, LREAD, LunChVF, LWORK, mTTvec, mTvec, MUSED, mvec, NAp, NAq, nBatch, nPorb(8), nRS, NUMV, &
                     nVec, nVrs
real(kind=wp) :: TCM1, TCM2, TCM3, TCM4, TCR1, TCR2, TCR3, TCR4, tmotr1(2), tmotr2(2), TOTCPU, TOTCPU1, TOTCPU2, TOTWALL, &
                 TOTWALL1, TOTWALL2, tread(2), TWM1, TWM2, TWM3, TWM4, TWR1, TWR2, TWR3, TWR4
logical(kind=iwp) :: DoRead
character(len=50) :: CFmt
type(SBA_Type) :: Laq(1)
real(kind=wp), allocatable :: Lrs(:,:)
real(kind=wp), allocatable, target :: Lpq(:)
real(kind=wp), pointer :: pLpq(:,:,:) => null()
character(len=*), parameter :: SECNAM = 'CHO_CC_drv'
integer(kind=iwp), external :: isFreeUnit

LunChVF = 80
LunChVF = isfreeunit(LunChVF)
call DaName_mf_wa(LunChVF,'CD1tmp')
idisk = 1
DoRead = .false.
IREDC = -1  ! unknown reduced set in core

iSwap = 0  ! Lpb,J are returned by cho_x_getVtra

call CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

! 1 --> CPU   2 --> Wall
tread(:) = zero   !time read/write vectors
tmotr1(:) = zero  !time 1st MO half-transf.
tmotr2(:) = zero  !time 2nd MO half-transf.

! --- Define MOs used in CC
! -----------------------------------
do i=1,nSym
  nPorb(i) = size(CMO%SB(i)%A2,1)
end do

! ======================================================================

! --- Various offsets & pointers
! ------------------------------

iLoc = 3 ! use scratch location in reduced index arrays

! *************** BIG LOOP OVER VECTORS SYMMETRY ***********************

do jSym=1,nSym

  if (NumCho(jSym) < 1) cycle

  ! --------------------------------------------------------------------

  ! **************     MEMORY MANAGEMENT SECTION    ********************
  !---------------------------------------------------------------------
  ! --- compute memory needed to store at least 1 vector of JSYM
  ! --- and do all the subsequent calculations
  ! --------------------------------------------------------------------
  mTvec = 0  ! mem for storing half-transformed vec Laq,J
  mTTvec = 0  ! mem for storing transformed vec Lpq,J

  do l=1,nSym
    k = Mul(l,JSYM)
    mTvec = mTvec+nPorb(l)*nBas(k)
    mTTvec = max(mTTvec,nPorb(l)*nPorb(k))
  end do

  mvec = mTvec+mTTvec

  ! --------------------------------------------------------------------
  ! --------------------------------------------------------------------

  JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
  JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last vec

  do JRED=JRED1,JRED2

    call Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)

    if (nVrs == 0) cycle  ! no vectors in that (jred,jsym)

    if (nVrs < 0) then
      write(u6,*) SECNAM//': Cho_X_nVecRS returned nVrs<0. STOP!'
      call abend()
    end if

    call Cho_X_SetRed(irc,iLoc,JRED) !set index arrays at iLoc
    if (irc /= 0) then
      write(u6,*) SECNAM//'cho_X_setred non-zero return code. rc= ',irc
      call abend()
    end if

    IREDC = JRED

    nRS = nDimRS(JSYM,JRED)

    call mma_maxDBLE(LWORK)

    nVec = min(LWORK/(nRS+mvec),nVrs)

    if (nVec < 1) then
      write(u6,*) SECNAM//': Insufficient memory for batch'
      write(u6,*) 'LWORK= ',LWORK
      write(u6,*) 'min. mem. need= ',nRS+mTvec
      write(u6,*) 'reading ',nRS,' and transforming to ',mvec
      write(u6,*) 'of jsym= ',jsym,' and JRED= ',JRED
      rc = 33
      call Abend()
      nBatch = -9999  ! dummy assignment
    end if

    LREAD = nRS*nVec

    call mma_allocate(Lrs,nRS,nVec,Label='Lrs')
    call mma_allocate(Lpq,mTTVec*nVec,Label='Lpq')

    ! --- BATCH over the vectors ----------------------------

    nBatch = (nVrs-1)/nVec+1

    do iBatch=1,nBatch

      if (iBatch == nBatch) then
        JNUM = nVrs-nVec*(nBatch-1)
      else
        JNUM = nVec
      end if
      call Allocate_DT(Laq(1),nPorb,nBas,JNUM,jSym,nSym,iSwap)

      JVEC = nVec*(iBatch-1)+iVrs
      IVEC2 = JVEC-1+JNUM

      call CWTIME(TCR1,TWR1)

      call CHO_VECRD(Lrs,LREAD,JVEC,IVEC2,JSYM,NUMV,IREDC,MUSED)

      if ((NUMV <= 0) .or. (NUMV /= JNUM)) then
        rc = 77
        return
      end if

      call CWTIME(TCR2,TWR2)
      tread(1) = tread(1)+(TCR2-TCR1)
      tread(2) = tread(2)+(TWR2-TWR1)

      ! ----------------------------------------------------------------
      ! --- First half MO transformation  Lpb,J = sum_a  C(p,a) * Lab,J
      ! ----------------------------------------------------------------

      call CWTIME(TCM1,TWM1)

      call CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,JSYM,iSwap,IREDC,1,1,[CMO],Laq,DoRead)

      if (irc /= 0) then
        rc = irc
        return
      end if

      call CWTIME(TCM2,TWM2)
      tmotr1(1) = tmotr1(1)+(TCM2-TCM1)
      tmotr1(2) = tmotr1(2)+(TWM2-TWM1)

      ! ----------------------------------------------------------------
      ! --- 2nd half of MO transformation  Lpq,J = sum_b  Lpb,J * C(q,b)
      ! ----------------------------------------------------------------
      do iSymb=1,nSym

        iSymp = Mul(JSYM,iSymb)
        NAp = nPorb(iSymp)
        NAq = nPorb(iSymb) ! iSymb=iSymq
        iS = 1
        iE = NAp*NAq*JNUM

        pLpq(1:NAp,1:NAq,1:JNUM) => Lpq(iS:iE)

        call CWTIME(TCM3,TWM3)

        if (NAp*NAq /= 0) then

          do JVC=1,JNUM

            call DGEMM_('N','T',NAp,NAq,nBas(iSymb),One,Laq(1)%SB(iSymb)%A3(:,:,JVC),NAp,CMO%SB(iSymb)%A2,NAq,Zero,pLpq(:,:,JVC), &
                        NAp)

          end do

        end if

        call CWTIME(TCM4,TWM4)
        tmotr2(1) = tmotr2(1)+(TCM4-TCM3)
        tmotr2(2) = tmotr2(2)+(TWM4-TWM3)

        ! if u need to compute fock matrix elements this should be done probably here
        !     I can help you with that

        call CWTIME(TCR3,TWR3)
        ! --- WRITE transformed vectors to disk (each Jsym on a separate file!)

        call ddafile(LunChVF,1,Lpq,NAp*NAq*JNUM,idisk)

        ! --- remember that this is inside a batch over J, the vector index

        call CWTIME(TCR4,TWR4)
        tread(1) = tread(1)+(TCR4-TCR3)
        tread(2) = tread(2)+(TWR4-TWR3)

        pLpq => null()

      end do

      ! --------------------------------------------------------------------
      ! --------------------------------------------------------------------
      call Deallocate_DT(Laq(1))

    end do  ! end batch loop

    ! --- free memory
    call mma_deallocate(Lpq)
    call mma_deallocate(Lrs)

  end do   ! loop over red sets

end do   !loop over JSYM
call daclos(LunChVF)

call CWTIME(TOTCPU2,TOTWALL2)
TOTCPU = TOTCPU2-TOTCPU1
TOTWALL = TOTWALL2-TOTWALL1

!---- Write out timing information
if (timings) then

  CFmt = '(2x,A)'
  write(u6,*)
  write(u6,CFmt) 'Cholesky-CC timing from '//SECNAM
  write(u6,CFmt) '----------------------------------------'
  write(u6,*)
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,CFmt) 'MO transf. Cholesky vectors     CPU       WALL   '
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'

  write(u6,'(2x,A26,2f10.2)') 'READ/WRITE VECTORS                        ',tread(1),tread(2)
  write(u6,'(2x,A26,2f10.2)') '1st half-transf.                          ',tmotr1(1),tmotr1(2)
  write(u6,'(2x,A26,2f10.2)') '2nd half-transf.                          ',tmotr2(1),tmotr2(2)
  write(u6,*)
  write(u6,'(2x,A26,2f10.2)') 'TOTAL                                     ',TOTCPU,TOTWALL
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,*)

end if

rc = 0

return

end subroutine CHO_CC_drv
