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

subroutine CHO_TR_drv(rc,nIsh,nAsh,nSsh,Porb,BName,Do_int,ihdf5,Xint,lXint)
!*********************************************************************
!  a,b,g,d:  AO-index
!  p,q,r,s:  MO-indices belonging to all (fro and del excluded)
!*********************************************************************

#ifdef _HDF5_QCM_
use hdf5_utils, only: file_id, hdf5_close_cholesky, hdf5_init_wr_cholesky, hdf5_write_cholesky, HID_T
#endif
use ChoArr, only: nDimRS
use ChoSwp, only: InfVec
use Symmetry_Info, only: Mul
use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type, SBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp), intent(in) :: nIsh(*), nAsh(*), nSsh(*), ihdf5, lXint
type(DSBA_Type), intent(in) :: Porb(1)
character(len=6), intent(in) :: BName
logical(kind=iwp), intent(in) :: Do_int
real(kind=wp), intent(out) :: Xint(0:lXint-1)
integer(kind=iwp) :: i, iBatch, idisk, iLoc, iOffB(8), ipq, irc, IREDC, iSwap, iSymb, iSymp, IVEC2, iVrs, JNUM, JRED, JRED1, &
                     JRED2, jSym, JVC, JVEC, k, kMOs, kt, l, LREAD, LunChVF(8), LWORK, kOff(8), Mpq, mTTvec, mTvec, MUSED, mvec, &
                     NAp, nApq, nAq, nBatch, nMOs, nOB(8), nPorb(8), nRS, NUMV, nVec, nVrs
#ifdef _HDF5_QCM_
integer(kind=HID_T) :: choset_id, space_id
#endif
real(kind=wp) :: TCM1, TCM2, TCM3, TCM4, TCR1, TCR2, TCR3, TCR4, TOTCPU, TOTCPU1, TOTCPU2, TOTWALL, TOTWALL1, TOTWALL2, TWM1, &
                 TWM2, TWM3, TWM4, TWR1, TWR2, TWR3, TWR4, tmotr1(2), tmotr2(2), tread(2)
character(len=50) :: CFmt
character(len=7) :: Fnam
type(SBA_Type), target :: ChoT(1)
real(kind=wp), allocatable :: Lpq(:,:), Lpq_J(:), Lrs(:)
logical(kind=iwp), parameter :: DoRead = .false.
character(len=10), parameter :: SECNAM = 'CHO_TR_drv'
integer(kind=iwp), external :: IsFreeUnit
real(kind=wp), external :: ddot_
#include "chotime.fh"
#include "chotraw.fh"
#include "cholesky.fh"
#include "choorb.fh"

#ifndef _HDF5_QCM_
#include "macros.fh"
unused_var(ihdf5)
#endif

#ifdef _HDF5_QCM_
! Leon 13.6.2017: Avoid opening a regular file if HDF5 is used
if (ihdf5 /= 1) then
#endif
  do i=1,nSym
    LunChVF(i) = 80
    LunChVF(i) = isfreeunit(LunChVF(i))
    write(Fnam,'(A6,I1)') BName,i
    call DaName_mf_wa(LunChVF(i),Fnam)
  end do
#ifdef _HDF5_QCM_
end if
#endif

IREDC = -1  ! unknown reduced set in core

call CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

! 1 --> CPU   2 --> Wall
tread(:) = Zero   !time read/write vectors
tmotr1(:) = Zero  !time 1st MO half-transf.
tmotr2(:) = Zero  !time 2nd MO half-transf.

if (Do_int) call Fzero(Xint(0),lXint)

! Define MO space used
!---------------------
do i=1,nSym
  nPorb(i) = nIsh(i)+nAsh(i)+nSsh(i)
end do

! ======================================================================

iLoc = 3 ! use scratch location in reduced index arrays

!****************** BIG LOOP OVER VECTORS SYMMETRY *********************

Mpq = 0

do jSym=1,nSym

! Init the HDF5 section
#ifdef _HDF5_QCM_
  if (ihdf5 == 1) then
    ! max size of the cholesky vector for now is the
    ! max(nPorb)*(max(nPorb)+1)/2
    ! probably this can be chosen more efficiently,
    ! but would matter only if we use symmetry
    call hdf5_init_wr_cholesky(file_id(1),JSym,maxval(nPorb(1:nSym))*(maxval(nPorb(1:nSym))+1)/2,NumCho(JSym),choset_id,space_id)
  end if
#endif
  if (NumCho(jSym) < 1) cycle

  ! Set up the skipping flags + some initializations
  !-------------------------------------------------
  do i=1,nSym
    k = Mul(i,JSYM)
    if (i < k) then
      kOff(i) = Mpq
      nOB(i) = nPorb(i)*nPorb(k)*NumCho(jSym)
      Mpq = Mpq+nPorb(i)*nPorb(k)
    else if (k == i) then
      kOff(i) = Mpq
      nOB(i) = nPorb(i)*(nPorb(i)+1)/2*NumCho(jSym)
      Mpq = Mpq+nPorb(i)*(nPorb(i)+1)/2
    else
      nOB(i) = 0
    end if
    iOffB(i) = 0
  end do

  do i=2,nSym
    iOffB(i) = iOffB(i-1)+nOB(i-1)
  end do
  do i=1,nSym
    if (nOB(i) == 0) then
      k = Mul(i,JSYM)
      iOffB(i) = iOffB(k)
      kOff(i) = kOff(k)
    end if
  end do

  !********************* MEMORY MANAGEMENT SECTION *********************
  !---------------------------------------------------------------------
  ! compute memory needed to store at least 1 vector of JSYM
  ! and do all the subsequent calculations
  !---------------------------------------------------------------------
  mTvec = 0  ! mem for storing half-transformed vec Laq,J
  mTTvec = 0  ! mem for storing transformed vec Lpq,J

  do l=1,nSym
    k = Mul(l,JSYM)
    mTvec = mTvec+nPorb(k)*nBas(l)
    mTTvec = max(mTTvec,nPorb(k)*nPorb(l))
  end do

  mvec = mTvec+mTTvec

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------

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
      write(u6,*) SECNAM//': Cho_X_SetRed non-zero return code. rc= ',irc
      call abend()
    end if

    IREDC = JRED

    nRS = nDimRS(JSYM,JRED)

    call mma_maxDBLE(LWORK)

    nVec = min(LWORK/(nRS+mvec+1),nVrs)

    if (nVec < 1) then
      write(u6,*) SECNAM//': Insufficient memory for batch'
      write(u6,*) 'LWORK= ',LWORK
      write(u6,*) 'Min. mem. need= ',nRS+mvec+1
      write(u6,*) 'Reading ',nRS,' and then MO-transform.'
      write(u6,*) 'In jsym= ',jsym,' and JRED= ',JRED
      rc = 33
      call Abend()
      nBatch = -9999  ! dummy assignment
    end if

    LREAD = nRS*nVec

    call mma_allocate(Lrs,LREAD,Label='Lrs')
    call mma_allocate(Lpq_J,nVec,Label='Lpq_j')

    iSwap = 0  ! Lpb,J are returned by cho_x_getVtra
    call Allocate_DT(ChoT(1),nPorb,nBas,nVec,JSYM,nSym,iSwap)
    ChoT(1)%A0(:) = Zero

    ! BATCH over the vectors

    nBatch = (nVrs-1)/nVec+1

    do iBatch=1,nBatch

      if (iBatch == nBatch) then
        JNUM = nVrs-nVec*(nBatch-1)
      else
        JNUM = nVec
      end if

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

      !-----------------------------------------------------------------
      ! First half MO transformation  Lpb,J = sum_a  C(p,a) * Lab,J
      !-----------------------------------------------------------------

      call CWTIME(TCM1,TWM1)

      kMOs = 1
      nMOs = 1

      call CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,jSym,iSwap,IREDC,nMOs,kMOs,POrb,ChoT(1),DoRead)

      if (irc /= 0) then
        rc = irc
        return
      end if

      call CWTIME(TCM2,TWM2)
      tmotr1(1) = tmotr1(1)+(TCM2-TCM1)
      tmotr1(2) = tmotr1(2)+(TWM2-TWM1)

      !-----------------------------------------------------------------
      ! 2nd half of MO transformation  Lpq,J = sum_b  Lpb,J * C(q,b)
      !-----------------------------------------------------------------

      if (JSYM == 1) then   !  Lpq,J in LT-storage

        do iSymb=1,nSym

          NAp = nPorb(iSymb)
          NApq = NAp*(NAp+1)/2

          call CWTIME(TCM3,TWM3)

          if (NApq == 0) cycle

          call mma_allocate(Lpq,NApq,JNUM,Label='Lpq')

          do JVC=1,JNUM

            call DGEMM_Tri('N','T',NAp,NAp,nBas(iSymb),One,ChoT(1)%SB(iSymb)%A3(:,:,JVC),NAp,Porb(1)%SB(iSymb)%A2,NAp,Zero, &
                           Lpq(:,jVC),NAp)

          end do

          call CWTIME(TCM4,TWM4)
          tmotr2(1) = tmotr2(1)+(TCM4-TCM3)
          tmotr2(2) = tmotr2(2)+(TWM4-TWM3)

          call CWTIME(TCR3,TWR3)

          if (tv2disk == 'PQK') then
#           ifdef _HDF5_QCM_
            if (ihdf5 /= 1) then
#           endif
              call ddafile(LunChVF(jSym),1,Lpq,NApq*JNUM,iOffB(iSymb))
#           ifdef _HDF5_QCM_
            else
              ! this should never happen, this case should be caught in motra.f
              write(u6,*) ' Writing of Cholesky vectors in HDF5 format as (pq,k) is not supported.'
              call Abend()
            end if
#           endif
            if (Do_int) then
              do ipq=1,NApq
                kt = kOff(iSymb)+ipq-1
                Xint(kt) = Xint(kt)+ddot_(JNUM,Lpq(ipq,:),NApq,Lpq(ipq,:),NApq)
              end do
            end if
          else
            do ipq=1,NApq
              Lpq_J(1:JNUM) = Lpq(ipq,1:JNUM)
              if (Do_int) then
                kt = kOff(iSymb)+ipq-1
                Xint(kt) = Xint(kt)+ddot_(JNUM,Lpq_J,1,Lpq_J,1)
              end if
              idisk = iOffB(iSymb)+NumCho(jSym)*(ipq-1)
#           ifdef _HDF5_QCM_
              ! Leon 13.6.2017: Do not write Cholesky vectors to the regular file
              ! if the hdf5 file is written. It becomes counterproductive to write
              ! the same content twice for large basis sets
              if (ihdf5 /= 1) then
#             endif
                call ddafile(LunChVF(jSym),1,Lpq_J,JNUM,idisk)

#             ifdef _HDF5_QCM_
              else
                ! Write the transformed Cholesky batch to the hdf5 dataset
                ! The ordering in HDF5 is in column-major order, corresponding to
                ! the 'Kpq' storage
                ! This way all the elements needed to compute one integral can be
                ! read with one read operation.

                ! TODO: eventually row-major order storage + chunked dataset might
                ! improve the performance -- but probably it's irrelevant.
                ! Leon 22.4.2016 -- modified the write_cholesky call below to
                ! account for multiple reduced sets
                call hdf5_write_cholesky(choset_id,space_id,ipq-1,nVec*(iBatch-1)+iVrs-1,JNUM,Lpq_J)
              end if
#             endif

            end do
            iOffB(iSymb) = iOffB(iSymb)+JNUM
          end if

          call mma_deallocate(Lpq)

          call CWTIME(TCR4,TWR4)
          tread(1) = tread(1)+(TCR4-TCR3)
          tread(2) = tread(2)+(TWR4-TWR3)

        end do

      else

        do iSymb=1,nSym

          iSymp = Mul(JSYM,iSymb)
          NAp = nPorb(iSymp)
          NAq = nPorb(iSymb) ! iSymb=iSymq
          NApq = NAp*NAq

          call CWTIME(TCM3,TWM3)

          if (NApq == 0) cycle

          call mma_allocate(Lpq,NApq,JNUM,Label='Lpq')

          if (iSymp < iSymb) then
            do JVC=1,JNUM

              call DGEMM_('N','T',NAp,NAq,nBas(iSymb),One,ChoT(1)%SB(iSymp)%A3(:,:,JVC),NAp,Porb(1)%SB(iSymb)%A2,NAq,Zero, &
                          Lpq(:,JVC),NAp)

            end do
          else
            Lpq(:,:) = Zero
          end if

          call CWTIME(TCM4,TWM4)
          tmotr2(1) = tmotr2(1)+(TCM4-TCM3)
          tmotr2(2) = tmotr2(2)+(TWM4-TWM3)

          call CWTIME(TCR3,TWR3)

          if (iSymp < iSymb) then

            if (tv2disk == 'PQK') then
              call ddafile(LunChVF(jSym),1,Lpq,NApq*JNUM,iOffB(iSymp))
              if (Do_int) then
                do ipq=1,NApq
                  kt = kOff(iSymp)+ipq-1
                  Xint(kt) = Xint(kt)+ddot_(JNUM,Lpq(ipq,:),NApq,Lpq(ipq,:),NApq)
                end do
              end if

            else
              do ipq=1,NApq
                Lpq_J(1:JNUM) = Lpq(ipq,1:JNUM)
                if (Do_int) then
                  kt = kOff(iSymp)+ipq-1
                  Xint(kt) = Xint(kt)+ddot_(JNUM,Lpq_J,1,Lpq_J,1)
                end if
                idisk = iOffB(iSymp)+NumCho(jSym)*(ipq-1)
#               ifdef _HDF5_QCM_
                ! Write the transformed Cholesky batch to the hdf5 dataset
                ! The ordering in HDF5 is in column-major order, corresponding to
                ! the 'Kpq' storage
                ! See above for more explanation
                if (ihdf5 == 1) then
                  call hdf5_write_cholesky(choset_id,space_id,ipq-1,nVec*(iBatch-1),JNUM,Lpq_J)
                end if
#               endif

                call ddafile(LunChVF(jSym),1,Lpq_J,JNUM,idisk)
              end do
              iOffB(iSymp) = iOffB(iSymp)+JNUM
            end if

          end if

          call CWTIME(TCR4,TWR4)
          tread(1) = tread(1)+(TCR4-TCR3)
          tread(2) = tread(2)+(TWR4-TWR3)

          call mma_deallocate(Lpq)

        end do

      end if

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

    end do  ! end batch loop

    ! free memory
    call mma_deallocate(Lpq_J)
    call Deallocate_DT(ChoT(1))
    call mma_deallocate(Lrs)

  end do   ! loop over red sets

# ifdef _HDF5_QCM_
  ! close Cholesky HDF5 stuff
  if (ihdf5 == 1) then
    call hdf5_close_cholesky(choset_id,space_id)
  else
# endif
    call daclos(LunChVF(jSym))
# ifdef _HDF5_QCM_
  end if
# endif

end do   !loop over JSYM

call CWTIME(TOTCPU2,TOTWALL2)
TOTCPU = TOTCPU2-TOTCPU1
TOTWALL = TOTWALL2-TOTWALL1

!---- Write out timing information
if (timings) then

  CFmt = '(6x,A)'
  write(u6,*)
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,CFmt) 'Cholesky-MOTRA timings            CPU       WALL '
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'

  if (Do_int) then
    write(u6,'(6x,A28,2f10.2)') 'I/O vectors + diag ERIs step',tread(1),tread(2)
  else
    write(u6,'(6x,A28,2f10.2)') 'I/O vectors',tread(1),tread(2)
  end if
  write(u6,'(6x,A28,2f10.2)') '1st half-transf.',tmotr1(1),tmotr1(2)
  write(u6,'(6x,A28,2f10.2)') '2nd half-transf.',tmotr2(1),tmotr2(2)
  write(u6,*)
  write(u6,'(6x,A28,2f10.2)') 'TOTAL',TOTCPU,TOTWALL
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,*)

end if

rc = 0

write(u6,*)
if (tv2disk == 'PQK') then
  tv2disk(1:2) = 'pq'
  write(u6,*) '     Transformed Cholesky vectors stored as L(',tv2disk(1:2),',',tv2disk(3:3),')'
else
  tv2disk(2:3) = 'pq'
  write(u6,*) '     Transformed Cholesky vectors stored as L(',tv2disk(1:1),',',tv2disk(2:3),')'
end if
write(u6,*)

return

end subroutine CHO_TR_drv
