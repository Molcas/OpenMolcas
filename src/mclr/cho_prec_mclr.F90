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
! Copyright (C) Mickael G. Delcey                                      *
!***********************************************************************

subroutine CHO_Prec_MCLR(CMO,nIsh,nAsh,LuAChoVec,LuChoInt)
!***********************************************************************
!                                                                      *
!  Author : M. G. Delcey                                               *
!                                                                      *
!  Purpose: form 2-electron integrals needed for the preconditioner    *
!           those are (ii|pq) and (ip|iq) with i inactive              *
!           and p, q active+virtual                                    *
!           as well as (tu|pq) and (tp|uq) with t, u active            *
!           and p, q all molecular orbitals                            *
!                                                                      *
!           For small active spaces, this is much less than            *
!           the full list of integrals!                                *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use Cholesky, only: InfVec, nBas, nDimRS, nSym, NumCho
use Data_structures, only: Allocate_DT, Deallocate_DT, DSBA_Type, SBA_Type
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Constants, only: One, Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: CMO(*)
integer(kind=iwp), intent(in) :: nIsh(8), nAsh(8), LuAChoVec(8), LuChoInt(2)
#include "warnings.h"
integer(kind=iwp) :: i, iAdr, iAdr2, iAdrtu, iE, ii, iLoc, ioff, ioff2, ioff3, ip1, ip2, ip3, ipiaib, ipInt, ipIntj, ipLtu, ipMO, &
                     ipMO2, ipMOi, ipMOj, iptpuq, irc, IREDC, iS, ISTSQ(8), isum, isum2, iSwap, ISYM, iSyma, isymb, it, itt, itu, &
                     iu, iuu, IVEC2, iVrs, j, jBatch, JNUM, JRED, JRED1, JRED2, jsym, JVEC, k, kMOs, ksym, ksym2, labatch, &
                     Libatch, LREAD, lvec, LWORK, lWorke, maxpq, maxRS, maxtpq, MaxVecPerBatch, memneeded, MUSED, na2, nab, nab2, &
                     naleft, nAshb(8), nAshe(8), nBatch, nileft, nip, nIshb(8), nIshe(8), nMat, nMOs, npq, nRS, ntota, ntotae, &
                     ntoti, ntotie, ntp, ntue, NumCV, NUMV, nVec, nvirt, nvirt2, nVrs
real(kind=wp) :: TCR1, TCR2, TCstart1, TCstart2, tform(2), tform2(2), tforma(2), tforma2(2), tMO(2), TOTCPU, TOTWALL, tread(2), &
                 ttran(2), TWR1, TWR2, TWstart1, TWstart2
real(kind=wp), target :: Dum(1)
logical(kind=iwp) :: add, taskleft, timings
character(len=50) :: CFmt
type(DSBA_Type) :: CMOt, Tmp(1)
type(SBA_Type) :: Lpq(1)
real(kind=wp), allocatable :: iiab(:), Integral(:), Lrs(:,:), tupq(:)
real(kind=wp), allocatable, target :: iirs(:), Lii(:), Lij(:), turs(:)
real(kind=wp), pointer :: piirs(:,:), pLii(:,:), pLij(:,:), pturs(:,:)
logical(kind=iwp), parameter :: DoRead = .false.
integer(kind=iwp), external :: Cho_LK_MaxVecPerBatch

!                                                                      *
!***********************************************************************
!                                                                      *

timings = .false.
call CWTIME(TCstart1,TWstart1)
! 1 --> CPU   2 --> Wall
tread(:) = Zero   !time read vectors
ttran(:) = Zero   !time transform vectors
tform(:) = Zero   !time form integrals
tform2(:) = Zero  !time form integrals
tforma(:) = Zero  !time form integrals
tforma2(:) = Zero !time form integrals
tMO(:) = Zero     !time for final MO transform
MaxVecPerBatch = Cho_LK_MaxVecPerBatch()
iLoc = 3

! dummy association to keep the compiler happy
piirs(0:0,0:0) => Dum(:)
pturs(0:0,0:0) => Dum(:)

ISTSQ(1) = 0
do ISYM=2,NSYM
  ISTSQ(iSYM) = ISTSQ(iSYM-1)+NBAS(ISYM-1)**2
end do

! Start with big loop over symmetries

do jsym=1,nsym
  iAdr = 0

  ! Compute some sizes

  nip = 0
  ntp = 0
  npq = 0
  maxpq = 0
  maxtpq = 0
  ntoti = 0
  ntota = 0
  do isymb=1,nsym
    iSyma = Mul(iSymb,jsym)

    npq = npq+nBas(iSymb)**2
    maxpq = max(npq,nBas(iSymb)**2)
    maxtpq = max(maxtpq,nAsh(iSyma)*nBas(iSymb)**2)

    ! For inactive half-transformed Cho vector + Lii^J
    nip = nip+nIsh(iSyma)*(nBas(isymb)+1)
    ! For active half-transformed Cho vector + Lij^J
    ntp = ntp+nAsh(iSyma)*(nBas(isymb)+nAsh(iSyma))

    ntoti = ntoti+nIsh(isymb)
    ntota = ntota+nAsh(isymb)
  end do

  NumCV = NumCho(jSym)
  call GAIGOP_SCAL(NumCV,'max')
  if (NumCV < 1) cycle

  JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
  JRED2 = InfVec(NumCho(jSym),2,jSym) ! red set of the last vec
  call GAIGOP_SCAL(JRED1,'min')
  call GAIGOP_SCAL(JRED2,'max')
  maxRS = max(0,maxval(nDimRS(JSYM,JRED1:JRED2)))

  ! Check for maxmem to see if integrals would fit in memory
  ! or if batching is required

  call mma_MaxDBLE(LWORK)
  call GAIGOP_SCAL(LWORK,'min')

  nIshb(1:nSym) = 0
  nAshb(1:nSym) = 0

  ! Loop over i and t batches
  ! each batch will increase the i/o, therefore we preferably want to
  ! have the whole iiab and tupq in-core!

  ! First, do we have enough memory at all!

  memneeded = maxRS+max(nip,ntp)                ! for 1 Jbatch
  memneeded = max(memneeded,maxpq)              ! for MO transform
  if (ntota > 0) then
    memneeded = memneeded+maxtpq                ! for 1 (ta|ub)
    ! for 1 (tu|ab) and 1 reduced set
    if (jsym == 1) memneeded = memneeded+ntota*(maxpq+maxRS)
  else if (ntoti > 0) then
    memneeded = memneeded+maxpq                 ! for 1 (ia|ib)
    ! for 1 (ii|ab) and 1 reduced set
    if (jsym == 1) memneeded = memneeded+npq+maxRS
  end if

  if (memneeded > lWork) then
    write(u6,*) 'CHO_PREC_MCLR: Insufficient memory for I/T batch'
    write(u6,*) 'LWORK= ',LWORK
    write(u6,*) 'min. mem. need= ',memneeded
    write(u6,*) 'maxRS+nip/ntp= ',maxRS+max(nip,ntp)
    write(u6,*) 'maxpq        = ',maxpq
    if (ntota > 0) then
      write(u6,*) 'maxtpq       = ',maxtpq
      if (jsym == 1) write(u6,*) 'ntota*(maxpq,maxRS) = ',ntota*maxpq,ntota*maxRS
    else
      if (jsym == 1) write(u6,*) 'npq,maxRS = ',npq,maxRS
    end if
    call Quit(_RC_MEMORY_ERROR_)
  end if
  lWork = lWork-max(maxpq,maxRS+max(nip,ntp))

  ! How many nab + nrs can be stored in memory while stil having place for
  ! at least 1 nRS, 1 maxpq and 1 nIP
  ! (actually nIP will be smaller but that does not matter very much)
  !
  ! nIshe is the number of inactive orbitals dealt in the batch
  ! nIshb is inactive dealt with after the previous batches

  taskleft = .true.
  do while (taskleft)
    lWorke = lWork
    nIshe(1:nSym) = 0
    nAshe(1:nSym) = 0
    taskleft = .false.
    Libatch = 0
    do i=1,nsym
      k = Mul(i,jsym)

      nab = 0
      nRS = 0
      if (jsym == 1) then
        nab = npq
        nRS = maxRS
      end if
      nab2 = nBas(k)**2
      nileft = nIsh(i)-nIshb(i)
      if ((nab+nab2)*nileft > 0) then
        nIshe(i) = min(lWorke/(nrs+nab+nab2),nileft)
        lWorke = lWorke-(nrs+nab+nab2)*nIshe(i)
        libatch = libatch+(nrs+nab+nab2)*nIshe(i)
        if (nIshe(i) /= nileft) then
          taskleft = .true.
          exit
        end if
      end if

    end do

    ! Update nip and compute sum(nIshe)

    nip = 0
    ntotie = 0
    do i=1,nsym
      k = Mul(i,jsym)
      nip = nip+nIshe(i)*nBas(k)
      ntotie = ntotie+nIshe(i) ! For Lii^J
    end do
    nip = nip+ntotie ! for Lii^J

    call mma_allocate(iiab,libatch,Label='iiab')
    nab = 0
    if ((jsym == 1) .and. (ntotie > 0)) then
      nab = npq
      call mma_allocate(iirs,ntotie*maxRS,Label='iirs')
    end if
    ipiaib = 1+nab*ntotie
    iiab(:) = Zero

    if (taskleft) then
      ntotae = 0
      ntue = 0
      ntp = 0      ! do not allocate those
      labatch = 0
      iptpuq = -1
      !write(u6,*) 'Batching loop i'
    else

      ! Batching T loop

      labatch = 0
      do i=1,nsym
        k = Mul(i,jsym)
        nab = 0
        nRS = 0
        if (jsym == 1) then
          nab = ntota*npq
          nRS = ntota*maxRS
        end if
        nab2 = nAsh(i)*nBas(k)**2
        naleft = nAsh(i)-nAshb(i)
        if ((nab+nab2)*naleft > 0) then
          nAshe(i) = min(lWorke/(nrs+nab+nab2),naleft)
          lWorke = lWorke-(nrs+nab+nab2)*nAshe(i)
          labatch = labatch+(nrs+nab+nab2)*nAshe(i)
          if (nAshe(i) /= naleft) then
            taskleft = .true.
            exit
          end if
        end if
      end do

      ntotae = 0
      ntue = 0
      ntp = 0
      do i=1,nsym
        k = Mul(i,jsym)
        ntotae = ntotae+nAshe(i)
        ntue = ntue+nAshe(i)*nAsh(i)
        ntp = ntp+nAsh(i)*(nBas(k)+nAshe(i))
      end do

      if (ntotae > 0) call mma_allocate(tupq,labatch,Label='tupq')
      nab = 0
      if (jsym == 1) then
        nab = npq
        if (ntue > 0) call mma_allocate(turs,ntue*maxRS,Label='turs')
      end if
      iptpuq = 1+nab*ntue
      tupq(:) = Zero
      if (taskleft) write(u6,*) 'Batching loop a'
    end if

    ! Transpose CMO

    call Allocate_DT(CMOt,nIShe,nBas,nSym)

    ioff = 0
    do iSym=1,nsym
      do j=1,nIshe(iSym)
        ioff3 = ioff+nBas(iSym)*(nIshb(iSym)+j-1)
        CMOt%SB(iSym)%A2(j,:) = CMO(ioff3+1:ioff3+nBas(iSym))
      end do
      ioff = ioff+nBas(iSym)**2
    end do

    ! Loop over reduced sets

    do JRED=JRED1,JRED2
      call Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)

      if (nVrs == 0) cycle  ! no vectors in that (jred,isym)

      if (nVrs < 0) then
        write(u6,*) 'CHO_PREC_MCLR: Cho_X_nVecRS returned nVrs<0. STOP!'
        call Abend()
      end if

      call Cho_X_SetRed(irc,iLoc,JRED)
      ! set index arrays at iLoc
      if (irc /= 0) then
        write(u6,*) 'CHO_PREC_MCLR: cho_X_setred non-zero return code. rc= ',irc
        call Abend()
      end if

      IREDC = JRED

      nRS = nDimRS(JSYM,JRED)

      if (jSym == 1) then
        piirs(1:nRS,1:ntotie) => iirs(1:nRS*ntotie)
        piirs(:,:) = Zero
        pturs(1:nRS,1:ntue) => turs(1:nRS*ntue)
        pturs(:,:) = Zero
      end if

      call mma_MaxDBLE(LWORKe)
      nVec = min(LWORKE/(nRS+max(nip,ntp)),min(nVrs,MaxVecPerBatch))
      if (nVec < 1) then
        write(u6,*) 'CHO_PREC_MCLR: Insufficient memory for J batch'
        write(u6,*) 'That should not happen here'
        write(u6,*) 'Contact the developers'
        call Quit(_RC_MEMORY_ERROR_)
        nBatch = -9999  ! dummy assignment
      end if
      LREAD = nRS*nVec

      call mma_allocate(Lrs,nRS,nVec,Label='Lrs')

      nBatch = (nVrs-1)/nVec+1

      ! Loop over batches of Cholesky vectors

      do jBatch=1,nBatch

        if (jBatch == nBatch) then
          JNUM = nVrs-nVec*(nBatch-1)
        else
          JNUM = nVec
        end if
        JVEC = nVec*(jBatch-1)+iVrs
        IVEC2 = JVEC-1+JNUM

        iSwap = 1 ! Lqi,J are returned
        call Allocate_DT(Lpq(1),nIshe,nBas,JNUM,JSYM,nSym,iSwap)
        call mma_allocate(Lii,ntotie*nVec,Label='Lii')
        !***************************************************************
        !                                                              *
        !          Let's start the real work                           *
        !                                                              *
        !***************************************************************

        ! Read Cholesky vectors

        call CWTIME(TCR1,TWR1)

        call CHO_VECRD(Lrs,LREAD,JVEC,IVEC2,JSYM,NUMV,IREDC,MUSED)

        if ((NUMV <= 0) .or. (NUMV /= JNUM)) return

        call CWTIME(TCR2,TWR2)
        tread(1) = tread(1)+(TCR2-TCR1)
        tread(2) = tread(2)+(TWR2-TWR1)

        !***************************************************************
        ! MO Half-transformation
        ! Liq^J= sum_p Lpq^J Xip

        kMOs = 1  !
        nMOs = 1  ! Active MOs (1st set)

        call CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,JSYM,iSwap,IREDC,nMOs,kMOs,[CMOt],Lpq(1),DoRead)

        if (irc /= 0) return

        call CWTIME(TCR1,TWR1)
        ttran(1) = ttran(1)+(TCR1-TCR2)
        ttran(2) = ttran(2)+(TWR1-TWR2)

        !***************************************************************
        ! Integral formation
        ! (i p | i q) = sum_J Lip^J  Liq^J

        ip1 = ipiaib
        do isym=1,nsym
          ksym = Mul(iSym,jsym)
          do ii=1,nIshe(isym)

            call DGEMM_('N','T',nBas(kSym),nBas(kSym),JNUM,One,Lpq(1)%SB(kSym)%A3(:,ii,1),nBas(kSym)*nIshe(iSym), &
                        Lpq(1)%SB(kSym)%A3(:,ii,1),nBas(kSym)*nIshe(iSym),One,iiab(ip1:),nBas(kSym))
            ip1 = ip1+nBas(kSym)**2
          end do
        end do
        call CWTIME(TCR2,TWR2)
        tform(1) = tform(1)+(TCR2-TCR1)
        tform(2) = tform(2)+(TWR2-TWR1)

        !***************************************************************
        ! Second MO transformation
        ! Lii^J = sum_q Liq^J Xiq

        if (jSym == 1) then
          ipMO = 1

          iE = 0
          do isym=1,nsym
            iS = 1+iE
            iE = iE+JNUM*nIshe(iSym)

            pLii(1:JNUM,1:nIshe(iSym)) => Lii(iS:iE)

            do ii=1,nIshe(iSym)
              ipMO = ipMO+ISTSQ(iSym)
              ipMOi = ipMO+(nIshb(isym)+ii-1)*nBas(iSym)
              call dGeMV_('T',nBas(iSym),JNUM,One,Lpq(1)%SB(iSym)%A3(:,ii,1),nBas(iSym)*nIshe(iSym),CMO(ipMOi),1,Zero,pLii(:,ii),1)
            end do
            nullify(pLii)
          end do

          call CWTIME(TCR1,TWR1)
          ttran(1) = ttran(1)+(TCR1-TCR2)
          ttran(2) = ttran(2)+(TWR1-TWR2)

          !*************************************************************
          ! Integral formation
          ! (i i | p q) = sum_J Lii^J  Lpq^J

          pLii(1:JNUM,1:ntotie) => Lii(1:JNUM*ntotie)
          call DGEMM_('N','N',nRS,ntotie,JNUM,One,Lrs,nRS,pLii,JNUM,One,piirs,nRS)
          nullify(pLii)

          call CWTIME(TCR2,TWR2)
          tform2(1) = tform2(1)+(TCR2-TCR1)
          tform2(2) = tform2(2)+(TWR2-TWR1)

        end if  ! jSym

        call mma_deallocate(Lii)
        call Deallocate_DT(Lpq(1))

        !***************************************************************
        ! Read half-transformed active vectors

        iSwap = 0 ! Lvb,J
        call Allocate_DT(Lpq(1),nAsh,nBas,JNUM,JSYM,nSym,iSwap)

        call mma_allocate(Lij,ntue*nVec,Label='Lij')

        do i=1,nSym
          k = Mul(i,JSYM)
          lvec = nAsh(k)*nBas(i)*JNUM
          iAdr2 = (JVEC-1)*nAsh(k)*nBas(i)
          call DDAFILE(LuAChoVec(Jsym),2,Lpq(1)%SB(i)%A3,lvec,iAdr2)
        end do

        !***************************************************************
        ! Form (tp|uq) integrals

        do j=1,JNUM

          ioff = 0
          do i=1,nsym
            k = Mul(i,JSYM)

            do it=0,nAshe(k)-1
              itt = nAshb(k)+it+1

              do iu=0,nAshb(k)+it
                iuu = iu+1

                itu = nTri_Elem(it)+it*nAshb(k)+iu

                ipInt = iptpuq+ioff+itu*nBas(i)**2

                call DGER_(nBas(i),nBas(i),One,Lpq(1)%SB(i)%A2(itt:,j),nAsh(i),Lpq(1)%SB(i)%A2(iuu:,j),nAsh(i),tupq(ipInt),nBas(i))

              end do
            end do

            ioff = ioff+(nTri_Elem(nAshe(k))+nAshe(k)*nAshb(k))*nBas(i)**2
          end do
        end do
        call CWTIME(TCR1,TWR1)
        tforma(1) = tforma(1)+(TCR1-TCR2)
        tforma(2) = tforma(2)+(TWR1-TWR2)

        !***************************************************************
        ! Second MO transformation

        if (jsym == 1) then
          ipLtu = 1

          ioff = 0
          iE = 0
          do i=1,nsym
            iS = iE+1
            do j=1,JNUM

              ipMO = 1+ioff+nBas(i)*(nIsh(i)+nAshb(i))
              do k=0,nAshe(i)-1
                call dGeMV_('N',nAshb(i)+k+1,nBas(i),One,Lpq(1)%SB(i)%A3(:,1,j),nAsh(i),CMO(ipMO+k*nBas(i)),1,Zero,Lij(ipLtu),1)
                ipLtu = ipLtu+(nAshb(i)+k+1)
              end do
            end do
            ioff = ioff+nBas(i)**2
          end do

          call CWTIME(TCR2,TWR2)

          !*************************************************************
          ! Formation of the (tu|rs) integral

          ipInt = 1
          iE = 0
          do i=1,nsym
            na2 = nAshe(i)*nAshb(i)+nTri_Elem(nAshe(i))
            if (na2 == 0) cycle
            iS = iE+1
            iE = iE+na2*JNUM

            pLij(1:na2,1:JNUM) => Lij(iS:iE)

            call DGEMM_('N','T',nRS,na2,JNUM,One,Lrs,nRS,pLij,na2,One,turs(ipInt),nRS)
            ipInt = ipInt+nRS*na2
            nullify(pLij)
          end do
          call CWTIME(TCR1,TWR1)
          tforma2(1) = tforma2(1)+(TCR1-TCR2)
          tforma2(2) = tforma2(2)+(TWR1-TWR2)
        end if
        !***************************************************************
        !                                                              *
        !          Cholesky loop is over!                              *
        !                                                              *
        !***************************************************************

        call mma_deallocate(Lij)
        call Deallocate_DT(Lpq(1))
      end do ! jbatch

      ! Transform to full storage, use Lrs as temp storage

      if (jsym == 1) then
        add = .true.
        nMat = 1
        do i=1,ntotie
          call Allocate_DT(Tmp(1),nBas,nBas,nSym,aCase='TRI',Ref=iiab(1+nab*(i-1):))
          call swap_rs2full(irc,iLoc,nRS,nMat,JSYM,Tmp,piirs(:,i),add)
          call Deallocate_DT(Tmp(1))
        end do
        do i=1,ntue
          call Allocate_DT(Tmp(1),nBas,nBas,nSym,aCase='TRI',Ref=tupq(1+npq*(i-1):))
          call swap_rs2full(irc,iLoc,nRS,nMat,JSYM,Tmp,pturs(:,i),add)
          call Deallocate_DT(Tmp(1))
        end do
      end if

      call mma_deallocate(Lrs)
    end do ! reduced set JRED

    if (jsym == 1) then
      if (ntotie > 0) call mma_deallocate(iirs)
      if (ntue > 0) call mma_deallocate(turs)
      nullify(piirs,pturs)
    end if

    !MGD  Gather integrals from parallel runs

    call CWTIME(TCR1,TWR1)
    call mma_allocate(Integral,maxpq,Label='Integral')
    ip1 = ipiaib
    ip3 = iptpuq
    isum = 0
    isum2 = 0
    do isym=1,nSym
      ksym = Mul(iSym,jsym)
      nvirt = nBas(kSym)-nIsh(kSym)
      ipMO = 1+ISTSQ(kSym)+nBas(kSym)*nIsh(kSym)

      do ii=1,nIshe(isym)

        ! MO transform (i i | p q)

        if (jsym == 1) then
          isum = isum+1
          ip2 = 1+nab*(isum-1)
          ioff2 = 0
          do ksym2=1,nsym
            ipMO2 = 1+ioff2+nBas(kSym2)*nIsh(kSym2)
            nvirt2 = nBas(kSym2)-nIsh(kSym2)
            do j=1,nvirt2
              ipMOj = ipMO2+(j-1)*nBas(kSym2)
              ipIntj = 1+(j-1)*nBas(kSym2)
              call DSPMV_('U',nBas(kSym2),One,iiab(ip2),CMO(ipMOj),1,Zero,Integral(ipIntj),1)
            end do
            call DGEMM_('T','N',nvirt2,nvirt2,nBas(kSym2),One,Integral,nBas(kSym2),CMO(ipMO2),nBas(kSym2),Zero,iiab(ip2),nvirt2)

            call GADSum(iiab(ip2),nvirt2**2)
            call DDAFILE(LuChoInt(1),1,iiab(ip2),nvirt2**2,iAdr)

            ip2 = ip2+nBas(kSym2)**2
            ioff2 = ioff2+nBas(kSym2)**2
          end do
        else
          iAdr = iAdr+sum((nBas(1:nSym)-nIsh(1:nSym))**2)
        end if

        ! MO transform (i p | i q)

        call DGEMM_('N','N',nBas(kSym),nvirt,nBas(kSym),One,iiab(ip1:),nBas(kSym),CMO(ipMO),nBas(kSym),Zero,Integral,nBas(kSym))
        call DGEMM_('T','N',nvirt,nvirt,nBas(kSym),One,Integral,nBas(kSym),CMO(ipMO),nBas(kSym),Zero,iiab(ip1:),nvirt)

        iAdr = iAdr+sum((nBas(1:kSym-1)-nIsh(1:kSym-1))**2)
        call GADSum(iiab(ip1:),nvirt**2)
        call DDAFILE(LuChoInt(1),1,iiab(ip1:),nvirt**2,iAdr)
        iAdr = iAdr+sum((nBas(kSym+1:nSym)-nIsh(kSym+1:nSym))**2)

        ip1 = ip1+nBas(kSym)**2
      end do

      ! MO transform (t u | p q)

      ! compute address
      iAdrtu = 0
      do i=1,ksym-1
        na2 = nAshe(i)*nAshb(i)+nTri_Elem(nAshe(i))
        do j=1,na2
          iAdrtu = iAdrtu+npq+sum(nBas(1:nSym)**2)
        end do
      end do

      ipMO = 1+ISTSQ(iSym)
      na2 = nAshe(ksym)*nAshb(ksym)+nTri_Elem(nAshe(ksym))
      do itu=1,na2
        if (jsym == 1) then
          isum2 = isum2+1
          ip2 = 1+npq*(isum2-1)
          ioff2 = 0
          do ksym2=1,nsym
            ipMO2 = 1+ioff2
            do j=1,nBas(ksym2)
              ipMOj = ipMO2+(j-1)*nBas(kSym2)
              ipIntj = 1+(j-1)*nBas(kSym2)
              call DSPMV_('U',nBas(kSym2),One,tupq(ip2),CMO(ipMOj),1,Zero,Integral(ipIntj),1)
            end do
            if (nBas(ksym2) > 0) then
              call DGEMM_('T','N',nBas(ksym2),nBas(kSym2),nBas(kSym2),One,Integral,nBas(kSym2),CMO(ipMO2),nBas(kSym2),Zero, &
                          tupq(ip2),nBas(kSym2))
              call GADSum(tupq(ip2),nBas(kSym2)**2)
              call DDAFILE(LuChoInt(2),1,tupq(ip2),nBas(kSym2)**2,iAdrtu)

              ip2 = ip2+nBas(kSym2)**2
              ioff2 = ioff2+nBas(kSym2)**2
            end if
          end do
        else
          iAdrtu = iAdrtu+npq
        end if

        ! MO transform (t p | u q)

        call DGEMM_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),One,tupq(ip3),nBas(iSym),CMO(ipMO),nBas(iSym),Zero,Integral,nBas(iSym))
        call DGEMM_('T','N',nBas(iSym),nBas(iSym),nBas(iSym),One,Integral,nBas(iSym),CMO(ipMO),nBas(iSym),Zero,tupq(ip3),nBas(iSym))
        iAdrtu = iAdrtu+sum(nBas(1:iSym-1)**2)

        call GADSum(tupq(ip3),nBas(iSym)**2)
        call DDAFILE(LuChoInt(2),1,tupq(ip3),nBas(iSym)**2,iAdrtu)

        iAdrtu = iAdrtu+sum(nBas(iSym+1:nSym)**2)
        ip3 = ip3+nBas(iSym)**2

      end do

    end do
    call CWTIME(TCR2,TWR2)
    tMO(1) = tMO(1)+(TCR2-TCR1)
    tMO(2) = tMO(2)+(TWR2-TWR1)

    call mma_deallocate(Integral)

    nIshb(1:nSym) = nIshb(1:nSym)+nIshe(1:nSym)  ! now those are done!
    nAshb(1:nSym) = nAshb(1:nSym)+nAshe(1:nSym)  ! now those are done!
    call Deallocate_DT(CMOt)
    call mma_deallocate(iiab)
    if (ntotae > 0) call mma_deallocate(tupq)
  end do  ! loop over i/t batches

end do ! jsym

call CWTIME(TCstart2,TWstart2)
TOTCPU = TCstart2-TCstart1
TOTWALL = TWstart2-TWstart1
if (timings) then

  CFmt = '(2x,A)'
  write(u6,*)
  write(u6,CFmt) 'Cholesky MCLR timing from CHO_PREC_MCLR'
  write(u6,CFmt) '----------------------------------------'
  write(u6,*)
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,CFmt) 'Integral construction           CPU       WALL   '
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'

  write(u6,'(2x,A26,2f10.2)') 'READ VECTORS                              ',tread(1),tread(2)
  write(u6,'(2x,A26,2f10.2)') 'TRANSFORMATION                            ',ttran(1),ttran(2)
  write(u6,'(2x,A26,2f10.2)') '(IA|IB) FORMATION                         ',tform(1),tform(2)
  write(u6,'(2x,A26,2f10.2)') '(II|AB) FORMATION                         ',tform2(1),tform2(2)
  write(u6,'(2x,A26,2f10.2)') '(TP|UQ) FORMATION                         ',tforma(1),tforma(2)
  write(u6,'(2x,A26,2f10.2)') '(TU|PQ) FORMATION                         ',tforma2(1),tforma2(2)
  write(u6,'(2x,A26,2f10.2)') 'MO TRANSFORMATION                         ',tMO(1),tMO(2)
  write(u6,*)
  write(u6,'(2x,A26,2f10.2)') 'TOTAL                                     ',TOTCPU,TOTWALL
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,*)

end if

end subroutine CHO_Prec_MCLR
