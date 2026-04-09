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
! Copyright (C) 1996, Anders Bernhardsson                              *
!               1996, Bjorn O. Roos                                    *
!***********************************************************************

subroutine TRAMO_MCLR(LBUF,X1,n1,X2,n2,X3,n3,X4,n4,Buffer,MEMX,NBP,NBQ,NBR,NBS,iSP,iSQ,iSR,iSS,nAP,nAQ,nAR,nAS,CMP,CMQ,CMR,CMS, &
                      iAD13,iAD14,iAD23,iAD24,iAD34,IDAHLF2,IRLHLF2,IDAHLF3,IRLHLF3,LIOTAB)
!***********************************************************************
!                                                                      *
!     Purpose: two-electron transformation routine.                    *
!              Transformation is made in core if all                   *
!              half-transformed integrals fit.                         *
!              Otherwise sorted integrals are written onto             *
!              unit LUHALF.                                            *
!                                                                      *
!***********************************************************************
!                                                                      *
!  Stolen from MOTRA and destroyed by Anders to make it                *
!  possible to generate all sort of combination of transformed         *
!  and untransformed indexes                                           *
!                                                                      *
!***********************************************************************

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Index_Functions, only: nTri_Elem
use MCLR_Data, only: FnHlf2, FnHlf3, LuHlf2, LuHlf3, LuTri1, LuTri2, LuTri3, LuTri4, LuTri5, NoFile
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6, ItoB

implicit none
integer(kind=iwp), intent(in) :: LBUF, n1, n2, n3, n4, MEMX, NBP, NBQ, NBR, NBS, iSP, iSQ, iSR, iSS, nAP, nAQ, nAR, nAS, LIOTAB
real(kind=wp), intent(out) :: X1(n1), X2(n2), X3(n3), X4(n4), Buffer(MemX)
real(kind=wp), intent(in) :: CMP(nBP,nAP), CMQ(nBQ,nAQ), CMR(nBR,nAR), CMS(nBS,nAS)
integer(kind=iwp), intent(inout) :: iAD13, iAD14, iAD23, iAD24, iAD34
integer(kind=iwp), intent(out) :: IDAHLF2(LIOTAB), IRLHLF2(LIOTAB), IDAHLF3(LIOTAB), IRLHLF3(LIOTAB)
integer(kind=iwp) :: I1, I2, IAD2, IAD3, IBUF, IBUF2, IBUF3, IMAX, INCORE, iOne, IOPT, IP2, IP3, IP4, IP5, IPB, IPQ, IPQST, &
                     IPQUT2, IPQUT3, IPQUT4, IPX, IPY, IPZ, IRC, IRSST, IST, IST2, IST3, IVX, LPKREC, LPQ, MemXX, MEMXXX, NARS, &
                     NBPQ, NBRS, NBUF, NBYTES, NOUT, NPQ, NQM, NR2, NR3

call TRAMO_MCLR_INTERNAL(Buffer)

! This is to allow type punning without an explicit interface
contains

subroutine TRAMO_MCLR_INTERNAL(Buffer)

  real(kind=wp), target, intent(out) :: Buffer(MemX)
  integer(kind=iwp), pointer :: iBuffer(:)
  integer(kind=iwp) :: I, NP, NQ, NR, NS

  ione = 1
  if (nofile) ione = 0

  ! Set some constants

  MemXX = MemX-lBuf-1010
  if (ISS == isr) then
    NBPQ = nTri_Elem(NBP)
    NBRS = nTri_Elem(NBR)
    NARS = nTri_Elem(NAR)
  else
    NBPQ = NBP*NBQ
    NBRS = NBR*NBS
    NARS = NAR*NAS
  end if
  iBUF2 = 0
  iBUF3 = 0
  iAD2 = 0      ! Disk address for (kL|ij)
  iAD3 = 0      ! Disk address for (Kl|ij)
  iPQUT2 = 0    ! Index(ij)  in buffer.
  iPQUT3 = 0    ! Just one should be enough but...
  iPQUT4 = 0    ! ...easy is better than intelligent
  if (iSP == iSQ) then
    iMax = nTri_Elem(nBP)
  else
    iMax = nBP*nBQ  ! Size of buffer
  end if
  Incore = 0        ! First we try to  transform in core
  if (iss /= isr) then
    nBuf = nBR*nAS+nAR*nBS
  else
    nbuf = nbr*nas
  end if

  ! Buffer needed in the final step in the calculation of the exchange type
  ! integrals. The factor ten is arbitrary, however a 10 writes to disk per
  ! batch is acceptable before paging out integrals to disk to increase
  ! memory, however we dont want to write in to small batches 512 dWords (4k) is
  ! the lower limit
  !
  ! /*
  !    __________________________________
  !    |          |           |         |
  !    |          |           |         |
  !    ----------------------------------
  !    |          |           |         |
  !    |         -|-----------|---------|
  !    ---------------------------------|\
  !    |        | |           |         | \
  !    | ORDINT | |  (IjkL)   | (Ij|Kl) |  \
  !    |----------|           |         |   \ N
  !    ----------------------------------   /
  !    |          |           |         |  /
  !    | (ij|KL)  |  (iJ|kL)  | (iJ|Kl) | /
  !    |          |           |         |/
  !    ----------------------------------\
  !    |          |           |         | \
  !    | (ij|Kl)  |  (ij|Kl)  | (ij|Kl) |  \
  !    |          |           |         |   \
  !    ----------------------------------    nBuf*iMax
  !    |          |           |         |   /
  !    | (ij|kL)  |  (ij|kL)  | (ij|kL) |  /
  !    |          |           |         | /
  !    ----------------------------------/
  ! */
  i1 = nAS*(nAP*nBQ*nBR+nAQ*nBP*nBR)
  if (ISP == ISQ) i1 = nAS*nAP*nBQ*nBR
  i2 = nAR*(nAP*nBQ*nBS+nAQ*nBP*nBS)
  if (ISP == ISQ) i2 = nAR*nAP*nBQ*nBS
  ! Size of output buffer if all data should fit in memory
  NOUT = max(i1,i2,NARS*nBPQ)
  !NOUT = max(nAS*(nAP*nBQ*nBR+nAQ*nBP*nBR),nAR*(nAP*nBQ*nBS+nAQ*nBP*nBS),NAR*nAS*nBQ*nBP)

  if (nBuf*iMax+NOUT <= MemXX) then
    Incore = 0 ! We can do the calculation totally in core
  else
    ! But we dont have enough memory
    Incore = 1
    MemXXX = MemXX-nOut
    ! So we give each buffer same amount of memory
    iMax = MemXXX/nBuf
    ! But we need space to unpack the integrals
    iMax = (MemXXX-iMax)/nBuf
    iMax = (MemXXX-iMax)/nBuf
    iMax = (MemXXX-iMax)/nBuf
    if (iMax < 1) then
      write(u6,*) 'TraMO_MCLR: iMax < 1'
      write(u6,*) 'iMax=',iMax
      call Abend()
    end if
    call DAName(LUHLF2,FNHLF2)
    call DAName(LUHLF3,FNHLF3)
  end if

  !*********************************************************************
  !
  ! OK to start with with we need buffers but in the final step
  ! the only memory we need is a place to load (pq) integrals plus
  ! unpacking area plus a area of the size n to put the final result
  ! if we do not have enough memory, we will have a very long nose
  !
  !*********************************************************************

  ip2 = 1                   ! (ij|kL)
  ip3 = ip2+iMax*nBR*nAS+1  ! (ij|Kl)
  Buffer(ip3-1) = -99999.0_wp
  if (iss /= isr) then
    ip4 = ip3+iMax*nBS*nAR+1  ! (ij|KL)
  else
    ip4 = ip3+1
  end if
  Buffer(ip4-1) = -99999.0_wp
  ipB = ip4+nBPQ*nARS+1   ! Read from ORDINT
  Buffer(ipB-1) = -99999.0_wp

  !*********************************************************************
  !
  ! Start loop over sorted AO-integrals: npq pq-pairs in each buffer
  !
  !*********************************************************************

  IPQ = 0
  LPQ = 0
  NPQ = 0
  IRSST = 1-NBRS    ! startpoint for present batch of integrals
  iopt = 1

  if (.not. NOFILE) then
    do NP=1,NBP
      NQM = NBQ
      if (ISP == ISQ) NQM = NP
      do NQ=1,NQM
        IPQ = IPQ+1
        IPQUT2 = IPQUT2+1
        IPQUT3 = IPQUT3+1
        IPQUT4 = IPQUT4+1

        !***************************************************************
        !                                                              *
        !         Read in a block of integrals  NPQ pq values          *
        !                                                              *
        !***************************************************************

        if ((LPQ == NPQ) .and. (.not. NOFILE)) then
          call RDORD(IRC,IOPT,ISP,ISQ,ISR,ISS,Buffer(ipB),LBUF,NPQ)
          call GADSum(Buffer(ipB),LBuf)
          if (irc /= 0) then
            write(u6,*) 'TraMO_MCLR: error reading ORDINT!'
            call Abend()
          end if
          IOPT = 2
          LPQ = 0
          IRSST = 1-NBRS
        end if
        LPQ = LPQ+1
        IRSST = IRSST+NBRS

        !***************************************************************
        !                                                              *
        !       Start transformation of this pq pair -> (ij|kL) (ijKl) *
        !                                                              *
        !***************************************************************

        if ((nAS*nAR /= 0) .or. (nAS*(nAP+nAQ) /= 0)) then
          if (iSR == iSS) then
            call SQUARE(Buffer(ipB+IRSST-1),X1,1,NBS,NBS)
            call DGEMM_('T','N',NBR,NAS,NBS,One,X1,NBS,CMS,NBS,Zero,X2,NBR)
          else
            call DGEMM_('T','N',NBR,NAS,NBS,One,Buffer(ipB+IRSST-1),NBS,CMS,NBS,Zero,X2,NBR)
          end if
        end if

        !***************************************************************
        !                                                              *
        !     Sort and save integrals transformed in one index         *
        !     that should be used for exchange type integrals          *
        !                                                              *
        !***************************************************************

        if ((nAP+nAQ)*nAS /= 0) then
          if (IPQUT2 > iMax) then
            iPQUT2 = 1
            iST2 = 1
            do i=1,nAS*nBR
              iBuf2 = iBuf2+1
              if (iBuf2 > LIOTAB) then
                write(u6,*) 'TraMO_MCLR: iBuf2 > LIOTAB'
                write(u6,*) 'iBuf2,LIOTAB=',iBuf2,LIOTAB
                call Abend()
              end if
              call PKR8(0,iMax,NBYTES,Buffer(ip2+ist2-1),Buffer(ip2+ist2-1))
              LPKREC = (NBYTES+ItoB-1)/ItoB
              ! Save the length of this record
              IRLHLF2(iBUF2) = LPKREC
              ! Save the address of this record
              IDAHLF2(iBUF2) = iAD2
              call c_f_pointer(c_loc(Buffer(ip2+IST2-1)),iBuffer,[LPKREC])
              call iDAFILE(LUHLF2,1,iBuffer,LPKREC,IAD2)
              nullify(iBuffer)
              iST2 = iST2+iMax
            end do
          end if

          Buffer(ip2+iPQUT2-1:ip2+iPQUT2-1+nAs*nBR*iMax-1:iMax) = X2(1:nAs*nBR)

        end if
        !***************************************************************

        if (ISR /= ISS) then
          call DGEMM_('N','N',NBS,NAR,NBR,One,Buffer(ipB+IRSST-1),NBS,CMR,NBR,Zero,X3,NBS)

          !*************************************************************
          if ((nAP+nAQ)*nAR /= 0) then
            if (IPQUT3 > iMax) then
              IPQUT3 = 1
              IST3 = 1
              do I=1,nAR*nBS
                IBUF3 = IBUF3+1
                if (IBUF3 > LIOTAB) then
                  write(u6,*) 'TraMO_MCLR: iBuf3 > LIOTAB'
                  write(u6,*) 'iBuf3,LIOTAB=',iBuf3,LIOTAB
                  call Abend()
                end if
                call PKR8(0,IMAX,NBYTES,Buffer(ip3+IST3-1),Buffer(ip3+IST3-1))
                LPKREC = (NBYTES+ItoB-1)/ItoB
                IRLHLF3(IBUF3) = LPKREC
                IDAHLF3(IBUF3) = IAD3
                call c_f_pointer(c_loc(Buffer(ip3+IST3-1)),iBuffer,[LPKREC])
                call iDAFILE(LUHLF3,1,iBuffer,LPKREC,IAD3)
                nullify(iBuffer)
                IST3 = IST3+iMax
              end do
            end if
          end if

          Buffer(ip3+iPQUT3-1:ip3+iPQUT3-1+nAR*nBS*iMax:iMax) = X3(1:nAR*nBS)

        end if

        !***************************************************************
        !
        !  Transform to Coulomb type integrals
        !
        !***************************************************************

        if (NAR*NAS /= 0) then
          if (ISR == ISS) then
            call DGEMM_Tri('T','N',NAR,NAR,NBR,One,X2,NBR,CMR,NBR,Zero,X4,NAR)
          else
            call DGEMM_('T','N',NAS,NAR,NBR,One,X2,NBR,CMR,NBR,Zero,X4,NAS)
          end if
        end if

        !***************************************************************

        Buffer(ip4+iPQUT4-1:ip4+iPQUT4-1+nARS*nBPQ-1:nBPQ) = X4(1:nARS)

        !***************************************************************

      end do
    end do

    !*******************************************************************
    !
    !     Empty last buffers
    !
    !*******************************************************************

    if (INCORE == 1) then
      IST = 1
      do I=1,nBR*nAS
        IBUF2 = IBUF2+1
        if (IBUF2 > LIOTAB) then
          write(u6,*) 'TraMO_MCLR: iBuf2 > LIOTAB'
          write(u6,*) 'iBuf2,LIOTAB=',iBuf2,LIOTAB
          call Abend()
        end if
        call PKR8(0,iMax,NBYTES,Buffer(ip2+IST-1),Buffer(ip2+IST-1))
        LPKREC = (NBYTES+ItoB-1)/ItoB
        IDAHLF2(IBUF2) = IAD2
        IRLHLF2(IBUF2) = LPKREC
        call c_f_pointer(c_loc(Buffer(ip2+IST-1)),iBuffer,[LPKREC])
        call iDAFILE(LUHLF2,1,iBuffer,LPKREC,IAD2)
        nullify(iBuffer)
        IST = IST+IMAX
      end do
      if (iSS /= iSR) then
        IST = 1
        do I=1,nBS*nAR
          IBUF3 = IBUF3+1
          if (IBUF3 > LIOTAB) then
            write(u6,*) 'TraMO_MCLR: iBuf3 > LIOTAB'
            write(u6,*) 'iBuf3,LIOTAB=',iBuf3,LIOTAB
            call Abend()
          end if
          call PKR8(0,iMax,NBYTES,Buffer(ip3+IST-1),Buffer(ip3+IST-1))
          LPKREC = (NBYTES+ItoB-1)/ItoB
          IDAHLF3(IBUF3) = IAD3
          IRLHLF3(IBUF3) = LPKREC
          call c_f_pointer(c_loc(Buffer(ip3+IST-1)),iBuffer,[LPKREC])
          call iDAFILE(LUHLF3,1,iBuffer,LPKREC,IAD3)
          nullify(iBuffer)
          IST = IST+IMAX
        end do
      end if
    end if
  end if ! nofile
  if (NOFILE) ipqut4 = NBPQ
  if (nars /= 0) call dDafile(LUTRI1,iOne,Buffer(ip4),nARS*NBPQ,iAD34)
  if (Buffer(ip3-1) /= -99999.0_wp) then
    write(u6,*) 'TraMO_MCLR: Buffer(ip3-1) /= -99999.0'
    write(u6,*) 'Buffer(ip3-1)=',Buffer(ip3-1)
    call Abend()
  end if
  if (Buffer(ip4-1) /= -99999.0_wp) then
    write(u6,*) 'TraMO_MCLR: Buffer(ip4-1) /= -99999.0'
    write(u6,*) 'Buffer(ip4-1)=',Buffer(ip4-1)
    call Abend()
  end if
  if (Buffer(ipB-1) /= -99999.0_wp) then
    write(u6,*) 'TraMO_MCLR: Buffer(ipB-1) /= -99999.0'
    write(u6,*) 'Buffer(ipB-1)=',Buffer(ipB-1)
    call Abend()
  end if

  if (nAP+nAQ /= 0) then
    !*******************************************************************
    !
    !     second transformation
    !     =====================
    !
    !    The integrals are transformed and written to disk
    !    this is done so that the integrals are written to
    !    disk in batches of nB*nB*nA integrals. This is the
    !    Least number of integrals that have to match in buffer
    !    because we want to write the integrals sorted onto disk
    !
    !*******************************************************************

    ipX = ip4
    ipY = ipX+nBR*NAS*nBP*nAQ+1   ! (Ij|kL)
    ipz = ipy+1
    if (isp /= isq) ipZ = ipY+nBR*NAS*nBQ*nAP+1
    if (ipz > memx) then
      write(u6,*) 'TraMO_MCLR: ipz > memx'
      write(u6,*) 'ipz,memx=',ipz,memx
      call Abend()
    end if
    Buffer(ipX-1) = -99999.0_wp
    Buffer(ipY-1) = -99999.0_wp
    Buffer(ipZ-1) = -99999.0_wp

    IAD2 = 0
    IAD3 = 0
    IVX = 0
    do NS=1,nAS
      nr2 = (NS-1)*NBR*NBP*NAQ
      nr3 = (NS-1)*nBR*NBQ*NAP
      if (.not. NoFile) then
        do NR=1,nBR
          IPQST = 1+NBPQ*IVX
          IVX = IVX+1
          IBUF = IVX

          !*************************************************************
          !
          !    Read one buffer of integrals back into core if INCORE=1
          !
          !*************************************************************

          if (INCORE == 1) then
            ipq = 0
            ip2 = 1
            ip5 = ip2+((nBPQ-1)/iMax+1)*iMax
            ipX = ip5+iMax+1
            ipY = ipX+nBP*nAQ*nBR*NAS+1
            ipZ = ipY+1
            if (isp /= isq) ipZ = ipX+nAQ*NBP*nBR*NAS+1
            Buffer(ipX-1) = -99999.0_wp
            Buffer(ipY-1) = -99999.0_wp
            Buffer(ipZ-1) = -99999.0_wp
            ipq = 1
            do
              IAD2 = IDAHLF2(IBUF)
              LPKREC = IRLHLF2(IBUF)
              call c_f_pointer(c_loc(Buffer(ip5)),iBuffer,[LPKREC])
              call iDAFILE(LUHLF2,2,iBuffer,LPKREC,IAD2)
              nullify(iBuffer)
              call UPKR8(0,iMax,NBYTES,Buffer(ip5),Buffer(ip2+IPQ-1))
              IPQ = IPQ+iMax
              iBuf = iBuf+nAS*nBR
              if (IPQ > NBPQ) exit
            end do
            IPQST = 1
          end if

          !*************************************************************
          !                                                            *
          ! Transform (ij|kL) -> (Ij|kL) & (iJ|kL)                     *
          !                                                            *
          !*************************************************************

          if (ISP == ISQ) then
            if (NAQ /= 0) then
              call SQUARE(Buffer(ip2+IPQST-1),X1,1,NBP,NBQ)
              call DGEMM_('T','N',NBP,NAQ,NBQ,One,X1,NBQ,CMQ,NBQ,Zero,X2,NBP)
            end if
          else
            if (NAQ /= 0) call DGEMM_('T','N',NBP,NAQ,NBQ,One,Buffer(ip2+IPQST-1),NBQ,CMQ,NBQ,Zero,X2,NBP)
            if (NAP /= 0) call DGEMM_('N','N',NBQ,NAP,NBP,One,Buffer(ip2+IPQST-1),NBQ,CMP,NBP,Zero,X3,NBQ)
          end if

          !*************************************************************
          !
          !      Move integrals to output buffer
          !
          !*************************************************************

          do i=0,nAQ-1
            Buffer(ipX+nr2+i*nBR*nBP:ipX+nr2+(i*nBR+1)*nBP-1) = X2(i*nBP+1:(i+1)*nBP)
          end do
          nr2 = nr2+nBP

          if (isp /= isq) then
            do i=0,nAP-1
              Buffer(ipY+nr3+i*nBR*nBQ:ipY+nr3+(i*nBR+1)*nBQ-1) = X3(i*nBQ+1:(i+1)*nBQ)
            end do
            nr3 = nr3+nBQ
          end if

        end do

        !***************************************************************
        !
        !       Empty last buffers
        !
        !***************************************************************

      end if ! nofile
    end do
    if (nAQ /= 0) then
      !call GADSum(Buffer(ipX),nAQ*NAS*NBR*NBP)
      call dDafile(LUTRI2,ione,Buffer(ipX),nAQ*NAS*NBR*NBP,iAD24)
    end if
    if (iSP /= iSQ .and. nAP /= 0) then
      !call GADSum(Buffer(ipY),NAP*NAS*NBR*NBQ)
      call dDafile(LUTRI3,ione,Buffer(ipY),NAP*NAS*NBR*NBQ,iAD14)
    end if

    !*******************************************************************

    if (buffer(ipX-1) /= -99999.0_wp) then
      write(u6,*) 'TraMO_MCLR: buffer(ipX-1) /= -99999.0'
      write(u6,*) 'buffer(ipX-1)=',buffer(ipX-1)
      call Abend()
    end if
    if (buffer(ipY-1) /= -99999.0_wp) then
      write(u6,*) 'TraMO_MCLR: buffer(ipY-1) /= -99999.0'
      write(u6,*) 'buffer(ipY-1)=',buffer(ipY-1)
      call Abend()
    end if
    if (buffer(ipZ-1) /= -99999.0_wp) then
      write(u6,*) 'TraMO_MCLR: buffer(ipZ-1) /= -99999.0'
      write(u6,*) 'buffer(ipZ-1)=',buffer(ipZ-1)
      call Abend()
    end if
    if (iSS /= iSR) then
      ipX = ip4
      Buffer(ipX-1) = -99999.0_wp
      ipY = ipX+nAR*NBS*nBP*nAQ+1   ! (Ij|kL)
      Buffer(ipY-1) = -99999.0_wp
      ipZ = ipy+1
      if (isp /= isq) ipZ = ipY+nAR*NBS*nBQ*nAP+1
      Buffer(ipZ-1) = -99999.0_wp
      IVX = 0
      do NR=1,nAR
        nr2 = NBP*NBS*(nR-1)*NAQ
        nr3 = NBQ*NBS*(nR-1)*NAP
        if (.not. NoFile) then
          do NS=1,nBS
            IPQST = 1+NBPQ*IVX
            IVX = IVX+1
            IBUF = IVX

            !***********************************************************
            !
            !   Read one buffer of integrals back into core if INCORE=1
            !
            !***********************************************************

            if (INCORE == 1) then
              ip2 = 1
              ip3 = ip2
              ip5 = ip3+((nBPQ-1)/iMax+1)*iMax
              ipX = ip5+iMax+1
              ipY = ipX+nBP*nAQ*nAR*NBS+1
              ipZ = ipY+nBQ*nAP*NAR*NBS+1
              Buffer(ipX-1) = -99999.0_wp
              Buffer(ipY-1) = -99999.0_wp
              Buffer(ipZ-1) = -99999.0_wp
              IPQ = 1
              do
                IAD3 = IDAHLF3(IBUF)
                LPKREC = IRLHLF3(IBUF)
                call c_f_pointer(c_loc(Buffer(ip5)),iBuffer,[LPKREC])
                call iDAFILE(LUHLF3,2,iBuffer,LPKREC,IAD3)
                nullify(iBuffer)
                call UPKR8(0,iMax,NBYTES,Buffer(ip5),Buffer(ip3+IPQ-1))
                IPQ = IPQ+iMax
                iBuf = iBuf+nAR*nBS
                if (IPQ > NBPQ) exit
              end do
              IPQST = 1
            end if

            !***********************************************************
            !                                                          *
            ! Transform (ij|kL) -> (Ij|kL) & (iJ|kL)                   *
            !                                                          *
            !***********************************************************

            if (ISP == ISQ) then
              if (nBQ*nAP /= 0) then

                call SQUARE(Buffer(ip3+IPQST-1),X1,1,NBQ,NBQ)
                call DGEMM_('T','N',NBP,NAQ,NBQ,One,X1,NBQ,CMQ,NBQ,Zero,X2,NBP)
                !call xxDGEMUL(X1,NBQ,'N',CMP,NBP,'N',X3,NBQ,NBQ,NBP,NAP)
                !call RecPrt('PqRs',' ',X3,nBQ,nAP)
              end if
            else
              if (nBP*nAQ /= 0) call DGEMM_('T','N',NBP,NAQ,NBQ,One,Buffer(ip3+IPQST-1),NBQ,CMQ,NBQ,Zero,X2,NBP)
              if (nBQ*nAP /= 0) call DGEMM_('N','N',NBQ,NAP,NBP,One,Buffer(ip3+IPQST-1),NBQ,CMP,NBP,Zero,X3,NBQ)
            end if

            !***********************************************************
            !
            !      Move integrals to output buffer
            !
            !***********************************************************

            do i=0,nAQ-1
              Buffer(ipX+nr2+i*nBS*nBP:ipX+nr2+(i*nBS+1)*nBP-1) = X2(i*nBP+1:(i+1)*nBP)
            end do
            nr2 = nr2+nBP
            if (iSP /= iSQ) then
              do i=0,nAP-1
                Buffer(ipY+nr3+i*NBS*nBQ:ipY+nr3+(i*NBS-1)*nBQ-1) = X3(i*nBQ+1:(i+1)*nBQ)
              end do
              nr3 = nr3+nBQ
            end if

          end do
          if (buffer(ipX-1) /= -99999.0_wp) then
            write(u6,*) 'TraMO_MCLR: buffer(ipX-1) /= -99999.0'
            write(u6,*) 'buffer(ipX-1)=',buffer(ipX-1)
            call Abend()
          end if
          if (buffer(ipY-1) /= -99999.0_wp) then
            write(u6,*) 'TraMO_MCLR: buffer(ipY-1) /= -99999.0'
            write(u6,*) 'buffer(ipY-1)=',buffer(ipY-1)
            call Abend()
          end if
          if (buffer(ipZ-1) /= -99999.0_wp) then
            write(u6,*) 'TraMO_MCLR: buffer(ipZ-1) /= -99999.0'
            write(u6,*) 'buffer(ipZ-1)=',buffer(ipZ-1)
            call Abend()
          end if

          !*************************************************************
          !
          !      Write buffer to disk (iJ|Kl) & (Ij|Kl)
          !
          !*************************************************************

        end if ! nofile
      end do
      if (nAQ /= 0) then
        !call GADSum(Buffer(ipX),nAR*nAQ*nBS*nBP)
        call dDaFile(LuTri4,ione,Buffer(ipX),nAR*nAQ*nBS*nBP,iad23)
      end if
      if (iSP /= iSQ .and. nap /= 0) then
        !call GADSum(Buffer(ipY),nAR*nAP*nBS*nBQ)
        call dDaFile(LuTri5,ione,Buffer(ipY),nAR*nAP*nBS*nBQ,iad13)
      end if
    end if
  end if

  if (INCORE == 1) then
    call DACLOS(LUHLF2)
    call DACLOS(LUHLF3)
  end if

  return

end subroutine TRAMO_MCLR_INTERNAL

end subroutine TRAMO_MCLR
