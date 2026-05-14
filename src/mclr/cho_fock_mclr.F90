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

subroutine CHO_Fock_MCLR(DA,G2,W_JA,W_KA,W_FkA,CVa,W_CMO,nIsh,nAsh,LuAChoVec)
!***********************************************************************
!                                                                      *
!  Author : M. G. Delcey                                               *
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri
use Symmetry_Info, only: Mul
use Cholesky, only: InfVec, nBas, nBasSh, nDimRS, nShell, nSym, NumCho
use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: DA(*), G2(*), W_CMO(*)
real(kind=wp), intent(_OUT_) :: W_JA(*), W_KA(*), W_FkA(*)
type(DSBA_type), intent(in) :: CVa
integer(kind=iwp), intent(in) :: nIsh(8), nAsh(8), LuAChoVec(8)
#include "warnings.h"
integer(kind=iwp) :: i, ia, iabg, iAdr, iag, iaSh, ib, iBatch, ibg, ibSh, iLoc, ioff, ioffa, ioffb, ipLpq(8,3), ipLtvb, ipLvb, &
                     ipLvtw, ipLvw, ipLwb, ipLxy, ipVJ, irc, IREDC, is, iSym, iSyma, iSymb, iSymv, IVEC2, iVrs, JNUM, JRED, JRED1, &
                     JRED2, jS, jSym, JVC, JVEC, k, l, lChoa, LREAD, lvec, LWORK, MaxVecPerBatch, mTvec, MUSED, NAv, NAw, nBatch, &
                     nMat, nRS, NumCV, NUMV, nVec, nVrs
logical(kind=iwp) :: add
type(DSBA_type) :: JA(1), KA, Fka, CMO, Scr
integer(kind=iwp), allocatable :: kOffSh(:,:)
real(kind=wp), allocatable :: Fab(:), LF(:), Lrs(:)
real(kind=wp), parameter :: FactCI = -Two, FactXI = Half
integer(kind=iwp), external :: Cho_LK_MaxVecPerBatch

call Allocate_DT(JA(1),nBas,nBas,nSym,aCase='TRI',Ref=W_JA)
call Allocate_DT(KA,nBas,nBas,nSym,Ref=W_KA)
call Allocate_DT(FkA,nBas,nBas,nSym,Ref=W_FkA)
call Allocate_DT(CMO,nBas,nBas,nSym,Ref=W_CMO)

! Compute Shell Offsets (MOs and transformed vectors)

call mma_allocate(kOffSh,nShell,nSym,Label='kOffSh')
do iSyma=1,nSym
  kOffSh(1,iSyma) = 0
  do iaSh=2,nShell ! kOffSh(iSh,iSym)
    kOffSh(iaSh,iSyma) = kOffSh(iaSh-1,iSyma)+nBasSh(iSyma,iaSh-1)
  end do
end do

! memory for the Q matrices --- temporary array
call Allocate_DT(Scr,nBas,nBas,nSym)
Scr%A0(:) = Zero

MaxVecPerBatch = Cho_LK_MaxVecPerBatch()

! Start looping!

do jSym=1,nSym
  NumCV = NumCho(jSym)
  call GAIGOP_SCAL(NumCV,'max')
  if (NumCV < 1) cycle

  iLoc = 3 ! use scratch location in reduced index arrays

  ! Estimate memory need

  mTvec = 0
  do l=1,nSym
    mTvec = mTvec+nAsh(Mul(l,JSYM))*nBas(l)*3
  end do

  JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
  JRED2 = InfVec(NumCho(jSym),2,jSym) ! red set of the last vec

  ! entire red sets range for parallel run
  call GAIGOP_SCAL(JRED1,'min')
  call GAIGOP_SCAL(JRED2,'max')

  do JRED=JRED1,JRED2
    call Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)
    if (nVrs == 0) cycle
    if (nVrs < 0) then
      write(u6,*) 'CHO_FOCK_MCLR: Cho_X_nVecRS returned nVrs<0. STOP!',nVrs
      call Abend()
    end if
    call Cho_X_SetRed(irc,iLoc,JRED)
    ! set index arrays at iLoc
    if (irc /= 0) then
      write(u6,*) 'CHO_FOCK_MCLR: cho_X_setred non-zero return code. rc= ',irc
      call Abend()
    end if
    IREDC = JRED

    nRS = nDimRS(JSYM,JRED)
    if (jSym == 1) then
      call mma_allocate(Fab,nRS,Label='Fab')
      Fab(:) = Zero
    end if

    call mma_MaxDBLE(LWORK)
    nVec = min(LWORK/(nRS+mTvec+1),min(nVrs,MaxVecPerBatch))
    if (nVec < 1) then
      write(u6,*) 'CHO_FOCK_MCLR: Insufficient memory for batch'
      write(u6,*) 'LWORK= ',LWORK
      write(u6,*) 'min. mem. need= ',nRS+mTvec+1
      write(u6,*) 'nRS= ',nRS
      write(u6,*) 'mTvec= ',mTvec
      write(u6,*) 'jsym= ',jsym
      call Quit(_RC_MEMORY_ERROR_)
      nBatch = -9999  ! dummy assignment
    end if
    LREAD = nRS*nVec

    call mma_allocate(Lrs,LREAD,Label='Lrs')
    call mma_allocate(LF,mTvec*nVec,Label='LF')

    nBatch = (nVrs-1)/nVec+1

    do iBatch=1,nBatch

      if (iBatch == nBatch) then
        JNUM = nVrs-nVec*(nBatch-1)
      else
        JNUM = nVec
      end if
      !*****************************************************************
      !                                                                *
      !     START WORKING                                              *
      !                                                                *
      !*****************************************************************

      ! Read Cholesky vector

      JVEC = nVec*(iBatch-1)+iVrs
      IVEC2 = JVEC-1+JNUM

      call CHO_VECRD(Lrs,LREAD,JVEC,IVEC2,JSYM,NUMV,IREDC,MUSED)

      if ((NUMV <= 0) .or. (NUMV /= JNUM)) return

      !call CWTIME(TCINT1,TWINT1)

      ! Set up the skipping flags and the pointers ipLpq
      ! The memory used before for the full-dimension AO-vectors
      !     is now re-used to store half and full transformed
      !     vectors in the active space
      !-----------------------------------------------------------------
      lChoa = 0
      do i=1,nSym

        k = Mul(i,JSYM)

        ipLpq(k,1) = 1+lChoa   ! Lvb,J
        ipLpq(k,2) = ipLpq(k,1)+nAsh(k)*nBas(i)*JNUM ! Lvi,J i general MO index
        ipLpq(k,3) = ipLpq(k,2)+nAsh(k)*nBas(i)*JNUM ! L~vi,J ~ transformed index

        lChoa = lChoa+nAsh(k)*nBas(i)*3*JNUM

      end do

      !  Read half-transformed cho vectors

      ioff = 0
      do i=1,nSym
        k = Mul(i,JSYM)
        lvec = nAsh(k)*nBas(i)*JNUM
        iAdr = (JVEC-1)*nAsh(k)*nBas(i)+ioff
        call DDAFILE(LuAChoVec(Jsym),2,LF(ipLpq(k,1)),lvec,iAdr)
        ioff = ioff+nAsh(k)*nBas(i)*NumCho(jSym)
      end do
      !call CWTIME(TCINT2,TWINT2)
      !tint1(1) = tint1(1)+(TCINT2-TCINT1)
      !tint1(2) = tint1(2)+(TWINT2-TWINT1)
      !-----------------------------------------------------------------
      ! Active-Active transformation  Lvw,J = sum_b  Lvb,J * C2(w,b)
      !-----------------------------------------------------------------
      do iSymb=1,nSym

        iSymv = Mul(JSYM,iSymb)
        NAv = nAsh(iSymv)
        NAw = nAsh(iSymb)

        if (NAv*NBAS(iSymb) /= 0) then

          !call CWTIME(TCINT2,TWINT2)

          do JVC=1,JNUM
            ipLvb = ipLpq(iSymv,1)+NAv*NBAS(iSymb)*(JVC-1)
            ipLvw = ipLpq(iSymv,2)+NAv*Naw*(JVC-1)
            call DGEMM_('N','T',NAv,Naw,NBAS(iSymb),One,LF(ipLvb),NAv,CVa%SB(iSymb)%A2,Naw,Zero,LF(ipLvw),NAv)
          end do
          !call CWTIME(TCINT2,TWINT2)
          !tint1(1) = tint1(1)+(TCINT2-TCINT3)
          !tint1(2) = tint1(2)+(TWINT2-TWINT3)

          ! ********** EVALUATION OF THE ACTIVE FOCK MATRIX ************
          ! Coulomb term
          ipVJ = ipLpq(iSymv,3)
          ipLvtw = ipLpq(iSymv,2)

          call DGEMV_('T',Nav*Naw,JNUM,One,LF(ipLvtw),Nav*Naw,DA,1,Zero,LF(ipVJ),1)

          call DGEMV_('N',nRS,JNUM,-FactCI,Lrs,nRS,LF(ipVJ),1,One,Fab,1)

          !call CWTIME(TCINT2,TWINT2)
          !tact(1) = tact(1)+(TCINT2-TCINT3)
          !tact(2) = tact(2)+(TWINT2-TWINT3)
          !-------------------------------------------------------------
          ! Formation of the Q matrix Qpx = L~py Lvw Gxyvw
          !-------------------------------------------------------------
          do JVC=1,JNUM
            ! Lxy=Lvw Gxyvw
            !MGD probably additional nSym loop
            ipLvb = ipLpq(iSymv,1)+NAv*NBAS(iSymb)*(JVC-1)
            ipLvw = ipLpq(iSymv,2)+NAv*Naw*(JVC-1)
            ipLxy = ipLpq(iSymv,3)+NAv*Naw*(JVC-1)
            call DGEMV_('N',NAv*Naw,NAv*Naw,One,G2,NAv*Naw,LF(ipLvw),1,Zero,LF(ipLxy),1)
            ! Qpx=Lpy Lxy
            call DGEMM_('T','N',NBAS(iSymb),NAw,Nav,One,LF(ipLvb),NAv,LF(ipLxy),Naw,One,Scr%SB(iSymb)%A2,NBAS(iSymb))
          end do
          !call CWTIME(TCINT3,TWINT3)
          !tQmat(1) = tQmat(1)+(TCINT3-TCINT2)
          !tQmat(2) = tQmat(2)+(TWINT3-TWINT2)

          ! ********** EVALUATION OF THE ACTIVE FOCK MATRIX ************
          ! Exchange term
          do JVC=1,JNUM
            ipLvb = ipLpq(iSymv,1)+NAv*NBAS(iSymb)*(JVC-1)
            ipLwb = ipLpq(iSymv,2)+NAv*NBAS(iSymb)*(JVC-1)
            call DGEMM_('T','N',NBAS(iSymb),Nav,Nav,One,LF(ipLvb),Nav,DA,Nav,Zero,LF(ipLwb),NBAS(iSymb))
          end do
          !call CWTIME(TCINT2,TWINT2)
          !tact(1) = tact(1)+(TCINT2-TCINT3)
          !tact(2) = tact(2)+(TWINT2-TWINT3)
          do JVC=1,JNUM
            ipLwb = ipLpq(iSymv,2)+NAv*NBAS(iSymb)*(JVC-1)
            do is=1,NBAS(iSymb)
              ipLtvb = ipLpq(iSymv,1)+NAv*NBAS(iSymb)*(JVC-1)+Nav*(is-1)
              call DGEMV_('N',NBAS(iSymb),Nav,-FactXI,LF(ipLwb),NBAS(iSymb),LF(ipLtvb),1,One,KA%SB(iSymb)%A2(:,is),1)

            end do
          end do
          !call CWTIME(TCINT3,TWINT3)
          !tact(1) = tact(1)+(TCINT3-TCINT2)
          !tact(2) = tact(2)+(TWINT3-TWINT2)
        end if
      end do

      ! All good things come to an end

    end do  ! end batch loop

    ! backtransform fock matrix to full storage
    if (JSYM == 1) then
      add = .true.
      nMat = 1
      call swap_rs2full(irc,iLoc,nRS,nMat,JSYM,JA,Fab,add)
      call mma_deallocate(Fab)
    end if
    call mma_deallocate(Lrs)
    call mma_deallocate(LF)
  end do  ! loop over red sets
end do    ! loop over JSYM
!***********************************************************************
!                                                                      *
!     POST PROCESSING                                                  *
!                                                                      *
!***********************************************************************
!
! Accumulate Coulomb and Exchange contributions

do iSym=1,nSym

  do iaSh=1,nShell
    ioffa = kOffSh(iaSh,iSym)

    do ibSh=1,nShell

      ioffb = kOffSh(ibSh,iSym)

      do ib=1,nBasSh(iSym,ibSh)
        do ia=1,nBasSh(iSym,iaSh)
          !MGD warning with sym

          iag = ioffa+ia
          ibg = ioffb+ib
          iabg = iTri(iag,ibg)

          FkA%SB(iSym)%A2(iag,ibg) = JA(1)%SB(iSym)%A1(iabg)+KA%SB(iSym)%A2(iag,ibg)+KA%SB(iSym)%A2(ibg,iag)
        end do
      end do
    end do
  end do
end do

call Deallocate_DT(JA(1))
call Allocate_DT(JA(1),nBas,nBas,nSym,Ref=W_JA)

! Transform Fock and Q matrix to MO basis

do iS=1,nSym
  jS = iS
  if (nBas(iS) /= 0) then
    call DGEMM_('T','N',nBas(jS),nBas(iS),nBas(iS),One,FkA%SB(iS)%A2,nBas(iS),CMO%SB(iS)%A2,nBas(iS),Zero,JA(1)%SB(iS)%A2,nBas(jS))
    FkA%SB(is)%A2(:,:) = Zero
    call DGEMM_('T','N',nBas(jS),nIsh(jS),nBas(iS),One,JA(1)%SB(iS)%A2,nBas(iS),CMO%SB(jS)%A2,nBas(jS),Zero,FkA%SB(iS)%A2,nBas(jS))
    !ioff = nIsh(iS)+1
    iOff = 1+nIsh(iS)*nBas(iS)
    call DGEMM_('T','N',nBas(jS),nAsh(iS),nBas(jS),One,CMO%SB(iS)%A2,nBas(jS),Scr%SB(iS)%A2,nBas(jS),Zero,FkA%SB(iS)%A1(iOff:), &
                nBas(jS))
  end if
end do
!***********************************************************************
!                                                                      *
!     TERMINATING                                                      *
!                                                                      *
!***********************************************************************
call deallocate_DT(Scr)
call mma_deallocate(kOffSh)
call Deallocate_DT(CMO)
call Deallocate_DT(FkA)
call Deallocate_DT(KA)
call Deallocate_DT(JA(1))

return

end subroutine CHO_Fock_MCLR
