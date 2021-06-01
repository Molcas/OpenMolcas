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
! Copyright (C) 2007, Francesco Aquilante                              *
!***********************************************************************

subroutine Cho_SOSmp2_Energy(irc,EMP2,EOcc,EVir,Delete)
! Francesco Aquilante, May 2007.
!
! Purpose: compute "Scaled Opposite-Spin" MP2 energy correction from
!          MO Cholesky vectors of the matrix M(ai,bj)=(ai|bj)^2.

implicit real*8(a-h,o-z)
real*8 EMP2
real*8 EOcc(*), EVir(*)
integer irc
logical Delete
character*6 ThisNm
character*17 SecNam
parameter(SecNam='Cho_SOSmp2_Energy',ThisNm='Energy')
parameter(zero=0.0d0,one=1.0d0,two=2.0d0)
integer ipWrk, lWrk, iiSoff(8), iaSoff(8)
integer nEnrVec(8)
#include "cholesky.fh"
#include "chomp2.fh"
#include "chomp2_cfg.fh"
#include "WrkSpc.fh"
!****************************************************************
MulD2h(i,j) = ieor(i-1,j-1)+1
!****************************************************************

irc = 0

iTyp = 2
call iCopy(nSym,nMP2Vec,1,nEnrVec,1)

! Initialize SOS-MP2 energy correction.
! -------------------------------------

EMP2 = 0.0d0

! Some offsets
! ------------
nIt = nOcc(1)
nAt = nVir(1)
iiSoff(1) = 0
iaSoff(1) = 0
do iSym=2,nSym
  iiSoff(iSym) = nIt
  iaSoff(iSym) = nAt
  nIt = nIt+nOcc(iSym)
  nAt = nAt+nVir(iSym)
end do
MaxNVec = nIt*nAt

call GetMem('Ea-Ei','Allo','Real',ip_W,MaxNVec)
call GetMem('iD_bj','Allo','Inte',ID_bj,MaxNVec)

! Loop over Cholesky vector symmetries.
! -------------------------------------

do jSym=1,nSym

  Nai = nT1am(jSym)
  if ((Nai > 0) .and. (nEnrVec(jSym) > 0)) then

    nOV = 0
    do iiSym=1,nSym
      iaSym = MulD2h(iiSym,jSym)
      do ii=1,nOcc(iiSym)
        iiT = iiSoff(iiSym)+ii
        iaS = nOV+nVir(iaSym)*(ii-1)
        do ia=1,nVir(iaSym)
          iaT = iaSoff(iaSym)+ia
          iaiS = iaS+ia-1
          Work(ip_W+iaiS) = EVir(iaT)-EOcc(iiT)
        end do
      end do
      nOV = nOV+nVir(iaSym)*nOcc(iiSym) ! ... = Nai
    end do

    ! Cholesky decompsition of the Orbital Energy
    ! Denominators (OED)
    ! -------------------------------------------
    call CHO_GET_ORD_bj(nOV,MaxNVec,OED_Thr,Work(ip_W),iWork(ID_bj),NKVec,Dmax)

    if (Verbose .or. (NKVec < 1)) then
      write(6,'(A)') '---------------------------------------'
      write(6,'(A,I2,A)') 'Orbital energy denominators CD (sym=',jSym,')'
      write(6,'(A)') '---------------------------------------'
      write(6,'(1X,A,I3,A,I9,A,1P,D25.16)') 'Number of vectors needed: ',NKVec,'   ( nAocc x nAvir : ',nOV,' ), max residual:',Dmax
      call xFlush(6)
    end if

    if (NKVec > 0) then

      call GetMem('Yai_k','Allo','Real',ip_Y,nOV*NKVec)
      ! init to one the 1st col
      call dcopy_(nOV,[one],0,Work(ip_Y),1)

      call CHO_GET_OED_cd(.true.,nOV,Work(ip_W),NKVec,iWork(ID_bj),1,Work(ip_Y),Work(ip_Y))

      call GetMem('GetMax','Max ','Real',ipWrk,lWrk)
      call GetMem('GetMax','Allo','Real',ipWrk,lWrk)

      ! Set up batch over Cholesky vectors.
      ! -----------------------------------

      nVec = min(lWrk/(Nai+NKVec),nEnrVec(jSym))
      if (nVec < 1) then
        call ChoMP2_Quit(SecNam,'Insufficient memory','Batch setup')
      end if
      nBat = (nEnrVec(jSym)-1)/nVec+1

      kRead = ipWrk+NKVec*nVec

      ! Open Cholesky vector files.
      ! ---------------------------

      call ChoMP2_OpenF(1,iTyp,jSym)

      ! Start vector batch loop.
      ! ------------------------

      do iBat=1,nBat

        if (iBat == nBat) then
          NumVec = nEnrVec(jSym)-nVec*(nBat-1)
        else
          NumVec = nVec
        end if
        jVec = nVec*(iBat-1)+1

        ! Read vectors
        ! ------------
        lTot = nT1am(jSym)*NumVec
        iAdr = nT1am(jSym)*(jVec-1)+1
        call ddaFile(lUnit_F(jSym,iTyp),2,Work(kRead),lTot,iAdr)

        ! Compute   E(k,J) = sum_ai Y(ai,k) * R(ai,J)
        ! --------------------------------------------
        call DGEMM_('T','N',NKVec,NumVec,Nai,one,Work(ip_Y),Nai,Work(kRead),Nai,zero,Work(ipWrk),NKVec)

        ! Compute (unscaled) SOS-MP2 energy
        ! -----------------------------------

        EMP2 = EMP2+ddot_(NKVec*NumVec,Work(ipWrk),1,Work(ipWrk),1)

      end do ! Cholesky vector batch

      ! Close Cholesky vector files.
      ! ----------------------------
      call ChoMP2_OpenF(2,iTyp,jSym)

      call GetMem('GetMax','Free','Real',ipWrk,lWrk)
      call GetMem('Yai_k','Free','Real',ip_Y,nOV*NKVec)

    end if

  end if

end do

call GetMem('iD_bj','Free','Inte',ID_bj,MaxNVec)
call GetMem('Ea-Ei','Free','Real',ip_W,MaxNVec)

! If requested, delete vector files.
! ----------------------------------

if (Delete) then
  do jSym=1,nSym
    call ChoMP2_OpenF(1,iTyp,jSym)
    call ChoMP2_OpenF(3,iTyp,jSym)
  end do
end if

! Change sign and use proper factor on energy.
! --------------------------------------------

!-tbp, December 2012: removed factor 2
!-tbp  (wrong result for two-electron systems with factor 2)
!-tbp EMP2 = -two*EMP2
EMP2 = -EMP2

return

end subroutine Cho_SOSmp2_Energy
