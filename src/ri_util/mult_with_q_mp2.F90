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
! Copyright (C) Jonas Bostrom                                          *
!***********************************************************************

subroutine Mult_with_Q_MP2(nBas_aux,nBas,nIrrep)
!***********************************************************************
! Author: Jonas Bostrom
!
! Purpose: Multiply MP2 A~_sep and B~_sep with inverse cholesky factors.
!
!***********************************************************************

use Symmetry_Info, only: Mul
use RI_glob, only: nAuxVe
use Cholesky, only: nSym, NumCho
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nIrrep, nBas_Aux(1:nIrrep), nBas(1:nIrrep)
integer(kind=iwp) :: i, iA_in, iA_Out, iAdr, iAdrA, iAdrA_in(8), iAdrA_Out(8), iAdrQ, iOffQ1, iOpt, ip_B, iSym, iType, jSym, jVec, &
                     kSym, kVec, l_A, l_A_ht, l_A_t, l_B_t, l_Q, lTot, Lu_A(2), Lu_B(4), Lu_Q, MaxMem, nBas2, nLR, nLRb(8), &
                     NumAux, NumCV, NumVecJ, NumVecK, nVec
real(kind=wp) :: Fac, TotCPU1, TotWall1
character(len=6) :: FName, Name_Q
real(kind=wp), allocatable :: A(:), A_ht(:), A_t(:), B_t(:), QVec(:)
character(len=*), parameter :: SECNAM = 'Mult_with_Q_MP2'
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TotCPU1,TotWall1)

do i=1,2
  Lu_A(i) = IsFreeUnit(7)
  write(FName,'(A5,I1)') 'AMP2V',i
  call DaName_MF_WA(Lu_A(i),FName)
end do

do i=1,4
  Lu_B(i) = IsFreeUnit(7)
  write(FName,'(A5,I1)') 'BMP2V',i
  call DaName_MF_WA(Lu_B(i),FName)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
iA_in = 1
iA_Out = 1
do iSym=1,nSym
  NumCV = NumCho(iSym)
  NumAux = nBas_Aux(iSym)-1
  nLR = 0
  do jSym=1,nSym
    kSym = Mul(iSym,jSym)
    nLR = nLR+nBas(jSym)*nBas(kSym)
  end do
  nLRb(iSym) = nLR
  iAdrA_in(iSym) = iA_in
  iA_in = iA_in+NumCV*NumCV
  iAdrA_out(iSym) = iA_out
  iA_out = iA_out+NumAux*NumAux
end do
!                                                                      *
!***********************************************************************
!                                                                      *
do iSym=1,nSym

  nBas2 = nLRb(iSym)
  NumCV = NumCho(iSym)
  NumAux = nBas_Aux(iSym)-1
  nAuxVe = NumAux

  ! Get Q-vectors from disk
  ! -----------------------

  l_Q = NumCV*NumAux
  call mma_allocate(QVec,l_Q,Label='QVec')

  Lu_Q = IsFreeUnit(7)
  write(Name_Q,'(A4,I2.2)') 'QVEC',iSym-1
  call DaName_MF_WA(Lu_Q,Name_Q)

  iOpt = 2
  iAdrQ = 0
  call dDaFile(Lu_Q,iOpt,QVec,l_Q,iAdrQ)

  ! Get MP2 A-tilde vectors from disk
  ! ---------------------------------

  l_A_t = NumCV*NumCV
  l_A = NumAux*NumAux
  l_A_ht = NumAux*NumCV
  call mma_allocate(A_t,l_A_t,Label='A_t')
  call mma_allocate(A,l_A,Label='A')
  call mma_allocate(A_ht,l_A_ht,Label='A_ht')

  do iType=1,2

    iOpt = 2
    iAdrA = iAdrA_in(iSym)
    call dDaFile(Lu_A(iType),iOpt,A_t,l_A_t,iAdrA)
#   ifdef _DEBUGPRINT_
    write(u6,*) 'Q-vectors'
    do i=1,l_Q
      write(u6,*) QVec(i)
    end do

    write(u6,*) 'A-vectors'
    do i=1,l_A_t
      write(u6,*) A_t(i)
    end do
#   endif

    ! Make first halftransformation to cholesky-base
    ! ----------------------------------------------

    call dGemm_('N','N',NumAux,NumCV,NumCV,One,QVec,NumAux,A_t,NumCV,Zero,A_ht,NumAux)

    call dGemm_('N','T',NumAux,NumAux,NumCV,One,A_ht,NumAux,QVec,NumAux,Zero,A,NumAux)

    ! Put transformed A-vectors back on disk

    iOpt = 1
    iAdrA = iAdrA_out(iSym)
    call dDaFile(Lu_A(iType),iOpt,A,l_A,iAdrA)

  end do

  call mma_deallocate(A_t)
  call mma_deallocate(A)
  call mma_deallocate(A_ht)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call mma_maxDBLE(MaxMem)
  MaxMem = 9*(MaxMem/10)
  call mma_allocate(B_t,MaxMem,Label='B_t')

  nVec = MaxMem/(2*nLRb(iSym))
  nVec = min(max(NumCV,NumAux),nVec)
  if (nVec < 1) call SysAbendMsg(SecNam,'nVec is non-positive','[1]')

  l_B_t = nLRb(iSym)*nVec
  ip_B = 1+l_B_t

  do iType=1,2

    ! The B-vectors should be read one batch at the time
    ! --------------------------------------------------

    do kVec=1,NumAux,nVec
      NumVecK = min(nVec,NumAux-(kVec-1))

      do jVec=1,NumCV,nVec
        NumVecJ = min(nVec,NumCV-(jVec-1))

        iOpt = 2
        lTot = NumVecJ*nLRb(iSym)
        iAdr = 1+nLRb(iSym)*(jVec-1)
        call dDaFile(Lu_B(iType),iOpt,B_t,lTot,iAdr)

        Fac = Zero
        if (jVec /= 1) Fac = One
        iOffQ1 = kVec+NumAux*(jVec-1)
        call dGemm_('N','T',nBas2,NumVecK,NumVecJ,One,B_t,nBas2,QVec(iOffQ1),NumAux,Fac,B_t(ip_B),nBas2)
      end do

      iOpt = 1
      lTot = NumVecK*nBas2
      iAdr = 1+nBas2*(kVec-1)
      call dDaFile(Lu_B(iType+2),iOpt,B_t(ip_B),lTot,iAdr)

    end do

  end do

  call mma_deallocate(B_t)
  call mma_deallocate(QVec)

  call DaClos(Lu_Q)

end do ! iSym

do i=1,2
  call DaClos(Lu_A(i))
end do
do i=1,4
  call DaClos(Lu_B(i))
end do

return

end subroutine Mult_with_Q_MP2
