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

subroutine Mult_RijK_QKL(iSO,nBas_aux,nIrrep)
!*************************************************************************
!     Author: Jonas Bostrom
!
!     Purpose: Computation of the DF coefficient vectors in MO-basis.
!
!     Equation:   C(il,K)  =  sum_J  R(il,J) * Q(K,J)
!
!     Input:
!            iSO : alpha (iSO=1) or beta (iSO=2) orbitals
!            nBas_aux : number of aux bsfs in each irrep.
!            nIrrep : number of irreps.
!*************************************************************************

use RI_glob, only: iAdrCVec, nChOrb, nIJ1, nIJR, NumAuxVec
use Symmetry_Info, only: Mul
use pso_stuff, only: lSA
use Para_Info, only: Is_Real_Par
use Cholesky, only: InfVec, NumCho, timings
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iSO, nIrrep, nBas_Aux(1:nIrrep)
!#define _CD_TIMING_
#ifdef _CD_TIMING_
#include "temptime.fh"
#endif
integer(kind=iwp) :: i, iAdrC, iAdrQ, iAdrR, iAux, iFirstCho, iJBat, indx, indx2, iRest, iSeed, iSeed2, iSym, jSym, kSym, &
                     l_CVector, l_Q, l_QVector, l_RVec, l_RVector, lSym, Lu_Q, LuCVec, LuRVec, MaxCho, MaxLocCho, MaxMOprod, &
                     MaxMOProdR, MemMax, nJbat, njVec, nJvec1, nJvecLast, nTotCho, nTotFIorb, NumAux, NumCV, nVec(1:nIrrep)
real(kind=wp) :: TotCPU, TotCPU1, TotCPU2, TotWall, TotWall1, TotWall2
character(len=50) :: CFmt
character(len=6) :: Fname, Fname2, Name_Q
real(kind=wp), allocatable :: CVector(:), QVector(:), RVector(:)
character(len=*), parameter :: SECNAM = 'MULT_RIJK_QKL'
integer(kind=iwp), external :: IsFreeUnit

call CWTime(TotCPU1,TotWall1)

nVec(1:nIrrep) = NumCho(1:nIrrep)
call GAIGOP(nVec,nIrrep,'+')

nTotCho = 0
MaxCho = 0
MaxLocCho = 0
do jSym=1,nIrrep
  nTotCho = nTotCho+nVec(jSym)
  MaxCho = max(MaxCho,nVec(jSym))
  MaxLocCho = max(MaxLocCho,NumCho(jSym))
end do
iAdrCVec(1:nIrrep,1:nIrrep,iSO) = 0

! Loop over the first cholesky symmetry

do jSym=1,nIrrep

  ! Check so the symmetry contains vectors
  !---------------------------------------
  NumCV = NumCho(jSym)
  NumAux = nBas_Aux(jSym)
  if (jSym == 1) NumAux = NumAux-1
  ! Save the number of auxiliary basis functions to be accessed later
  NumAuxVec(jSym) = NumAux

  call GAIGOP_SCAL(NumCV,'max')
  if (NumCV < 1) cycle

  nTotFIorb = 0
  MaxMOprod = 0
  MaxMOprodR = 0
  do iSym=1,nIrrep
    kSym = Mul(JSym,iSym)

    nTotFIorb = nTotFIorb+nIJ1(iSym,kSym,iSO)
    MaxMOprod = max(MaxMOprod,nIJ1(iSym,kSym,iSO))
    MaxMOProdR = max(MaxMOprodR,nIJR(iSym,kSym,iSO))
  end do

  call mma_maxDBLE(MemMax)
  nJvec1 = (MemMax-NumAux*MaxMOprod)/(MaxMOprod+NumAux)
  if (nJvec1 < 1) then
    write(u6,*) 'Too little memory in:',SECNAM
    call Abend()
  end if
  nJvec1 = min(nJvec1,NumCho(jSym))
  nJbat = NumCho(jSym)/nJvec1
  iRest = mod(NumCho(jSym),nJvec1)
  if (iRest /= 0) then
    nJbat = nJbat+1
    nJvecLast = iRest
  else
    nJvecLast = nJvec1
  end if

  if (Is_Real_Par()) then
    iFirstCho = InfVec(1,5,jSym)
  else
    iFirstCho = 1
  end if

  l_QVector = nJVec1*NumAux
  l_RVector = MaxMOprod*nJVec1
  l_CVector = MaxMOprodR*NumAux

  call mma_allocate(QVector,l_QVector,Label='QVector')
  call mma_allocate(RVector,l_RVector,Label='RVector')
  call mma_allocate(CVector,l_CVector,Label='CVector')

  iSeed2 = 8
  LuCVec = IsFreeUnit(iSeed2)
  if (iSO == 1) then
    write(Fname2,'(A4,I1,I1)') 'CVEA',jSym
  else if (iSO == 2) then
    write(Fname2,'(A4,I1,I1)') 'CVEB',jSym
  end if
  call DANAME_MF_WA(LuCVec,Fname2)
  iAdrC = 0

  ! Get Q Vectors from Disk
  !------------------------
  do iSym=1,nIrrep
    lSym = Mul(iSym,jSym)

    if (nIJ1(iSym,lSym,iSO) >= 1) then
      CVector(:) = Zero

      do iJBat=1,nJBat
        if (iJBat == nJBat) then
          njVec = nJVecLast
        else
          nJvec = nJvec1
        end if

        iSeed = 55+jSym-1
        Lu_Q = IsFreeUnit(iSeed)
        write(Name_Q,'(A4,I2.2)') 'QVEC',jSym-1
        call DaName_MF_WA(Lu_Q,Name_Q)
        l_Q = nJvec*NumAux
        iAdrQ = (iFirstCho-1)*NumAux+(iJBat-1)*nJVec*NumAux
        call dDaFile(Lu_Q,2,Qvector,l_Q,iAdrQ)

#       ifdef _DEBUGPRINT_
        call RecPrt('Q-vectors',' ',QVector,nJVec,NumAux)
#       endif

        iSeed = 7
        LuRVec = IsFreeUnit(iSeed)
        if (iSO == 1) then
          write(Fname,'(A4,I1,I1)') 'CHTA',iSym,lSym
        else if (iSO == 2) then
          write(Fname,'(A4,I1,I1)') 'CHTB',iSym,lSym
        end if
        call DANAME_MF_WA(LuRVec,Fname)

        ! Loop over all cholesky vectors on all nodes
        !--------------------------------------------

        ! Get R-Vectors from disk
        !------------------------

        iAdrR = nIJ1(iSym,lSym,iSO)*nJVec1*(iJBat-1)
        l_RVec = nJvec*nIJ1(iSym,lSym,iSO)
        call dDaFile(LuRVec,2,RVector,l_RVec,iAdrR)

        call dGemm_('N','T',nIJ1(iSym,lSym,iSO),NumAux,nJVec,One,RVector,nIJ1(iSym,lSym,iSO),QVector,NumAux,Zero,CVector, &
                    nIJ1(iSym,lSym,iSO))
      end do

#     ifdef _DEBUGPRINT_
      write(u6,*) 'jSym=',jSym
      call RecPrt('R-Vectors',' ',RVector,nIJ1(iSym,lSym,iSO),NumAux)
      call RecPrt('C-Vectors',' ',CVector,nIJ1(iSym,lSym,iSO),NumAux)
#     endif
      if ((.not. lSA) .and. (iSym == lSym)) then
        do iAux=1,NumAux
          indx = -1
          do i=1,nChOrb(iSym-1,iSO)
            indx = indx+i
            indx2 = indx+(iAux-1)*nIJ1(iSym,lSym,iSO)
            CVector(1+indx2) = CVector(1+indx2)/sqrt(Two)
          end do
        end do
      end if

      call DaClos(Lu_Q)
      call DACLOS(LuRVec)

    end if

    call GADGOP(CVector,l_CVector,'+')
    iAdrCVec(jSym,iSym,iSO) = iAdrC
    call dDaFile(LuCVec,1,CVector,nIJ1(iSym,lSym,iSO)*NumAux,iAdrC)
    if (nIJ1(iSym,lSym,iSO) < nIJR(iSym,lSym,iSO)) then
      iAdrC = iAdrC+(nIJR(iSym,lSym,iSO)-nIJ1(iSym,lSym,iSO))*NumAux
    end if

  end do !iSym

  call mma_deallocate(CVector)
  call mma_deallocate(RVector)
  call mma_deallocate(QVector)

  call DACLOS(LuCVec)

end do  ! jSym

call CWTime(TotCPU2,TotWall2)
TotCPU = TotCPU2-TotCPU1
TotWall = TotWall2-TotWall1
#ifdef _CD_TIMING_
rMult_CPU = TOTCPU
rMult_Wall = TOTWALL
#endif

if (timings) then

  CFmt = '(2x,A)'
  write(u6,*)
  write(u6,CFmt) 'Cholesky Gradients timing from '//SECNAM
  write(u6,CFmt) '----------------------------------------'
  write(u6,*)
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,CFmt) '                                CPU       WALL   '
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,'(2x,A26,2f10.2)') 'TOTAL                                     ',TOTCPU,TOTWALL
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,*)

end if

return

end subroutine Mult_RijK_QKL
