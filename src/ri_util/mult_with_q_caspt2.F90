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

subroutine Mult_with_Q_CASPT2(nBas_aux,nBas,nIrrep,SubAux)
!***********************************************************************
! Author: Jonas Bostrom
!
! Purpose: Multiply MP2 A~_sep and B~_sep with inverse cholesky factors.
!
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use Cholesky, only: nSym, NumCho
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nIrrep, nBas_Aux(1:nIrrep), nBas(1:nIrrep)
logical(kind=iwp), intent(in) :: SubAux
integer(kind=iwp) :: i, iAdrQ, id, iOffQ1, iOpt, iost, ip_B, ip_B2, iSym, j, jSym, jVec, kSym, kVec, l_A, l_A_ht, l_A_t, l_B_t, &
                     l_Q, lRealName, Lu_Q, LUGAMMA, LuGamma2, LUAPT2, lVec, MaxMem, nBas2, nBasTri, nLR, nLRb(8), nseq, NumAux, &
                     NumCV, NumVecJ, NumVecK, nVec
real(kind=wp) :: aaa, Fac, TotCPU0, TotCPU1, TotWall0, TotWall1
logical(kind=iwp) :: is_error
character(len=4096) :: RealName
character(len=6) :: Name_Q
real(kind=wp), allocatable :: A(:), A_ht(:), A_t(:), B_t(:), QVec(:)
character(len=*), parameter :: SECNAM = 'Mult_with_Q_MP2'
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TotCPU0,TotWall0)
!                                                                      *
!***********************************************************************
!                                                                      *
do iSym=1,nSym
  NumCV = NumCho(iSym)
  NumAux = nBas_Aux(iSym)-1
  if (SubAux) NumAux = nBas_Aux(iSym)-1
  nLR = 0
  do jSym=1,nSym
    kSym = Mul(iSym,jSym)
    nLR = nLR+nBas(jSym)*nBas(kSym)
  end do
  nLRb(iSym) = nLR
end do
!                                                                      *
!***********************************************************************
!                                                                      *
do iSym=1,nSym

  nBas2 = nLRb(iSym)
  !write(u6,*) 'nBas2 = ',nBas2
  NumCV = NumCho(iSym)
  NumAux = nBas_Aux(iSym)-1
  if (SubAux) NumAux = nBas_Aux(iSym)-1

  ! Get Q-vectors from disk
  ! -----------------------

  l_Q = NumCV*NumAux
  call mma_allocate(QVec,l_Q,Label='Q_Vector')

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

  !LUCMOPT2 = 61
  !call PrgmTranslate('CMOPT2',RealName,lRealName)
  !call MOLCAS_Open_Ext2(LuCMOPT2,RealName(1:lRealName),'DIRECT','UNFORMATTED',iost,.false.,1,'OLD',is_error)

  !read(LuCMOPT2,iostat=iost) A_t
  !if (iost < 0) then
  !  write (6,*) 'Maybe, you did not add GRAD or GRDT keyword in &CASPT2?'
  !  write (6,*) 'Please add either one, if this is single-point gradient calculation.'
  !  write (6,*) 'Otherwise, something is wrong...'
  !  call abend()
  !end if

  ! Read A_PT2 from LUAPT2
  LUAPT2 = 77
  call daname_mf_wa(LUAPT2,'A_PT2')
  id = 0
  call ddafile(LUAPT2,2,A_t,l_A_t,id)

  ! Symmetrized A_PT2
  do i=1,NumCV
    do j=1,i
      aaa = Half*(A_t((j-1)*NumCV+i)+A_t((i-1)*NumCV+j))
      A_t((j-1)*NumCV+i) = aaa
      A_t((i-1)*NumCV+j) = aaa
    end do
  end do

# ifdef _DEBUGPRINT_
  write(u6,*) 'Q-vectors'
  do i=1,l_Q
    write(u6,*) QVec(i)
  end do

  write(u6,*) 'A-vectors'
  do i=1,l_A_t
    write(u6,*) A_t(i)
  end do
# endif

  ! Make first halftransformation to cholesky-base
  ! ----------------------------------------------

  call dGemm_('N','N',NumAux,NumCV,NumCV,One,QVec,NumAux,A_t,NumCV,Zero,A_ht,NumAux)

  call dGemm_('N','T',NumAux,NumAux,NumCV,One,A_ht,NumAux,QVec,NumAux,Zero,A,NumAux)

  ! Put transformed A-vectors back on disk

  !rewind(LuCMOPT2)
  !write(LuCMOPT2) A
  !close(LuCMOPT2)

  ! write A_PT2 to LUAPT2
  id = 0
  call ddafile(LUAPT2,1,A,l_A,id)

  call mma_deallocate(A_t)
  call mma_deallocate(A)
  call mma_deallocate(A_ht)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call mma_maxDBLE(MaxMem)
  MaxMem = 9*(MaxMem/10)
  call mma_allocate(B_t,MaxMem,Label='B_t')

  ! Applicable to C1 only
  nBasTri = nTri_Elem(nBas(1))
  nVec = MaxMem/(2*nBasTri+1) ! MaxMem/(2*nLRb(iSym)+1)
  nVec = min(max(NumCV,NumAux),nVec)
  if (nVec < 1) call SysAbendMsg(SecNam,'nVec is non-positive','[1]')

  l_B_t = nBasTri*nVec ! nLRb(iSym)*nVec
  ip_B = 1+l_B_t
  ip_B2 = ip_B+l_B_t

  LUGAMMA = 60
  call PrgmTranslate('GAMMA',RealName,lRealName)
  call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),'DIRECT','UNFORMATTED',iost,.true.,nBas2*8,'OLD',is_error)

  LuGamma2 = 62
  call PrgmTranslate('GAMMA2',RealName,lRealName)
  call MOLCAS_Open_Ext2(LuGamma2,RealName(1:lRealName),'DIRECT','UNFORMATTED',iost,.true.,NumAux*8,'REPLACE',is_error)

  ! The B-vectors should be read one batch at the time
  ! --------------------------------------------------

  do kVec=1,NumAux,nVec
    NumVecK = min(nVec,NumAux-(kVec-1))

    do jVec=1,NumCV,nVec
      NumVecJ = min(nVec,NumCV-(jVec-1))

      do lVec=1,NumVecJ
        read(Unit=LuGAMMA,Rec=jVec+lVec-1) B_t(ip_B2:ip_B2+nBas2-1)
        ! symmetrize (mu nu | P)
        ! only the lower triangle part is used
        nseq = 0
        do i=1,nBas(1)
          do j=1,i
            nseq = nseq+1
            aaa = B_t(ip_B2+i-1+nBas(1)*(j-1))+B_t(ip_B2+j-1+nBas(1)*(i-1))
            B_t(nseq+nBasTri*(lVec-1)) = aaa
          end do
        end do
      end do

      Fac = Zero
      if (jVec /= 1) Fac = One
      iOffQ1 = kVec+NumAux*(jVec-1)
      call dGemm_('N','T',nBasTri,NumVecK,NumVecJ,One,B_t,nBasTri,QVec(iOffQ1),NumAux,Fac,B_t(ip_B),nBasTri)
    end do

    ! Because scaling with 0.5 is omitted when symmetrized
    call DScal_(nBasTri*NumVecK,Half,B_t(ip_B),1)

    ! (mu nu | P) --> (P | mu nu)
    if (max(NumCV,NumAux) == nVec) then
      nseq = 0
      do lVec=1,nBas(1)
        do jVec=1,lVec
          nseq = nseq+1
          call DCopy_(NumAux,B_t(ip_B+nseq-1),nBasTri,B_t,1)
          write(unit=LuGAMMA2,rec=nseq) B_t(1:NumAux)
        end do
      end do
    else
      nseq = 0
      do lVec=1,nBas(1)
        do jVec=1,lVec
          nseq = nseq+1
          if (kVec /= 1) read(unit=LuGAMMA2,rec=nseq) B_t(1:kVec-1)
          call DCopy_(NumVecK,B_t(ip_B+nseq-1),nBasTri,B_t(kVec),1)
          write(unit=LuGAMMA2,rec=nseq) B_t(1:kVec+NumVecK-1)
        end do
      end do
    end if
  end do

  close(LuGAMMA2)

  call mma_deallocate(B_t)
  call mma_deallocate(QVec)

  call DaClos(Lu_Q)
  call DaClos(LUAPT2)
  close(LuGAMMA)

end do ! iSym

call CWTime(TotCPU1,TotWall1)
!write(u6,*) 'CPU/Wall Time for mult_with_q_caspt2:',totcpu1-totcpu0,totwall1-totwall0

return

end subroutine Mult_with_Q_CASPT2
