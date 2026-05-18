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

subroutine MEMORY_ESTIMATE(JSYM,LBGRP,NBGRP,NCHOBUF,NPIQK,NADDBUF)

use Symmetry_Info, only: Mul
use definitions, only: iwp
use CHOVEC_IO, only: NVLOC_CHOBATCH
use caspt2_global, only: iParRHS, iPrGlb, iStpGrd
use PrintLevel, only: VERBOSE
use stdalloc, only: mma_MaxDBLE
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use caspt2_module, only: NSYM, NASHT, NISUP, NISH, NASH, NSSH, NTU, NTUV, NASH, NIGEJ, NIGTJ, NAGEB, NAGTB, NTGEU, NTGTU, NBTCHES, &
                         NBTCH

implicit none
integer(kind=iwp), intent(in) :: JSYM
integer(kind=iwp), intent(out) :: NBGRP
integer(kind=iwp), intent(out) :: LBGRP(2,*)
integer(kind=iwp), intent(out) :: NCHOBUF, NPIQK, NADDBUF
integer(kind=iwp) IB, IB1, IB2, IBGRP, ICASE, ISYI, ISYK, ISYP, ISYQ, MAXBUFF, MAXCHOL, MAXPIQK, MINBUFF, MINCHOL, MINGOOD, &
                  MINNICE, MINPIQK, MINSLOW, MXAVAIL, MXBATCH, MXCHOVEC, MXNPITOT, MXRHS, NCHOVEC, NCHUNK, NI, NK, NP, NPI, NQ, &
                  NQK, NV, NVECTOT
integer(kind=iwp), external :: iPARDIV
integer(kind=iwp), parameter :: Inactive = 1, Active = 2, Virtual = 3
integer(kind=iwp) nSh(8,3)
logical(kind=iwp) :: call_from_grad
integer(kind=iwp) :: ITYPE(4,9) = reshape([Inactive,Active ,Active  ,Active , &
                                           Inactive,Active ,Inactive,Active , &
                                           Inactive,Virtual,Active  ,Active , &
                                           Inactive,Virtual,Inactive,Virtual, &
                                           Active  ,Virtual,Active  ,Active , &
                                           Active  ,Virtual,Active  ,Virtual, &
                                           Active  ,Virtual,Inactive,Active , &
                                           Active  ,Virtual,Inactive,Virtual, &
                                           Inactive,Virtual,Inactive,Active ],[4,9])
integer(kind=iwp) ISYM

call_from_grad = .false.
if (iStpGrd == -1) call_from_grad = .true.

call ICopy(NSYM,NISH,1,nSh(1,Inactive),1)
call ICopy(NSYM,NASH,1,nSh(1,Active),1)
call ICopy(NSYM,NSSH,1,nSh(1,Virtual),1)

!SVC: compute the maximum RHS size that will be allocated per node, this
!     will be later allocated when addrhs routines are called. If global
!     arrays are used, this may actually be done outside of GETMEM, but
!     let's just reserve it anyway.
ISYM = JSYM
MXRHS = 2*NTU(ISYM)*NISUP(ISYM,5)
do ISYI=1,NSYM
  ISYM = ISYI
  ! case A
  MXRHS = max(MXRHS,NTUV(ISYM)*NISH(ISYM))
  ! case GP,GM
  MXRHS = max(MXRHS,NASH(ISYM)*NISUP(ISYM,10)+NASH(ISYM)*NISUP(ISYM,11))

  ISYM = Mul(JSYM,ISYI)
  ! case C
  MXRHS = max(MXRHS,NTUV(ISYM)*NSSH(ISYM))
  ! case EP,EM
  MXRHS = max(MXRHS,NASH(ISYM)*NISUP(ISYM,6)+NASH(ISYM)*NISUP(ISYM,7))
  do ISYK=1,NSYM
    ISYM = Mul(ISYI,ISYK)
    ! case BP,BM
    MXRHS = max(MXRHS,NTGEU(ISYM)*NIGEJ(ISYM))
    MXRHS = max(MXRHS,NTGTU(ISYM)*NIGTJ(ISYM))
    ! case HP,HM
    MXRHS = max(MXRHS,NAGEB(ISYM)*NIGEJ(ISYM))
    MXRHS = max(MXRHS,NAGTB(ISYM)*NIGTJ(ISYM))
    ! case F
    MXRHS = max(MXRHS,NTGEU(ISYM)*NAGEB(ISYM))
    MXRHS = max(MXRHS,NTGTU(ISYM)*NAGTB(ISYM))

    ISYM = Mul(ISYI,Mul(JSYM,ISYK))
    ! case D1,D2
    MXRHS = max(MXRHS,2*NTU(ISYM)*NISUP(ISYM,5))
  end do
end do
MXRHS = iPARDIV(MXRHS,0)

!SVC: determine maximum pair index size per symmetry and in total.
!     MXBFSZ=0
MXNPITOT = 0
do ISYI=1,NSYM
  ISYP = Mul(ISYI,JSYM)
  NI = max(NISH(ISYI),NASH(ISYI))
  NP = max(NASH(ISYP),NSSH(ISYP))
  MXNPITOT = MXNPITOT+NP*NI
  !MXBFSZ = MXBFSZ+NP*NI*NJSCT
end do

!SVC: determine maximum and minimum size to hold integral matrix.
!     cases here correspond to A, B, D1, H, C, F, D2, G, E.
!     with number icase =      1, 2,  3, 4, 5, 6,  7, 8, 9.
!     both case G and H can use blocking, currently this is done over
!     the number of inactive orbitals in one of the pair indices.
!     for case G, the minimum needed is NAU*NC, because of blocking NL.
!     for case H, the minimum needed is NA*NCL, because of blocking NJ.
!     to start, reserve space for TUVX integrals (NASHT**4)
MAXPIQK = NASHT**4
MINPIQK = NASHT**4
do ICASE=1,9
  do ISYI=1,NSYM
    NI = NSH(ISYI,ITYPE(1,ICASE))
    ISYP = Mul(ISYI,JSYM)
    NP = NSH(ISYP,ITYPE(2,ICASE))
    NPI = NP*NI
    do ISYK=1,NSYM
      NK = NSH(ISYK,ITYPE(3,ICASE))
      ISYQ = Mul(ISYK,JSYM)
      NQ = NSH(ISYQ,ITYPE(4,ICASE))
      NQK = NQ*NK
      MAXPIQK = max(MAXPIQK,NPI*NQK)
      if (ICASE == 4) then
        MINPIQK = max(MINPIQK,NP*NQK)
      else if (ICASE == 8) then
        MINPIQK = max(MINPIQK,NPI*NQ)
      else
        MINPIQK = max(MINPIQK,NPI*NQK)
      end if
    end do
  end do
end do
MAXBUFF = nint(sqrt(dble(MAXPIQK)))
MINBUFF = nint(sqrt(dble(MINPIQK)))

! In some cases, we do not use the buffer arrays
if (call_from_grad .or. (iParRHS == 2)) then
  MAXBUFF = 1
  MINBUFF = 1
end if

!SVC: total number of cholesky vectors
MXBATCH = 0
NVECTOT = 0
IB1 = NBTCHES(JSYM)+1
IB2 = NBTCHES(JSYM)+NBTCH(JSYM)
do IB=IB1,IB2
  NV = NVLOC_CHOBATCH(IB)
  MXBATCH = max(MXBATCH,NV)
  NVECTOT = NVECTOT+NV
end do
MAXCHOL = MXNPITOT*NVECTOT
MINCHOL = MXNPITOT*MXBATCH

!SVC: can we fit this all in memory?
call mma_MaxDBLE(MXAVAIL)

MINNICE = MXRHS+MAXPIQK+2*MAXBUFF+2*MAXCHOL
MINGOOD = MXRHS+MINPIQK+2*MINBUFF+2*MAXCHOL
MINSLOW = MXRHS+MINPIQK+2*MINBUFF+2*MINCHOL

if (call_from_grad) then
# ifdef _MOLCAS_MPP_
  if (is_real_par()) then
    !! One more vector is needed for buffer
    MAXCHOL = MAXCHOL+MXNPITOT
    MINNICE = 2*MXRHS+MAXPIQK+2*MAXBUFF+4*MAXCHOL
    MINGOOD = 2*MXRHS+MINPIQK+2*MINBUFF+4*MAXCHOL
    MINSLOW = 2*MXRHS+MINPIQK+2*MINBUFF+4*MINCHOL
  else
# endif
    MINNICE = MXRHS+MAXPIQK+2*MAXBUFF+4*MAXCHOL
    MINGOOD = MXRHS+MINPIQK+2*MINBUFF+4*MAXCHOL
    MINSLOW = MXRHS+MINPIQK+2*MINBUFF+4*MINCHOL
# ifdef _MOLCAS_MPP_
  end if
# endif
end if

if (IPRGLB > VERBOSE) then
  write(6,*)
  write(6,'(A,I1)') '  Memory estimates in RHSALL, SYM ',JSYM
  write(6,'(A,2X,I16)') '   allocatable:    ',MXAVAIL
  write(6,'(A,2X,I16)') '   recommended:    ',MINNICE
  write(6,'(A,2X,I16)') '   convenient:     ',MINGOOD
  write(6,'(A,2X,I16)') '   minimum:        ',MINSLOW
  write(6,*)
  if (MXAVAIL >= MINNICE) then
    write(6,*) ' I can use all cholesky vectors at once'
    write(6,*) ' as well as the whole integral matrixi.'
  else if (MXAVAIL >= MINGOOD) then
    write(6,*) ' I will group batches of cholesky vectors'
    write(6,*) ' and then maximize use of the integral matrix.'
  else if (MXAVAIL >= MINSLOW) then
    write(6,*) ' I will at least try to group batches.'
  else
    write(6,*) ' Do you see my problem?'
  end if
end if

if (MXAVAIL >= MINNICE) then
  ! group all batches, take the maximum needed for integrals and try
  ! to max out the buffer size, check it is larger than the minimum.
  NCHOBUF = MAXCHOL
  NBGRP = 1
  LBGRP(1,1) = IB1
  LBGRP(2,1) = IB2
  NADDBUF = MAXBUFF
  NPIQK = MAXPIQK
else if (MXAVAIL >= MINGOOD) then
  ! group all batches, take smaller buffer size, and try to max out
  ! integrals, and check they are lager than minimum needed
  NCHOBUF = MAXCHOL
  NBGRP = 1
  LBGRP(1,1) = IB1
  LBGRP(2,1) = IB2
  NADDBUF = MINBUFF
  NCHUNK = (MXAVAIL-MXRHS-2*NADDBUF-2*NCHOBUF)/MINPIQK
  if (call_from_grad) NCHUNK = (MXAVAIL-2*MXRHS-2*NADDBUF-4*NCHOBUF)/MINPIQK
  NPIQK = MINPIQK*NCHUNK
else if (MXAVAIL >= MINSLOW) then
  NADDBUF = MINBUFF
  NPIQK = MINPIQK
  NCHOBUF = (MXAVAIL-MXRHS-NPIQK-2*NADDBUF)/2
  if (call_from_grad) NCHOBUF = (MXAVAIL-2*MXRHS-NPIQK-4*NADDBUF)/2
  ! create batch groups that have at most MXCHOVEC cholesky vectors
  MXCHOVEC = max(NCHOBUF/MXNPITOT,1)
  NCHOVEC = 0
  IBGRP = 1
  LBGRP(1,IBGRP) = IB1
  do IB=IB1,IB2
    NV = NVLOC_CHOBATCH(IB)
    NCHOVEC = NCHOVEC+NV
    if (NCHOVEC > MXCHOVEC) then
      LBGRP(2,IBGRP) = IB-1
      IBGRP = IBGRP+1
      LBGRP(1,IBGRP) = IB
      NCHOVEC = NV
    end if
  end do
  LBGRP(2,IBGRP) = IB2
  NBGRP = IBGRP
else
  write(6,*)
  write(6,*) '  Not enough memory in RHSLL2...'
  write(6,'(A,I16)') '   allocatable:    ',MXAVAIL
  write(6,'(A,I16)') '   minimum:        ',MINSLOW
  call AbEnd()
end if

!SVC: sanity check, should not happen.
if (MXRHS > MXAVAIL-2*NCHOBUF-NPIQK-2*NADDBUF) then
  write(6,*)
  write(6,*) 'RHSALL2: RHS allocation will starve.'
  write(6,*) 'Possible bug in memory estimate.'
  write(6,*) 'This should not happen, please report.'
  write(6,*)
  write(6,'(2X,A8,2X,I14)') 'MXAVAIL ',MXAVAIL
  write(6,'(2X,A8,2X,I14)') 'MXRHS   ',MXRHS
  write(6,'(2X,A8,2X,I14)') 'NCHOBUF ',NCHOBUF
  write(6,'(2X,A8,2X,I14)') 'NPIQK   ',NPIQK
  write(6,'(2X,A8,2X,I14)') 'NADDBUF ',NADDBUF
  call AbEnd()
end if

if (IPRGLB > VERBOSE) then
  write(6,*)
  write(6,'(A16,2A16)') '  Buffer sizes:','           used','          ideal'
  write(6,'(A16,2I16)') '  ChoVecs:  ',2*NCHOBUF,2*MAXCHOL
  write(6,'(A16,2I16)') '  Integral: ',NPIQK,MAXPIQK
  write(6,'(A16,2I16)') '  Scatter:  ',2*NADDBUF,2*MAXBUFF
  write(6,*)
end if

end subroutine MEMORY_ESTIMATE
