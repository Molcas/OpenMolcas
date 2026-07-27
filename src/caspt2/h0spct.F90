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

subroutine H0SPCT()
! Write pertinent warnings and statistics for the energy
! denominators, i.e. the spectrum of (H0(diag)-E0).

use PrintLevel, only: VERBOSE
use SC_NEVPT2, only: SC_prop
#ifdef _MOLCAS_MPP_
use allgather_wrapper, only: allgather
use Para_Info, only: Is_Real_Par
use GA_Wrapper, only: DBL_MB, GA_NodeId
#endif
use EQSOLV, only: IDBMAT, IRHS, IVECX
use fake_GA, only: GA_Arrays
use caspt2_global, only: cmpThr, cntThr, dnmThr, iPrGlb, LUSBT
use caspt2_module, only: NSYM, NASUP, NISUP, NINDEP, CASES, ORBNAM
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: IAEND, IAS, IASTA, IBUF, ICASE, IIEND, IIS, IISTA, IP, IQ, IR, IS, ISYM, JD, lg_RHS, lg_VEC, MAXBUF, NAS, &
                     NBUF, NIN, NIS
real(kind=wp) :: COEF, DNOM, ECNT, RHS
character(len=80) :: LINE
integer(kind=iwp), allocatable, target :: IDXBUF(:,:)
integer(kind=iwp), pointer :: IDX(:,:)
real(kind=wp), allocatable :: BD(:), ID(:)
real(kind=wp), allocatable, target :: VALBUF(:,:)
real(kind=wp), pointer :: VAL(:,:)
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: LD, myRank, mRHS, mVEC
integer(kind=iwp), allocatable, target :: IDX_H(:,:)
real(kind=wp), allocatable, target :: VAL_H(:,:)
#endif

! While it is possible to do a similar wf inspection,
! the denominators etc. are very differently formulated for SC-NEVPT2, so we just skip here.
if (SC_prop) then
  write(u6,*) '  (Skip denominator check etc. with SC-NEVPT2 wavefunction)'
  return
end if

write(u6,*)
call CollapseOutput(1,'Denominators, etc.')
write(u6,'(A)') repeat('-',110)
write(u6,'(A)') ' Report on small energy denominators, large coefficients, and large energy contributions.'

if (IPRGLB >= VERBOSE) then
  write(u6,'(A)') '  The ACTIVE-MIX index denotes linear combinations which gives ON expansion functions'
  write(u6,'(A)') '  and makes H0 diagonal within type.'
  write(u6,'(A)') '  DENOMINATOR: The (H0_ii - E0) value from the above-mentioned diagonal approximation.'
  write(u6,'(A)') '  RHS VALUE  : Right-Hand Side of CASPT2 Eqs.'
  write(u6,'(A)') '  COEFFICIENT: Multiplies each of the above ON terms in the first-order wave function.'
  write(u6,'(A)') ' Thresholds used:'
  write(u6,'(a,f7.4)') '         Denominators:',DNMTHR
  write(u6,'(a,f7.4)') '         Coefficients:',CMPTHR
  write(u6,'(a,f7.4)') ' Energy contributions:',CNTTHR
  write(u6,*)
end if

write(u6,'(A)') 'CASE  SYMM ACTIVE-MIX  NON-ACTIVE INDICES          DENOMINATOR     RHS VALUE       COEFFICIENT     CONTRIBUTION'

!SVC: initial buffer size, will be reallocated on the fly
MAXBUF = 1024
call mma_allocate(IDXBUF,2,MAXBUF,LABEL='IDXBUF')
call mma_allocate(VALBUF,4,MAXBUF,LABEL='VALBUF')

! Very long loop over symmetry and case:
do ICASE=1,13
  do ISYM=1,NSYM
    NAS = NASUP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    if (NIS == 0) cycle
    NIN = NINDEP(ISYM,ICASE)
    if (NIN == 0) cycle
    LINE(1:12) = CASES(ICASE)//'    '
    write(LINE(10:10),'(i1)') ISYM

    ! Remember: NIN values in BDIAG, but must read NAS for correct
    ! positioning.
    call mma_allocate(BD,NAS,LABEL='BD')
    call mma_allocate(ID,NIS,LABEL='ID')
    JD = IDBMAT(ISYM,ICASE)
    call DDAFILE(LUSBT,2,BD,NAS,JD)
    call DDAFILE(LUSBT,2,ID,NIS,JD)

    call RHS_ALLO(NIN,NIS,lg_RHS)
    call RHS_ALLO(NIN,NIS,lg_VEC)
    call RHS_READ_SR(lg_RHS,ICASE,ISYM,IRHS)
    call RHS_READ_SR(lg_VEC,ICASE,ISYM,IVECX)
    IBUF = 0
#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      ! Get the superindex ranges of this process's block. If no elements are
      ! owned by a process, then ilo=0 and ihi=-1 such that the loop further
      ! down will just be skipped.
      call GA_Sync()
      myRank = GA_NodeID()
      call GA_Distribution(lg_RHS,myRank,IASTA,IAEND,IISTA,IIEND)
      if ((IASTA /= 0) .and. (IAEND-IASTA+1 /= NIN)) then
        write(u6,*) 'RHSOD: mismatch in range of the superindices'
        call AbEnd()
      end if
      ! if the block is non-empty, loop over its elements
      if ((IASTA > 0) .and. (IISTA > 0)) then
        call GA_Access(lg_RHS,IASTA,IAEND,IISTA,IIEND,mRHS,LD)
        call GA_Access(lg_VEC,IASTA,IAEND,IISTA,IIEND,mVEC,LD)
        if (LD /= NIN) then
          write(u6,*) 'RHSOD: assumption NAS=LDW wrong, abort'
          call AbEnd()
        end if
      end if
    else
#   endif
      IASTA = 1
      IAEND = NIN
      IISTA = 1
      IIEND = NIS
#   ifdef _MOLCAS_MPP_
    end if
#   endif

    !*******************************************************************
    ! inner loop over RHS elements in symmetry ISYM
    !*******************************************************************
    do IIS=IISTA,IIEND
      do IAS=IASTA,IAEND
        DNOM = BD(IAS)+ID(IIS)
#       ifdef _MOLCAS_MPP_
        if (Is_Real_Par()) then
          RHS = DBL_MB(mRHS+IAS-1+NIN*(IIS-IISTA))
          COEF = DBL_MB(mVEC+IAS-1+NIN*(IIS-IISTA))
        else
#       endif
          RHS = GA_Arrays(lg_RHS)%A(IAS+NIN*(IIS-IISTA))
          COEF = GA_Arrays(lg_VEC)%A(IAS+NIN*(IIS-IISTA))
#       ifdef _MOLCAS_MPP_
        end if
#       endif
        ECNT = COEF*RHS
        if ((abs(DNOM) < DNMTHR) .or. (abs(COEF) > CMPTHR) .or. (abs(ECNT) > CNTTHR)) then
          if (IBUF < MAXBUF) then
            IBUF = IBUF+1
            IDXBUF(1,IBUF) = IAS
            IDXBUF(2,IBUF) = IIS
            VALBUF(1,IBUF) = DNOM
            VALBUF(2,IBUF) = RHS
            VALBUF(3,IBUF) = COEF
            VALBUF(4,IBUF) = ECNT
          end if
        end if
      end do
    end do

#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      if ((IASTA > 0) .and. (IISTA > 0)) then
        call GA_Release(lg_RHS,IASTA,IAEND,IISTA,IIEND)
        call GA_Release(lg_VEC,IASTA,IAEND,IISTA,IIEND)
      end if
    end if
#   endif
    NBUF = IBUF
#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      call GAIGOP_SCAL(NBUF,'+')
      call mma_allocate(IDX_H,2,NBUF,LABEL='IDX_H')
      call mma_allocate(VAL_H,4,NBUF,LABEL='VAL_H')
      call allgather(IDXBUF,2*IBUF,IDX_H,2*NBUF)
      call allgather(VALBUF,4*IBUF,VAL_H,4*NBUF)
      IDX => IDX_H
      VAL => VAL_H
    else
#   endif
      IDX => IDXBUF
      VAL => VALBUF
#   ifdef _MOLCAS_MPP_
    end if
#   endif

    do IBUF=1,NBUF
      IAS = IDX(1,IBUF)
      IIS = IDX(2,IBUF)
      DNOM = VAL(1,IBUF)
      RHS = VAL(2,IBUF)
      COEF = VAL(3,IBUF)
      ECNT = VAL(4,IBUF)
      if ((ICASE == 12) .or. (ICASE == 13)) then
        call EXCIND(IAS,IIS,ISYM,ICASE,IP,IQ,IR,IS)
        LINE(13:20) = ORBNAM(IP)
        LINE(21:28) = ORBNAM(IQ)
        LINE(29:36) = ORBNAM(IR)
        LINE(37:44) = ORBNAM(IS)
        LINE(45:46) = '  '
      else
        write(LINE(13:22),'(A2,I1,A1,I4.4)') 'Mu',ISYM,'.',IAS
        call NSIND(IIS,ISYM,ICASE,IP,IQ,IR)
        LINE(23:30) = ORBNAM(IP)
        LINE(31:46) = '                '
        if (IQ > 0) LINE(31:38) = ORBNAM(IQ)
        if (IR > 0) LINE(39:46) = ORBNAM(IR)
      end if
      write(u6,'(A,4F16.8)') LINE(1:46),DNOM,RHS,COEF,ECNT
    end do

#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      call mma_deallocate(IDX_H)
      call mma_deallocate(VAL_H)
    end if
#   endif

    call RHS_FREE(lg_RHS)
    call RHS_FREE(lg_VEC)

    call mma_deallocate(BD)
    call mma_deallocate(ID)

    ! End of very long loop over symmetry and case:
  end do
end do

call mma_deallocate(IDXBUF)
call mma_deallocate(VALBUF)
nullify(IDX,VAL)

call CollapseOutput(0,'Denominators, etc.')

end subroutine H0SPCT
