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
! Copyright (C) 2011, Steven Vancoillie                                *
!***********************************************************************
!***********************************************************************
! Written by Steven Vancoillie, May 2011
! A set of subroutines that can handle RHS arrays in either a serial or
! parallel environment, depending on the situation.
!***********************************************************************
! --> when running serially, the RHS arrays are stored on LUSOLV and are
! loaded into the WORK array when needed.
! --> when running in parallel, the RHS arrays are stored on disk as
! disk resident arrays (DRAs) with filename RHS_XX_XX_XX, where XX is a
! number referring to the case, symmetry, and RHS vector respectively,
! and are loaded onto a global array when needed.
!***********************************************************************

subroutine RHS_SR2C(ITYP,IREV,NAS,NIS,NIN,lg_V1,lg_V2,ICASE,ISYM)
!SVC: this routine transforms the RHS arrays from SR format (V1) to C
!     format (V2) (IREV=0) and back (IREV=1), with ITYP specifying if
!     only the T matrix is used (ITYP=0) or the product of S and T
!     (ITYP=1).

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use caspt2_global, only: LUSBT
use EQSOLV, only: IDSTMAT, IDTMAT
use fake_GA, only: GA_Arrays
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ITYP, IREV, NAS, NIS, NIN, lg_V1, lg_V2, ICASE, ISYM
integer(kind=iwp) :: IDT
real(kind=wp), allocatable :: T(:)
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: iHiV1, iHiV2, iLoV1, iLoV2, jHiV1, jHiV2, jLoV1, jLoV2, LDV1, LDV2, lg_T, mV1, mV2, myRank, NCOL1, NCOL2, &
                     NROW1, NROW2
logical(kind=iwp) :: bStat
#include "global.fh"
#include "mafdecls.fh"

if (Is_Real_Par()) then
  if ((ICASE == 1) .or. (ICASE == 4)) then
    !-SVC: if case is A or C, the S/ST matrices are loaded as global arrays,
    !      then use the dgemm from GA to operate.
    call GA_CREATE_STRIPED('H',NAS,NIN,'TMAT',lg_T)
    if (ITYP == 0) then
      call PSBMAT_READ('T',iCase,iSym,lg_T,NAS*NIN)
    else if (ITYP == 1) then
      call PSBMAT_READ('M',iCase,iSym,lg_T,NAS*NIN)
    else
      write(u6,*) 'RHS_SR2C: invalid type = ',ITYP
      call AbEnd()
    end if

    if (IREV == 0) then
      call GA_DGEMM('N','N',NAS,NIS,NIN,One,lg_T,lg_V1,Zero,lg_V2)
    else
      call GA_DGEMM('T','N',NIN,NIS,NAS,One,lg_T,lg_V2,Zero,lg_V1)
    end if
    bStat = GA_Destroy(lg_T)
  else
    !-SVC: if case is not A or C, the S/ST matrices are stored in replicate
    !      fashion, and the RHS are stored as vertical stripes, so use dgemm
    !      on local memory, after accessing the local patch of the vector.
    call mma_allocate(T,NAS*NIN,Label='T')
    if (ITYP == 0) then
      IDT = IDTMAT(ISYM,ICASE)
    else if (ITYP == 1) then
      IDT = IDSTMAT(ISYM,ICASE)
    else
      write(u6,*) 'RHS_SR2C: invalid type = ',ITYP
      call AbEnd()
    end if
    call DDAFILE(LUSBT,2,T,NAS*NIN,IDT)
    !-SVC: get the local vertical stripes of the V1 and V2 vectors
    call GA_Sync()
    myRank = GA_NodeID()
    call GA_Distribution(lg_V1,myRank,iLoV1,iHiV1,jLoV1,jHiV1)
    call GA_Distribution(lg_V2,myRank,iLoV2,iHiV2,jLoV2,jHiV2)
    if ((jLoV1 /= 0) .and. (jLoV2 /= 0)) then
      NROW1 = iHiV1-iLoV1+1
      NROW2 = iHiV2-iLoV2+1
      NCOL1 = jHiV1-jLoV1+1
      NCOL2 = jHiV2-jLoV2+1
      if ((NCOL1 /= NCOL2) .or. (NROW1 /= NIN) .or. (NROW2 /= NAS)) then
        write(u6,*) 'RHS_SR2C: inconsistent stripe size'
        write(u6,'(A,I3)') 'ICASE = ',ICASE
        write(u6,'(A,I3)') 'ISYM  = ',ISYM
        write(u6,'(A,2I6)') 'NCOL1, NCOL2 = ',NCOL1,NCOL2
        write(u6,'(A,2I6)') 'NROW1, NIN   = ',NROW1,NIN
        write(u6,'(A,2I6)') 'NROW2, NAS   = ',NROW2,NAS
        call AbEnd()
      end if
      call GA_Access(lg_V1,iLoV1,iHiV1,jLoV1,jHiV1,mV1,LDV1)
      call GA_Access(lg_V2,iLoV2,iHiV2,jLoV2,jHiV2,mV2,LDV2)
      if (IREV == 0) then
        call DGEMM_('N','N',NAS,NCOL1,NIN,One,T,NAS,DBL_MB(mV1),LDV1,Zero,DBL_MB(mV2),LDV2)
      else
        call DGEMM_('T','N',NIN,NCOL1,NAS,One,T,NAS,DBL_MB(mV2),LDV2,Zero,DBL_MB(mV1),LDV1)
        !write(u6,*) 'Fingerprint =',RHS_DDOT(NAS,NIN,lg_V1,lg_V1)
      end if
      call GA_Release_Update(lg_V1,iLoV1,iHiV1,jLoV1,jHiV1)
      call GA_Release_Update(lg_V2,iLoV2,iHiV2,jLoV2,jHiV2)
    end if
    call mma_deallocate(T)
    call GA_Sync()
  end if
else
#endif
  call mma_allocate(T,NAS*NIN,Label='T')
  if (ITYP == 0) then
    IDT = IDTMAT(ISYM,ICASE)
  else if (ITYP == 1) then
    IDT = IDSTMAT(ISYM,ICASE)
  else
    write(u6,*) 'RHS_SR2C: invalid type = ',ITYP
    call AbEnd()
  end if
  call DDAFILE(LUSBT,2,T,NAS*NIN,IDT)
  if (IREV == 0) then
    call DGEMM_('N','N',NAS,NIS,NIN,One,T,NAS,GA_Arrays(lg_V1)%A,NIN,Zero,GA_Arrays(lg_V2)%A,NAS)
  else
    call DGEMM_('T','N',NIN,NIS,NAS,One,T,NAS,GA_Arrays(lg_V2)%A,NAS,Zero,GA_Arrays(lg_V1)%A,NIN)
  end if
  call mma_deallocate(T)
#ifdef _MOLCAS_MPP_
end if
#include "macros.fh"
unused_var(bStat)
#endif

end subroutine RHS_SR2C
