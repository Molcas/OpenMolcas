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

subroutine RHS_STRANS(NAS,NIS,ALPHA,lg_V1,lg_V2,ICASE,ISYM)
!SVC: this routine transforms RHS array V1 by multiplying on the left
!     with the S matrix and adds the result in V2: V2 <- V2 + alpha S*V1

use EQSOLV, only: IDSMAT
use fake_GA, only: GA_Arrays
use caspt2_global, only: LUSBT
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
use Constants, only: One
use Definitions, only: u6
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NAS, NIS, lg_V1, lg_V2, ICASE, ISYM
real(kind=wp), intent(in) :: ALPHA
integer(kind=iwp) :: IDS, NS
real(kind=wp), allocatable :: S(:)
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: iHiV1, iHiV2, iLoV1, iLoV2, jHiV1, jHiV2, jLoV1, jLoV2, LDV1, LDV2, lg_S, mV1, mV2, myRank, NCOL1, NCOL2, &
                     NROW1, NROW2
logical(kind=iwp) :: bStat
#include "global.fh"
#include "mafdecls.fh"

if (Is_Real_Par()) then
  if ((ICASE == 1) .or. (ICASE == 4)) then
    !-SVC: if case is A or C, the S/ST matrices are loaded as global arrays,
    !      then use the dgemm from GA to operate.
    call PSBMAT_GETMEM('S',lg_S,NAS)
    call PSBMAT_READ('S',iCase,iSym,lg_S,NAS)
    call GA_DGEMM('N','N',NAS,NIS,NAS,ALPHA,lg_S,lg_V1,One,lg_V2)
    bStat = GA_Destroy(lg_S)
  else
    !-SVC: if case is not A or C, the S/ST matrices are stored in replicate
    !      fashion, and the RHS are stored as vertical stripes, so use
    !      trimul on local memory, after accessing the local patch of the
    !      vector.
    NS = (NAS*(NAS+1))/2
    call mma_allocate(S,NS,Label='S')
    IDS = IDSMAT(ISYM,ICASE)
    call DDAFILE(LUSBT,2,S,NS,IDS)
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
      if ((NCOL1 /= NCOL2) .or. (NROW1 /= NROW2) .or. (NROW1 /= NAS)) then
        write(u6,*) 'RHS_STRANS: inconsistent stripe size'
        write(u6,'(A,I3)') 'ICASE = ',ICASE
        write(u6,'(A,I3)') 'ISYM  = ',ISYM
        write(u6,'(A,2I6)') 'NCOL1, NCOL2 = ',NCOL1,NCOL2
        write(u6,'(A,2I6)') 'NROW1, NROW2 = ',NROW1,NROW2
        call AbEnd()
      end if
      call GA_Access(lg_V1,iLoV1,iHiV1,jLoV1,jHiV1,mV1,LDV1)
      call GA_Access(lg_V2,iLoV2,iHiV2,jLoV2,jHiV2,mV2,LDV2)
      call TRIMUL(NAS,NCOL1,ALPHA,S,DBL_MB(mV1),LDV1,DBL_MB(mV2),LDV2)
      call GA_Release_Update(lg_V1,iLoV1,iHiV1,jLoV1,jHiV1)
      call GA_Release_Update(lg_V2,iLoV2,iHiV2,jLoV2,jHiV2)
    end if
    call GA_Sync()
    call mma_deallocate(S)
  end if
else
#endif
  NS = (NAS*(NAS+1))/2
  call mma_allocate(S,NS,Label='S')
  IDS = IDSMAT(ISYM,ICASE)
  call DDAFILE(LUSBT,2,S,NS,IDS)
  call TRIMUL(NAS,NIS,ALPHA,S,GA_Arrays(lg_V1)%A,NAS,GA_Arrays(lg_V2)%A,NAS)
  call mma_deallocate(S)
#ifdef _MOLCAS_MPP_
end if
#include "macros.fh"
unused_var(bStat)
#endif

end subroutine RHS_STRANS
