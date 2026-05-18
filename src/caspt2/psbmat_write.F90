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

! WRAPPER FOR PARALLEL S AND B MATRIX ROUTINES
subroutine PSBMAT_WRITE(cNAME,iCase,iSym,lg_M,nSize)
!SVC20100902: write the global array lg_M to disk using DRA interface,
! or if replicate or serial, write WORK(lg_M) to LUSBT

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
use caspt2_global, only: LUH0T
#endif
use caspt2_global, only: LUSBT
use EQSOLV, only: IDBMAT, IDSMAT, IDSTMAT, IDTMAT
use fake_ga, only: GA_arrays
use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: cNAME
integer(kind=iwp), intent(in) :: iCase, iSym, lg_M, nSize
integer(kind=iwp) :: IDISK, nBlock
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: IEND, ISTA, JEND, JSTA, LDM, LU, mpt_M, myRank
#include "global.fh"
#include "mafdecls.fh"
#endif

if (CNAME == 'S') then
# ifdef _MOLCAS_MPP_
  LU = LUH0T(1)
# endif
  IDISK = IDSMAT(iSym,iCase)
  nBlock = (nSize*(nSize+1))/2
else if (CNAME == 'B') then
# ifdef _MOLCAS_MPP_
  LU = LUH0T(2)
# endif
  IDISK = IDBMAT(iSym,iCase)
  nBlock = (nSize*(nSize+1))/2
else if (CNAME == 'T') then
# ifdef _MOLCAS_MPP_
  LU = LUH0T(3)
# endif
  IDISK = IDTMAT(iSym,iCase)
  nBlock = nSize
else if (CNAME == 'M') then
# ifdef _MOLCAS_MPP_
  LU = LUH0T(4)
# endif
  IDISK = IDSTMAT(iSym,iCase)
  nBlock = nSize
end if

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  call GA_Sync()
  myRank = GA_NodeID()
  call GA_Distribution(lg_M,myRank,ISTA,IEND,JSTA,JEND)
  if ((ISTA > 0) .and. (JSTA > 0)) then
    call GA_Access(lg_M,ISTA,IEND,JSTA,JEND,mpt_M,LDM)
    NBLOCK = LDM*(JEND-JSTA+1)
    call DDAFILE(LU,1,DBL_MB(mpt_M),NBLOCK,IDISK)
    call GA_Release(lg_M,ISTA,IEND,JSTA,JEND)
  end if
  call GA_Sync()
else
#endif
  call DDAFILE(LUSBT,1,GA_Arrays(lg_M)%A(:),nBlock,IDISK)
#ifdef _MOLCAS_MPP_
end if
#endif

end subroutine PSBMAT_WRITE
