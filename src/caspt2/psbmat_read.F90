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
subroutine PSBMAT_READ(cNAME,iCase,iSym,lg_M,nSize)
!SVC20100902: read the disk array stored as cName+iSym using DRA
! interface into global array lg_M, or if replicate or serial, read from
! LUSBT into WORK(lg_M)

use Index_Functions, only: nTri_Elem
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
use GA_Wrapper, only: DBL_MB, GA_NodeId
use caspt2_global, only: LUH0T
#endif
use caspt2_global, only: LUSBT
use EQSOLV, only: IDBMAT, IDSMAT, IDSTMAT, IDTMAT
use fake_ga, only: GA_arrays
use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: cNAME
integer(kind=iwp), intent(in) :: iCASE, iSym, lg_M, nSize
integer(kind=iwp) :: IDISK, nBlock
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: IEND, iSTA, JEND, JSTA, LDM, LU, mpt_M, myRank
#endif

select case (CNAME)
  case ('S')
#   ifdef _MOLCAS_MPP_
    LU = LUH0T(1)
#   endif
    IDISK = IDSMAT(iSym,iCase)
    nBlock = nTri_Elem(nSize)
  case ('B')
#   ifdef _MOLCAS_MPP_
    LU = LUH0T(2)
#   endif
    IDISK = IDBMAT(iSym,iCase)
    nBlock = nTri_Elem(nSize)
  case ('T')
#   ifdef _MOLCAS_MPP_
    LU = LUH0T(3)
#   endif
    IDISK = IDTMAT(iSym,iCase)
    nBlock = nSize
  case ('M')
#   ifdef _MOLCAS_MPP_
    LU = LUH0T(4)
#   endif
    IDISK = IDSTMAT(iSym,iCase)
    nBlock = nSize
end select

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  call GA_Sync()
  myRank = GA_NodeID()
  call GA_Distribution(lg_M,myRank,ISTA,IEND,JSTA,JEND)
  if ((ISTA > 0) .and. (JSTA > 0)) then
    call GA_Access(lg_M,ISTA,IEND,JSTA,JEND,mpt_M,LDM)
    NBLOCK = LDM*(JEND-JSTA+1)
    call DDAFILE(LU,2,DBL_MB(mpt_M),NBLOCK,IDISK)
    call GA_Release_Update(lg_M,ISTA,IEND,JSTA,JEND)
  end if
  call GA_Sync()
else
#endif
  call DDAFILE(LUSBT,2,GA_Arrays(lg_M)%A(:),nBlock,IDISK)
#ifdef _MOLCAS_MPP_
end if
#endif

end subroutine PSBMAT_READ
