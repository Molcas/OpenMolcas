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

subroutine JoinLvec(wrk,wrksize,PosV1,PosV2,NaGrpR,LunAux)
! This routine does:
! Join L0-L2 files (Global, dimensioned as nc)
! from local L0-L2 files (dimensioned as ncLoc)
! N.B. This file have sense only for parallel run
!
! Structure of Cholesky vector files
!
! L0(m,IJ)    L0vctr  I>=J
!
! L1(m,I ,A') L1vcxx xx - Group of A'
!@@        kokot som, ze som to takto urobil, prerobit na L(m,a,i) to treba
!
! L2(m,A'B')  L2xxyy xx - Group of A', A'>=B'
!                    yy - Group of B'
!
! Memory requirements:
!  real:
! V1 - max {ov'm; v'v'm;  oom}
! V2 - max {v'v'm, ov'm; oom}
!  actual:
! reord routine requirements are used (DistMemReord)

#ifdef _MOLCAS_MPP_
use Index_Functions, only: nTri_Elem
use chcc_global, only: DimGrpaR, L0Name, L1Name, L2Name, nc, NChLoc, no
use Para_Info, only: MyRank
use Constants, only: Zero
#endif
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, PosV1, PosV2, NaGrpR, LunAux
real(kind=wp), intent(inout) :: wrk(wrksize)
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: aGrp, bGrp, dim_1, dima, dimab, dimb, ncLoc, ncOff
character(len=6) :: LunName

!7.0 calculate ncOffset for given node
ncLoc = NChLoc(myRank)
ncOff = sum(NChLoc(0:myRank-1))

!7.1 Make global L0
LunName = L0Name

!7.1.1 Read Local L0 = V2(ml,ij) on proper place in V2
dim_1 = ncLoc*nTri_Elem(no)
call GetX(wrk(PosV2),dim_1,LunAux,LunName,1,1)

!7.1.2 Vanish V1(m,ij)
dim_1 = nc*nTri_Elem(no)
wrk(PosV1:PosV1+dim_1-1) = Zero

!7.1.3 Insert V2(ml,ij) -> V1(m,ij)
dim_1 = nTri_Elem(no)
call InsL(wrk(PosV2),wrk(PosV1),ncLoc,nc,ncOff,dim_1)

!## Synchronizacny bod
!7.1.4 Allreduce V1
dim_1 = nc*nTri_Elem(no)
call gadgop(wrk(PosV1),dim_1,'+')

!7.1.5 Save L0 (Global)
dim_1 = nc*nTri_Elem(no)
call SaveX(wrk(PosV1),dim_1,LunAux,LunName,1,1)

!7.2 make global L1 files

do aGrp=1,NaGrpR
  dima = DimGrpaR(aGrp)
  LunName = L1Name(aGrp)

  !7.2.1 read L1 = V2(ml,i,a') back into file
  dim_1 = no*dima*ncLoc
  call GetX(wrk(PosV2),dim_1,LunAux,LunName,1,1)

  !7.2.2 Vanish V1(m,i,a')
  dim_1 = nc*no*dima
  wrk(PosV1:PosV1+dim_1-1) = Zero

  !7.2.3 Insert V2(ml,i,a') -> V1(m,i,a')
  dim_1 = no*dima
  call InsL(wrk(PosV2),wrk(PosV1),ncLoc,nc,ncOff,dim_1)

  !## Synchronizacny bod
  !7.2.4 Allreduce V1
  dim_1 = nc*no*dima
  call gadgop(wrk(PosV1),dim_1,'+')

  !7.2.5 Save V1(m,i,a') (Global)
  dim_1 = nc*no*dima
  call SaveX(wrk(PosV1),dim_1,LunAux,LunName,1,1)

end do

!2.3 Make glogal L2(m,a'b') files

do aGrp=1,NaGrpR
  dima = DimGrpaR(aGrp)
  do bGrp=1,aGrp
    dimb = DimGrpaR(bGrp)
    if (aGrp == bGrp) then
      dimab = nTri_Elem(dima)
    else
      dimab = dima*dimb
    end if
    LunName = L2Name(aGrp,bGrp)

    !7.3.1 read L2 = V2(ml,a'b') back into file
    dim_1 = dimab*ncLoc
    call GetX(wrk(PosV2),dim_1,LunAux,LunName,1,1)

    !7.3.2 Vanish V1(m,a'b')
    dim_1 = nc*dimab
    wrk(PosV1:PosV1+dim_1-1) = Zero

    !7.3.3 Insert V2(ml,a',b') -> V1(m,a',b')
    dim_1 = dimab
    call InsL(wrk(PosV2),wrk(PosV1),ncLoc,nc,ncOff,dim_1)

    !## Synchronizacny bod
    !7.3.4 Allreduce V1
    dim_1 = nc*dimab
    call gadgop(wrk(PosV1),dim_1,'+')

    !7.3.5 Save V1(ml,a',b') (Global)
    dim_1 = nc*dimab
    call SaveX(wrk(PosV1),dim_1,LunAux,LunName,1,1)

  end do
end do

#else
#include "macros.fh"
unused_var(wrk)
unused_var(PosV1)
unused_var(PosV2)
unused_var(NaGrpR)
unused_var(LunAux)
#endif

return

end subroutine JoinLvec
