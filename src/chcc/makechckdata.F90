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

subroutine MakeChckData(wrk,wrksize,LunAux)
! this routine generates checkeroo data
! Use id possible only when NaGrp=NaSGRp=1
!
! assumption nc>=nv>=no (inac preverit dimenzovanie)
! dimensions: V1 - no*nv*nv2, nv2*nv2
!             V2 - nc*nv2
!             V3 - nc*nv*no

use Index_Functions, only: nTri_Elem
use chcc_global, only: I0Name, I1Name, I2Name, I3Name, L0Name, L1Name, L1k, L2Name, nc, no, nv, PosOE, Q21
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, LunAux
real(kind=wp), intent(inout) :: wrk(wrksize)
character(len=6) :: LunName
integer(kind=iwp) :: dim_1, PosT, PosV1, PosV2, PosV3

call DistMemChck(PosV1,PosV2,PosV3,PosT)

!1 make Q0
LunName = I0name
dim_1 = nTri_Elem(no)**2
call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)
call MkQ0(wrk(PosV1))

!2 make Q1
LunName = I1name(1)
dim_1 = nTri_Elem(no)*no*nv
call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)
call MkQ1(wrk(PosV1))

!3 Read Q21
LunName = I2name(1,1)
dim_1 = no*no*nv*nv
call GetX(Q21(1,1,1,1),dim_1,LunAux,LunName,1,1)

!4make Q22
LunName = I3name(1,1)
dim_1 = nTri_Elem(no)*nTri_Elem(nv)
call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)
call MkQ22(wrk(PosV1))

!5 make Q4, L2k
LunName = L2name(1,1)
dim_1 = nc*nTri_Elem(nv)
call GetX(wrk(PosV2),dim_1,LunAux,LunName,1,1)
call MkL2_chcc(wrk(PosV2))
dim_1 = nTri_Elem(nv)**2
wrk(PosV1:PosV1+dim_1-1) = Zero
dim_1 = nTri_Elem(nv)
call mc0c1at3b(nc,dim_1,nc,dim_1,dim_1,dim_1,dim_1,nc,dim_1,wrk(PosV2),wrk(PosV2),wrk(PosV1))
call MkQ4(wrk(PosV1))

! make Q3, L1k
LunName = L1name(1)
dim_1 = nc*no*nv
call GetX(wrk(PosV3),dim_1,LunAux,LunName,1,1)
call dcopy_(dim_1,wrk(PosV3),1,L1k,1)
dim_1 = no*nTri_Elem(nv)
wrk(PosV1:PosV1:dim_1-1) = Zero
dim_1 = nTri_Elem(nv)
call mc0c1at3b(nc,dim_1,nc,no*nv,dim_1,no*nv,dim_1,nc,no*nv,wrk(PosV2),wrk(PosV3),wrk(PosV1))
call MkQ3(wrk(PosV1))

! make L0
LunName = L0name
dim_1 = nc*nTri_Elem(no)
call GetX(wrk(PosV3),dim_1,LunAux,LunName,1,1)
call MkL0(wrk(PosV3))

! make OEo,OEv
call MkOE(wrk(PosOE))

call MkT1T2()

return

end subroutine MakeChckData
