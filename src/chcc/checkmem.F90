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

subroutine checkMem(NvGrp,NvSGrp,NchBlk,Jal1,Jal2,wrksize,maxdim)
! Check memory consumption according to the specified
! orbital segmentation

use chcc_global, only: PosFree
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: NvGrp, NvSGrp, NchBlk
integer(kind=iwp), intent(out) :: Jal1, Jal2, wrksize, maxdim
integer(kind=iwp) :: Jal3, mdGrpa, mdGrpbe, mdSGrpa, mdSGrpbe, PosH1, PosH2, PosH3, PosH4, PosH5, PosK, PosL11, PosL12, PosM1, &
                     PosM2, PosM3, PosM4, PosM5, PosMax, PosQ, PosT, PosV1, PosV2, PosV3, PosV4

!2 Distribute memory (and redefine wrksize)

!*.1 Distribute Permanent arrays

PosT = 1
call DistMemPerm(PosT)
PosMax = PosT

!*.2 Distribute Work arrays for Reord

call DefParReord(NvGrp,maxdim)
! also
call DefParo2v4(NvGrp,NvGrp,NvSGrp,NvSGrp,mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe)
PosT = PosFree
call DistMemReord(maxdim,mdSGrpa,NchBlk,PosV1,PosV2,PosV3,PosV4,PosM1,PosM2,PosT)
if (PosT > PosMax) PosMax = PosT

!*.3 Distribute Work arrays for o3v3

call DefParo3v3(NvGrp,maxdim)

PosT = PosFree
call DistMemo3v3jk(maxdim,PosV1,PosV2,PosV3,PosV4,PosH1,PosH2,PosH3,PosH4,PosH5,PosK,PosQ,PosT)
if (PosT > PosMax) PosMax = PosT

PosT = PosFree
call DistMemo3v3t2(maxdim,PosV1,PosV2,PosV3,PosV4,PosH1,PosH2,PosH3,PosH4,PosK,PosQ,PosT)
if (PosT > PosMax) PosMax = PosT

PosT = PosFree
call DistMemo3v3chol(maxdim,PosV1,PosV2,PosV3,PosV4,PosH1,PosH2,PosH3,PosH4,PosM1,PosM2,PosM3,PosM4,PosM5,PosK,PosQ,PosT)
if (PosT > PosMax) PosMax = PosT

!*.4 Distribute Work arrays for o2v4

call DefParo2v4(NvGrp,NvGrp,NvSGrp,NvSGrp,mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe)
PosT = PosFree
! Nazvy premennych tu davam ine, lebo ide len o vypocet narokov na Mem
call DistMemo2v4(NvGrp,NvGrp,NvSGrp,NvSGrp,mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe,PosV1,PosV2,PosV3,PosV4,PosL11,PosL12,PosH1,PosH2, &
                 PosH3,PosH4,PosQ,Jal1,Jal2,PosM1,PosM2,PosM3,PosM4,PosM5,Jal3,PosT,PosK)
if (PosT > PosMax) PosMax = PosT

!*.5 Memory requirements for Summary step are upperbounded by o3v3
!    and also o2v4, this it is skipped

wrksize = PosMax

return

end subroutine checkMem
