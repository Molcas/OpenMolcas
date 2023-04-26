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

subroutine DistMemReord(maxdim,maxdimSG,NchBlk,PosV1,PosV2,PosV3,PosV4,PosM1,PosM2,PosT)
! This routine does:
! define initial positions of OE,V1-V4,M1,2 arrays
! described in Reord routine
!
! Memory requirements:
!  intkey=0
! V1   - max {ov'ov'; nbas.nbas.m'; ov'm; v'v'm;  oom}
! V2   - max {ov'ov'; v'v'm, ov'm; oom}
! V3   - max {ov'm; oom}
! V4   - oom
!  intkey=0
! V1   - max {ov'ov'; nbas.nbas.m'; ov'm; v'v'm; oom; V"V"V"V"}
! V2   - max {ov'ov'; v'v'm, ov'm; oom}
! V3   - max {ov'm; oom; V'V'M}
! V4   - oom
! M1   - V"V"m
! M2   - max {V"V"M; OV"M)
!
! I/O parameter description:
! NxGrp    - # of groups in a,b,be,ga set (I)
! maxdim   - maximal dimension of V' (I)
! NChBlk   - # of Cholesky vectors in one Block - m' (I)
! Posx     - initial positions of arrays (O-all)
! PosT     - initial and last position (I/O)

use chcc_global, only: intkey, nc, no, nv, printkey
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: maxdim, maxdimSG, NchBlk
integer(kind=iwp), intent(out) :: PosV1, PosV2, PosV3, PosV4, PosM1, PosM2
integer(kind=iwp), intent(inout) :: PosT
integer(kind=iwp) :: length, nbas

!2 V1 file
!  V1 - max {ov'ov'; nbas.nbas.m'; ov'm; v'v'm;  oom, V"V"V"V"}

nbas = no+nv

length = maxdim*maxdim*no*no
if ((nbas*nbas*NChBlk) > length) length = nbas*nbas*NChBlk
if ((no*maxdim*nc) > length) length = no*maxdim*nc
if ((maxdim*maxdim*nc) > length) length = maxdim*maxdim*nc
if ((no*no*nc) > length) length = no*no*nc
if ((intkey == 1) .and. (length <= maxdimSG**4)) length = maxdimSG**4

PosV1 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM V1 ',PosV1,length

!3 V2 file
!  V2 - max {ov'ov'; v'v'm ; ov'm; oom}

length = maxdim*maxdim*no*no
if ((maxdim*maxdim*nc) > length) length = maxdim*maxdim*nc
if ((no*maxdim*nc) > length) length = no*maxdim*nc
if ((no*no*nc) > length) length = no*no*nc

PosV2 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM V2 ',PosV2,length

!4 V3 file
!  V3 - max {ov'm; oom, V'V'M}

if ((no*no*nc) > length) length = no*no*nc
if ((intkey == 1) .and. (length <= maxdim*maxdim*nc)) length = maxdim*maxdim*nc

PosV3 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM V3 ',PosV3,length

!5 V4 file
!  V4 - oom

length = no*no*nc

PosV4 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM V4 ',PosV4,length

!6 M1 - V"V"m

length = maxdimSG*maxdimSG*nc
if (intkey == 0) length = 0
PosM1 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM M1 ',PosM1,length

!7 M2 - max {V"V"M; OV"M)

length = maxdimSG*maxdimSG*nc
if (length < no*nc*maxdimSG) length = no*nc*maxdimSG
if (intkey == 0) length = 0
PosM2 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM M2 ',PosM2,length

return

end subroutine DistMemReord
