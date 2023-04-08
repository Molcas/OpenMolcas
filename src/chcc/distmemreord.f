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
        subroutine DistMemReord (NaGrpR,maxdim,maxdimSG,NchBlk,         &
     &        PossV1,PossV2,PossV3,PossV4,PossM1,PossM2,                &
     &        PossT)

!
!       This routine do:
!       define initial possitions of OE,V1-V4,M1,2 arrays
!       described in Reord routine
!
!       Memory requirements:
!        intkey=0
!       V1   - max {ov'ov'; nbas.nbas.m'; ov'm; v'v'm;  oom}
!       V2   - max {ov'ov'; v'v'm, ov'm; oom}
!       V3   - max {ov'm; oom}
!       V4   - oom
!        intkey=0
!       V1   - max {ov'ov'; nbas.nbas.m'; ov'm; v'v'm; oom; V"V"V"V"}
!       V2   - max {ov'ov'; v'v'm, ov'm; oom}
!       V3   - max {ov'm; oom; V'V'M}
!       V4   - oom
!        M1   - V"V"m
!        M2   - max {V"V"M; OV"M)
!
!       I/O parameter description:
!       NxGrp    - # of groups in a,b,be,ga set (I)
!       maxdim   - maximal dimension of V' (I)
!       NChBlk   - # of Cholesky vectors in one Block - m' (I)
!       Possx    - initial possitinos of arrays (O-all)
!       PossT    - initial and last possition (I/O)
!
        implicit none
#include "chcc1.fh"
!
        integer NaGrpR,maxdim,maxdimSG,NchBlk
        integer PossV1,PossV2,PossV3,PossV4,PossM1,PossM2
        integer PossT
!
!       help variables
        integer length,nbas
!
!
!2      V1 file
!       V1   - max {ov'ov'; nbas.nbas.m'; ov'm; v'v'm;  oom, V"V"V"V"}
!
        nbas=no+nv
!
        length=maxdim*maxdim*no*no
        if ((nbas*nbas*NChBlk).gt.length) then
          length=nbas*nbas*NChBlk
        end if
        if ((no*maxdim*nc).gt.length) then
          length=no*maxdim*nc
        end if
        if ((maxdim*maxdim*nc).gt.length) then
          length=maxdim*maxdim*nc
        end if
        if ((no*no*nc).gt.length) then
          length=no*no*nc
        end if
        if ((intkey.eq.1).and.(length.le.maxdimSG**4)) then
          length=maxdimSG**4
        end if
!
        PossV1=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V1 ',PossV1,length
        end if
!
!3      V2 file
!       V2   - max {ov'ov'; v'v'm ; ov'm; oom}
!
        length=maxdim*maxdim*no*no
        if ((maxdim*maxdim*nc).gt.length) then
          length=maxdim*maxdim*nc
        end if
        if ((no*maxdim*nc).gt.length) then
          length=no*maxdim*nc
        end if
        if ((no*no*nc).gt.length) then
          length=no*no*nc
        end if
!
        PossV2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V2 ',PossV2,length
        end if
!
!
!4      V3 file
!       V3   - max {ov'm; oom, V'V'M}
!
        length=no*maxdim*nc
        if ((no*no*nc).gt.length) then
          length=no*no*nc
        end if
        if ((intkey.eq.1).and.(length.le.maxdim*maxdim*nc)) then
          length=maxdim*maxdim*nc
        end if
!
        PossV3=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V3 ',PossV3,length
        end if
!
!5      V4 file
!       V4   - oom
!
        length=no*no*nc
!
        PossV4=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V4 ',PossV4,length
        end if
!
!6        M1   - V"V"m
!
        length=maxdimSG*maxdimSG*nc
        if (intkey.eq.0) then
          length=0
        end if
        PossM1=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM M1 ',PossM1,length
        end if
!
!7        M2   - max {V"V"M; OV"M)
!
        length=maxdimSG*maxdimSG*nc
        if (length.lt.no*nc*maxdimSG) then
          length=no*nc*maxdimSG
        end if
        if (intkey.eq.0) then
          length=0
        end if
        PossM2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM M2 ',PossM2,length
        end if
!
        return
! Avoid unused argument warnings
        if (.false.) Call Unused_integer(NaGrpR)
        end
