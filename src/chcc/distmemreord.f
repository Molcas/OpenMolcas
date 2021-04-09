************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
        subroutine DistMemReord (NaGrpR,maxdim,maxdimSG,NchBlk,
     c        PossV1,PossV2,PossV3,PossV4,PossM1,PossM2,
     c        PossT)

c
c       This routine do:
c       define initial possitions of OE,V1-V4,M1,2 arrays
c       described in Reord routine
c
c       Memory requirements:
c        intkey=0
c       V1   - max {ov'ov'; nbas.nbas.m'; ov'm; v'v'm;  oom}
c       V2   - max {ov'ov'; v'v'm, ov'm; oom}
c       V3   - max {ov'm; oom}
c       V4   - oom
c        intkey=0
c       V1   - max {ov'ov'; nbas.nbas.m'; ov'm; v'v'm; oom; V"V"V"V"}
c       V2   - max {ov'ov'; v'v'm, ov'm; oom}
c       V3   - max {ov'm; oom; V'V'M}
c       V4   - oom
c        M1   - V"V"m
c        M2   - max {V"V"M; OV"M)
c
c       I/O parameter description:
c       NxGrp    - # of groups in a,b,be,ga set (I)
c       maxdim   - maximal dimension of V' (I)
c       NChBlk   - # of Cholesky vectors in one Block - m' (I)
c       Possx    - initial possitinos of arrays (O-all)
c       PossT    - initial and last possition (I/O)
c
        implicit none
#include "chcc1.fh"
c
        integer NaGrpR,maxdim,maxdimSG,NchBlk
        integer PossV1,PossV2,PossV3,PossV4,PossM1,PossM2
        integer PossT
c
c       help variables
        integer length,nbas
c
c
c2      V1 file
c       V1   - max {ov'ov'; nbas.nbas.m'; ov'm; v'v'm;  oom, V"V"V"V"}
c
        nbas=no+nv
c
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
c
        PossV1=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V1 ',PossV1,length
        end if
c
c3      V2 file
c       V2   - max {ov'ov'; v'v'm ; ov'm; oom}
c
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
c
        PossV2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V2 ',PossV2,length
        end if
c
c
c4      V3 file
c       V3   - max {ov'm; oom, V'V'M}
c
        length=no*maxdim*nc
        if ((no*no*nc).gt.length) then
          length=no*no*nc
        end if
        if ((intkey.eq.1).and.(length.le.maxdim*maxdim*nc)) then
          length=maxdim*maxdim*nc
        end if
c
        PossV3=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V3 ',PossV3,length
        end if
c
c5      V4 file
c       V4   - oom
c
        length=no*no*nc
c
        PossV4=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V4 ',PossV4,length
        end if
c
c6        M1   - V"V"m
c
        length=maxdimSG*maxdimSG*nc
        if (intkey.eq.0) then
          length=0
        end if
        PossM1=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM M1 ',PossM1,length
        end if
c
c7        M2   - max {V"V"M; OV"M)
c
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
c
        return
c Avoid unused argument warnings
        if (.false.) Call Unused_integer(NaGrpR)
        end
