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
        subroutine DistMemo3v3chol (NvGrp,maxdim,
     c        PossV1,PossV2,PossV3,PossV4,
     c        PossH1,PossH2,PossH3,PossH4,
     c        PossM1,PossM2,PossM3,PossM4,PossM5,
     c        PossK,PossQ,
     c        PossT)

c
c       This routine do:
c       define initial possitions of T,L,M and W arrays,
c       used in routine o3v3chol
c
c
c       I/O parameter description:
c       NvGrp    - # of groups in a,b,be,ga set (I)
c       maxdim   - maximal dimension of V'
c       Possx    - initial possitinos of arrays (O-all)
c       PossT    - initial and last possition (I/O)
c
c        requirements of o3v3chol step
c        K  - max{v'v'oo}
c        Q  - max{v'v'oo}
c        V1 - max{v'oo2, mv'v', v'v'oo, mv'o, moo, vv}
c        V2 - max{v'ooo, mv'v', v'v'oo, vo}
c        V3 - max{mv'o, v'v'oo, v'ooo}
c        V4 - max{v'ooo}
c        H1 - max{v'o}
c        H2 - max{v'o}
c        H3 - max{v'o}
c        H4 - max{v'o}
c        M1 - max{moo}
c        M2 - max{mv'o}
c        M3 - max{mv'o}
c        M4 - max{moo}
c        M5 - max{mv'o}
c
        implicit none
#include "chcc1.fh"
c
        integer NvGrp,maxdim
        integer PossV1,PossV2,PossV3,PossV4
        integer PossH1,PossH2,PossH3,PossH4
        integer PossM1,PossM2,PossM3,PossM4,PossM5
        integer PossK,PossQ
        integer PossT
c
c       help variables
        integer length
c
c
c1      Q,K (used also as X,Y)
c
        length=no*no*maxdim*maxdim
c
        PossQ=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM Q  ',PossQ,length
        end if
        PossK=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM K  ',PossK,length
        end if
c
c2.1    V1 file max{v'oo2, mv'v', v'v'oo, mv'o, moo, vv}
        length=maxdim*no*no*(no+1)/2
        if (maxdim*maxdim*nc.gt.length) then
          length=maxdim*maxdim*nc
        end if
        if (no*no*maxdim*maxdim.gt.length) then
          length=no*no*maxdim*maxdim
        end if
        if (no*nc*maxdim.gt.length) then
          length=no*nc*maxdim
        end if
        if (no*nc*no.gt.length) then
          length=no*no*nc
        end if
        if (nv*nv.gt.length) then
          length=nv*nv
        end if
c
        PossV1=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V1 ',PossV1,length
        end if
c
c2.2    V2 file - max{v'ooo, mv'v', v'v'oo, vo}
c
        length=no*no*maxdim*maxdim
        if (no*no*no*maxdim.gt.length) then
          length=no*no*no*maxdim
        end if
        if (nc*maxdim*maxdim.gt.length) then
          length=nc*maxdim*maxdim
        end if
        if (no*nv.gt.length) then
          length=no*nv
        end if
c
        PossV2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V2 ',PossV2,length
        end if
c
c2.2    V3 files - max{mv'o, v'v'oo, v'ooo}
c
        length=no*no*maxdim*maxdim
        if (no*no*no*maxdim.gt.length) then
          length=no*no*no*maxdim
        end if
        if (no*nc*maxdim.gt.length) then
          length=nc*no*maxdim
        end if
c
        PossV3=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V3 ',PossV3,length
        end if
c
c2.4    V4 files - max{v'ooo}
c
        length=no*no*no*maxdim
c
        PossV4=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V4 ',PossV4,length
        end if
c
c3      H1-4 files - max{v'o}
c
        length=no*maxdim
c
        PossH1=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H1 ',PossH1,length
        end if
        PossH2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H2 ',PossH2,length
        end if
        PossH3=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H3 ',PossH3,length
        end if
        PossH4=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H4 ',PossH4,length
        end if
c
c4      M1-5 files
c        M1 - max{moo}
c        M2 - max{mv'o}
c        M3 - max{mv'o}
c        M4 - max{moo}
c        M5 - max{mv'o}
c
        length=no*no*nc
        PossM1=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM M1 ',PossM1,length
        end if
        length=no*nc*maxdim
        PossM2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM M2 ',PossM2,length
        end if
        length=no*nc*maxdim
        PossM3=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM M3 ',PossM3,length
        end if
        length=no*nc*no
        PossM4=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM M4 ',PossM4,length
        end if
        length=no*nc*maxdim
        PossM5=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM M5 ',PossM5,length
        end if
c
c
        if (printkey.ge.10) then
        write (6,*) 'PossT ',PossT
        end if
c
        return
c Avoid unused argument warnings
        if (.false.) Call Unused_integer(NvGrp)
        end
