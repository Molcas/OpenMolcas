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
        subroutine DistMemo3v3jk (NvGrp,maxdim,
     c        PossV1,PossV2,PossV3,PossV4,
     c        PossH1,PossH2,PossH3,PossH4,PossH5,
     c        PossK,PossQ,
     c        PossT)

c
c       This routine do:
c       define initial possitions of H,V,XY
c       described in o3v3jk routine
c
c
c       I/O parameter description:
c       NvGrp    - # of groups in a,b,be,ga set (I)
c       maxdim   - maximal dimension of V'
c       Possx    - initial possitinos of arrays (O-all)
c       PossT    - initial and last possition (I/O)
c
c        requirements for o3v3jk:
c        H1 - max {v'o}
c        H2 - max {v'o}
c        H3 - max {v'v'}
c        H4 - max {v'o}
c        H5 - max {v'o}
c        V1 - max {v'ooo, v'v'oo, o2oo}
c        V2 - max {oooo, v'v'oo, v'ooo}
c        V3 - max {v'v'oo}
c        V4 - max {v'ooo}
c        PX - max {v'v'oo}
c        QY - max {v'v'oo}
c
        implicit none
#include "chcc1.fh"
c
        integer NvGrp,maxdim
        integer PossV1,PossV2,PossV3,PossV4
        integer PossH1,PossH2,PossH3,PossH4,PossH5
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
c2.1    V1 file - max {v'ooo, v'v'oo, o2oo}
c
        length=no*no*maxdim*maxdim
        if (no*maxdim*nc.gt.length) then
          length=no*maxdim*nc
        end if
        if (no*no*no*maxdim.gt.length) then
          length=no*no*no*maxdim
        end if
        if (no*no*no*(no+1)/2.gt.length) then
          length=no*no*no*(no+1)/2
        end if
c
        PossV1=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V1 ',PossV1,length
        end if
c
c2.2    V2 files - max {oooo, v'v'oo, v'ooo}
c
        length=no*no*maxdim*maxdim
        if (no*no*no*maxdim.gt.length) then
          length=no*no*no*maxdim
        end if
        if (no*no*no*no.gt.length) then
          length=no*no*no*no
        end if
c
        PossV2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V2 ',PossV2,length
        end if
c
c2.3    V3 file - max {v'v'oo}
c
        length=no*no*maxdim*maxdim
        PossV3=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V3 ',PossV3,length
        end if
c
c2.4    V4 file - max {v'ooo}
c
        length=no*no*no*maxdim
        PossV4=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V4 ',PossV4,length
        end if
c
c3      H1,2 files
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
c
c3.2    H3 file
c
        length=maxdim*maxdim
        if (no*maxdim.gt.length) then
          length=no*maxdim
        end if
c
        PossH3=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H3 ',PossH3,length
        end if
c
c3.3        H4,H5 file
c
        length=no*maxdim
c
        PossH4=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H4 ',PossH4,length
        end if
        PossH5=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H5 ',PossH5,length
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
