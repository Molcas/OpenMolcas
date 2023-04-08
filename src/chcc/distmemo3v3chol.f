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
        subroutine DistMemo3v3chol (NvGrp,maxdim,                       &
     &        PossV1,PossV2,PossV3,PossV4,                              &
     &        PossH1,PossH2,PossH3,PossH4,                              &
     &        PossM1,PossM2,PossM3,PossM4,PossM5,                       &
     &        PossK,PossQ,                                              &
     &        PossT)

!
!       This routine do:
!       define initial possitions of T,L,M and W arrays,
!       used in routine o3v3chol
!
!
!       I/O parameter description:
!       NvGrp    - # of groups in a,b,be,ga set (I)
!       maxdim   - maximal dimension of V'
!       Possx    - initial possitinos of arrays (O-all)
!       PossT    - initial and last possition (I/O)
!
!        requirements of o3v3chol step
!        K  - max{v'v'oo}
!        Q  - max{v'v'oo}
!        V1 - max{v'oo2, mv'v', v'v'oo, mv'o, moo, vv}
!        V2 - max{v'ooo, mv'v', v'v'oo, vo}
!        V3 - max{mv'o, v'v'oo, v'ooo}
!        V4 - max{v'ooo}
!        H1 - max{v'o}
!        H2 - max{v'o}
!        H3 - max{v'o}
!        H4 - max{v'o}
!        M1 - max{moo}
!        M2 - max{mv'o}
!        M3 - max{mv'o}
!        M4 - max{moo}
!        M5 - max{mv'o}
!
        implicit none
#include "chcc1.fh"
!
        integer NvGrp,maxdim
        integer PossV1,PossV2,PossV3,PossV4
        integer PossH1,PossH2,PossH3,PossH4
        integer PossM1,PossM2,PossM3,PossM4,PossM5
        integer PossK,PossQ
        integer PossT
!
!       help variables
        integer length
!
!
!1      Q,K (used also as X,Y)
!
        length=no*no*maxdim*maxdim
!
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
!
!2.1    V1 file max{v'oo2, mv'v', v'v'oo, mv'o, moo, vv}
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
!
        PossV1=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V1 ',PossV1,length
        end if
!
!2.2    V2 file - max{v'ooo, mv'v', v'v'oo, vo}
!
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
!
        PossV2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V2 ',PossV2,length
        end if
!
!2.2    V3 files - max{mv'o, v'v'oo, v'ooo}
!
        length=no*no*maxdim*maxdim
        if (no*no*no*maxdim.gt.length) then
          length=no*no*no*maxdim
        end if
        if (no*nc*maxdim.gt.length) then
          length=nc*no*maxdim
        end if
!
        PossV3=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V3 ',PossV3,length
        end if
!
!2.4    V4 files - max{v'ooo}
!
        length=no*no*no*maxdim
!
        PossV4=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V4 ',PossV4,length
        end if
!
!3      H1-4 files - max{v'o}
!
        length=no*maxdim
!
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
!
!4      M1-5 files
!        M1 - max{moo}
!        M2 - max{mv'o}
!        M3 - max{mv'o}
!        M4 - max{moo}
!        M5 - max{mv'o}
!
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
!
!
        if (printkey.ge.10) then
        write (6,*) 'PossT ',PossT
        end if
!
        return
! Avoid unused argument warnings
        if (.false.) Call Unused_integer(NvGrp)
        end
