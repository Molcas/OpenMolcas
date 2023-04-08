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
        subroutine DistMemo3v3t2 (NvGrp,maxdim,                         &
     &        PossV1,PossV2,PossV3,PossV4,                              &
     &        PossH1,PossH2,PossH3,PossH4,                              &
     &        PossK,PossQ,                                              &
     &        PossT)

!
!       This routine do:
!       define initial possitions of H,V, QK
!
!
!       I/O parameter description:
!       NvGrp    - # of groups in a,b,be,ga set (I)
!       maxdim   - maximal dimension of V'
!       Possx    - initial possitinos of arrays (O-all)
!       PossT    - initial and last possition (I/O)
!
!        requirements for o3v3t2:
!        H1 - max {v'o}
!        H2 - max {v'o}
!        H3 - max {v'v',ooo}
!        H4 - max {v'o}
!        V1 - max {v'ov'o, vv'}
!        V2 - max {v'ov'o}
!        V3 - max {v'ov'o, o2oo}
!        V4 - max {v'ov'o}
!        PX - max {v'ov'o}
!        QY - max {v'ov'o}
!
        implicit none
#include "chcc1.fh"
!
        integer NvGrp,maxdim
        integer PossV1,PossV2,PossV3,PossV4
        integer PossH1,PossH2,PossH3,PossH4
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
!2.1    V1 file - max {v'ov'o, vv'}
!
        length=no*no*maxdim*maxdim
        if (nv*maxdim.gt.length) then
          length=maxdim*nv
        end if
!
        PossV1=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V1 ',PossV1,length
        end if
!
!2.2    V2-4 files
!        V2 - max {v'ov'o}
!        V3 - max {v'ov'o, o2oo}
!        V4 - max {v'ov'o}
!
        length=no*no*maxdim*maxdim
        if (no*no*no*maxdim.gt.length) then
          length=no*no*no*maxdim
        end if
!
        PossV2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V2 ',PossV2,length
        end if
!
        length=no*no*maxdim*maxdim
        if (no*(no+1)*no*no/2.gt.length) then
          length=no*no*no*(no+1)/2
        end if
!
        PossV3=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V3 ',PossV3,length
        end if
!
        length=no*no*maxdim*maxdim
        PossV4=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V4 ',PossV4,length
        end if
!
!3      H1,2 files
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
!
!3.2    H3 file
!
        length=maxdim*maxdim
        if (no*maxdim.gt.length) then
          length=no*maxdim
        end if
        if (no*no*no.gt.length) then
          length=no*no*no
        end if
!
        PossH3=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H3 ',PossH3,length
        end if
!
!3.3        H4 file
!
        length=no*maxdim
!
        PossH4=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H4 ',PossH4,length
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
