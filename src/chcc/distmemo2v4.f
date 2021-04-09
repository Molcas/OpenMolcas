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
        subroutine DistMemo2v4 (NaGrp,NbeGrp,NaSGrp,NbeSgrp,
     c                 mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe,
     c                 PossTau,PossT2n1,PossT2n2,PossT2w,
     c                 PossL11,PossL12,
     c                 PossL21,PossL22,PossL23,PossL24,PossL2W,
     c                 PossH1,PossH2,
     c                 PossM1,PossM2,PossW1,PossW2,PossWw,PossWx,
     c                 PossT,NL2)
c
c       This routine do:
c       define initial possitions of T,L,M and W arrays,
c       described in o2v4ctl routine
c
c
c       I/O parameter description:
c       NaGrp    - # of groups in a set (I)
c       NbeGrp   - # of groups in be set (I)
c       NaSGrp   - # of subgroups in each (a)' group (I)
c       NbeSGrp  - # of subgroups in each (be)' group (I)
c       mdGrpa   - # maximal dimension od (a)' group (I)
c       mdGrpbe  - # maximal dimension od (be)' group (I)
c       mdSGrpa  - # maximal dimension od (a)" subgroup (I)
c       mdSGrpbe - # maximal dimension od (be)" subgroup (I)
c       PossX    - initial possitinos of arrays (O-all)
c       PossT    - next free possition (O)
c       NL2      - # of L2 vectors (O)
c
        implicit none
#include "chcc1.fh"
c
        integer NaGrp,NbeGrp,NaSGrp,NbeSgrp
        integer mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe
        integer PossTau,PossT2n1,PossT2n2,PossT2w
        integer PossL11,PossL12
        integer PossH1,PossH2
        integer PossL21,PossL22,PossL23,PossL24,PossL2W
        integer PossM1,PossM2,PossW1,PossW2,PossWw,PossWx
        integer PossT,NL2
c
c       help variables
        integer length,lenab,lenbega
c
c1      Tau
c
        PossTau=possT
        if (NaGrp.eq.1) then
          length=no*no*nv*(nv+1)/2
        else
          length=no*no*mdGrpa*mdGrpa
        end if
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,99) 'DM Tau',PossTau
        end if
c
c2      T2n files
c
        possT2n1=PossT
        if ((NbeGrp.eq.1).and.(NbeSGrp.eq.1)) then
          length=no*(no+1)*nv*(nv+1)/4
        else
          length=no*(no+1)*mdSGrpbe*mdSGrpbe/2
        end if
        PossT=PossT+length
c
        possT2n2=PossT
        if ((NbeGrp.eq.1).and.(NbeSGrp.eq.1)) then
          length=no*(no-1)*nv*(nv-1)/4
        else
          length=no*(no-1)*mdSGrpbe*mdSGrpbe/2
        end if
        PossT=PossT+length
c
        possT2w=PossT
        if ((NbeGrp.eq.1).and.(NbeSGrp.eq.1)) then
          length=no*(no+1)*nv*(nv+1)/4
        else
          length=no*(no+1)*mdSGrpbe*mdSGrpbe/2
        end if
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,99) 'DM T2 ',PossT2n1,PossT2n2,PossT2w
        end if
c
c       L1 files
c
        if (intkey.eq.1) then
c        integral based
          length=0
        else
c        cholesky based
          length=nc*mdGrpa*no
        end if
c
        if (NaGrp.eq.1) then
          PossL11=PossT
          PossL12=PossT
          PossT=PossT+length
        else
          PossL11=PossT
          PossT=PossT+length
          PossL12=PossT
          PossT=PossT+length
        end if
        if (printkey.ge.10) then
        write (6,99) 'DM L1 ',PossL11,PossL12
        end if
c
c       L2 files
c
        if (intkey.eq.1) then
c        integral based
          length=0
        else
c        cholesky based
          length=nc*mdGrpa*mdGrpbe
        end if
c
        if (NaGrp.eq.1) then
          if (NbeGrp.eq.1) then
c         All L2 are identical
            PossL21=PossT
            PossL22=PossT
            PossL23=PossT
            PossL24=PossT
            PossT=PossT+length
            NL2=1
          else
c         L21=L23 and L22=L24
            PossL21=PossT
            PossL23=PossT
            PossT=PossT+length
            PossL22=PossT
            PossL24=PossT
            PossT=PossT+length
            NL2=2
          end if
        else
          if (NbeGrp.eq.1) then
c         L21=L22 and L23=L24
            PossL21=PossT
            PossL22=PossT
            PossT=PossT+length
            PossL23=PossT
            PossL24=PossT
            PossT=PossT+length
            NL2=2
          else
c         all L2 files are different
            PossL21=PossT
            PossT=PossT+length
            PossL22=PossT
            PossT=PossT+length
            PossL23=PossT
            PossT=PossT+length
            PossL24=PossT
            PossT=PossT+length
            NL2=4
          end if
        end if
c
        PossL2W=PossT
        if (intkey.eq.1) then
c        integral based
          length=0
        else
c        cholesky based
        if ((NaGrp.eq.1).and.(NbeGrp.eq.1)) then
          length=nc*nv*(nv+1)/2
        else
          length=nc*mdGrpa*mdGrpbe
        end if
        if ((no*mdGrpbe).gt.length) then
          length=no*mdGrpbe
        end if
        if (nc*no*mdGrpa.gt.length) then
          length=nc*no*mdGrpa
        end if
        end if
c
        PossT=PossT+length
c
        if (printkey.ge.10) then
        write (6,99) 'DM L2 ',PossL21,PossL22,PossL23,PossL24,
     c                       PossL2W
        end if
c
c       H files
c
        length=no*mdSGrpbe
        if (intkey.eq.1) then
          PossH1=PossT
          PossT=PossT+length
          PossH2=PossT
          PossT=PossT+length
        else
          PossH1=PossT
          PossH2=PossT
        end if
        if (printkey.ge.10) then
        write (6,99) 'DM H  ',PossH1,PossH2
        end if
c
c       M files
c
        length=nc*mdSGrpa*mdSGrpbe
c
        if ((NaGrp.eq.1).and.(NaSGrp.eq.1).and.
     c      (NbeGrp.eq.1).and.(NbeSGrp.eq.1)) then
c       M1=M2
          PossM1=PossT
          PossM2=PossT
          PossT=PossT+length
        else
c       M files are different
          PossM1=PossT
          PossT=PossT+length
          PossM2=PossT
          PossT=PossT+length
        end if
        if (printkey.ge.10) then
        write (6,99) 'DM M  ',PossM1,PossM2
        end if
c
c       W1,W2 files
c
        if (NaGrp*NaSGrp.eq.1) then
c       W1 and W2 are identical
          length=mdSGrpa*mdSGrpbe*mdSGrpa*mdSGrpbe
          if ((no.gt.mdSGrpbe).and.(intkey.eq.1)) then
          length=mdSGrpa*mdSGrpbe*mdSGrpa*no
          end if
          PossW1=PossT
          PossW2=PossT
          PossT=PossT+length
        else
c       W1 and W2 are different
          length=mdSGrpa*mdSGrpbe*mdSGrpa*mdSGrpbe
          if ((no.gt.mdSGrpbe).and.(intkey.eq.1)) then
          length=mdSGrpa*mdSGrpbe*mdSGrpa*no
          end if
          PossW1=PossT
          PossT=PossT+length
          PossW2=PossT
          PossT=PossT+length
        end if
c
c       Ww file
c
        PossWw=PossT
c
        if ((NaGrp.eq.1).and.(NaSGrp.eq.1)) then
          lenab=nv*(nv+1)/2
        else
          lenab=mdSGrpa*mdSGrpa
        end if
        if ((NbeGrp.eq.1).and.(NbeSGrp.eq.1)) then
          lenbega=nv*(nv+1)/2
        else
          lenbega=mdSGrpbe*mdSGrpbe
        end if
c
        length=lenab*lenbega
        if (intkey.eq.1) then
          length=mdSGrpa*mdSGrpbe*mdSGrpa*mdSGrpbe
          if (no.gt.mdSGrpbe) then
          length=mdSGrpa*mdSGrpbe*mdSGrpa*no
          end if
        end if
        PossT=PossT+length
c
c        Wx file
c
        PossWx=PossT
        length=mdSGrpa*mdSGrpbe*mdSGrpa*no
        if (intkey.eq.0) then
          length=0
        end if
        PossT=PossT+length
c
        if (printkey.ge.10) then
        write (6,99) 'DM W  ',PossW1,PossW2,PossWw,PossWx
c
        write (6,99) 'PossT ',PossT
        end if
c
99        format (a7,10(i10,1x))
c
        return
        end
