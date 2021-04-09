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
        subroutine DefParo2v4 (NaGrp,NbeGrp,NaSGrp,NbeSgrp,
     c                         mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe)
c
c       This routine do:
c       define parameters in o2v4.fh using NaGrp,NbeGrp,NaSGrp,NbeSgrp
c
c       I/O parameter description:
c       NaGrp    - # of groups in a set (I)
c       NbeGrp   - # of groups in be set (I)
c       NaSGrp   - # of subgroups in each (a)' group (I)
c       NbeSGrp  - # of subgroups in each (be)' group (I)
c       mdGrpa   - # maximal dimension od (a)' group (O)
c       mdGrpbe  - # maximal dimension od (be)' group (O)
c       mdSGrpa  - # maximal dimension od (a)" subgroup (O)
c       mdSGrpbe - # maximal dimension od (be)" subgroup (O)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
#include "chcc_files.fh"
c
        integer NaGrp,NbeGrp,NaSGrp,NbeSgrp
        integer mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe
c
c       help variables
c
        real*8 rdim
        integer i,j,ij
c
c
c1      define parameters of Groups of a(b) set
c
        rdim=1.0d0*nv/(1.0d0*NaGrp)
c
        Grpalow(1)=1
        Grpaup(1)=NaSGrp
c
        do i=1,NaGrp
c
           if (i.eq.1) then
             UpGrpa(i)=int(rdim*i)
             LowGrpa(i)=1
           else if (i.eq.NaGrp) then
             UpGrpa(i)=nv
             LowGrpa(i)=UpGrpa(i-1)+1
           else
             UpGrpa(i)=int(rdim*i)
             LowGrpa(i)=UpGrpa(i-1)+1
           end if
c
           DimGrpa(i)=(UpGrpa(i)-LowGrpa(i))+1
c
           if (i.gt.1) then
             Grpalow(i)=Grpalow(i-1)+NaSGrp
             Grpaup(i)=Grpaup(i-1)+NaSGrp
           end if
c
        end do
c
c
c2      define parameters of Groups of be(ga) set
c
        rdim=1.0d0*nv/(1.0d0*NbeGrp)
c
        Grpbelow(1)=1
        Grpbeup(1)=NbeSGrp
c
        do i=1,NbeGrp
c
           if (i.eq.1) then
             UpGrpbe(i)=int(rdim*i)
             LowGrpbe(i)=1
           else if (i.eq.NbeGrp) then
             UpGrpbe(i)=nv
             LowGrpbe(i)=UpGrpbe(i-1)+1
           else
             UpGrpbe(i)=int(rdim*i)
             LowGrpbe(i)=UpGrpbe(i-1)+1
           end if
c
           DimGrpbe(i)=(UpGrpbe(i)-LowGrpbe(i))+1
c
           if (i.gt.1) then
             Grpbelow(i)=Grpbelow(i-1)+NbeSGrp
             Grpbeup(i)=Grpbeup(i-1)+NbeSGrp
           end if
c
        end do
c
c
c3      define parameters of SubGroups of a(b) set
c
        ij=0
        do j=1,NaGrp

        rdim=1.0d0*DimGrpa(j)/(1.0d0*NaSGrp)
c
        do i=1,NaSGrp
        ij=ij+1
c
           if (i.eq.1) then
             UpSGrpa(ij)=(LowGrpa(j)-1)+int(rdim*i)
             LowSGrpa(ij)=LowGrpa(j)
           else if (i.eq.NaSGrp) then
             UpSGrpa(ij)=UpGrpa(j)
             LowSGrpa(ij)=UpSGrpa(ij-1)+1
           else
             UpSGrpa(ij)=(LowGrpa(j)-1)+int(rdim*i)
             LowSGrpa(ij)=UpSGrpa(ij-1)+1
           end if
c
          DimSGrpa(ij)=(UpSGrpa(ij)-LowSGrpa(ij))+1
c
        end do
        end do
c
c
c4      define parameters of SubGroups of be(ga) set
c
        ij=0
        do j=1,NbeGrp

        rdim=1.0d0*DimGrpbe(j)/(1.0d0*NbeSGrp)
c
        do i=1,NbeSGrp
        ij=ij+1
c
           if (i.eq.1) then
             UpSGrpbe(ij)=(LowGrpbe(j)-1)+int(rdim*i)
             LowSGrpbe(ij)=LowGrpbe(j)
           else if (i.eq.NbeSGrp) then
             UpSGrpbe(ij)=UpGrpbe(j)
             LowSGrpbe(ij)=UpSGrpbe(ij-1)+1
           else
             UpSGrpbe(ij)=(LowGrpbe(j)-1)+int(rdim*i)
             LowSGrpbe(ij)=UpSGrpbe(ij-1)+1
           end if
c
          DimSGrpbe(ij)=(UpSGrpbe(ij)-LowSGrpbe(ij))+1
c
        end do
        end do
c
c
c5      find maximal dimensions mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe
c
        mdGrpa=DimGrpa(1)
        do i=1,NaGrp
          if (DimGrpa(i).gt.mdGrpa) then
          mdGrpa=DimGrpa(i)
          end if
        end do
c
        mdGrpbe=DimGrpbe(1)
        do i=1,NbeGrp
          if (DimGrpbe(i).gt.mdGrpbe) then
          mdGrpbe=DimGrpbe(i)
          end if
        end do
c
        mdSGrpa=DimSGrpa(1)
        do ij=1,NaGrp*NaSGrp
          if (DimSGrpa(ij).gt.mdSGrpa) then
          mdSGrpa=DimSGrpa(ij)
          end if
        end do
c
        mdSGrpbe=DimSGrpbe(1)
        do ij=1,NbeGrp*NbeSGrp
          if (DimSGrpbe(ij).gt.mdSGrpbe) then
          mdSGrpbe=DimSGrpbe(ij)
          end if
        end do
c
c
c6      def L2Names
c
        do i=1,NaGrp
        do j=1,NbeGrp
          call DefParo3v3Hlp1(i,j,'L2',L2Name(i,j))
        end do
        end do
c
c7      def Tmp3Names
c
        do i=1,MaxSGrp
        do j=1,MaxSGrp
          call DefParo3v3Hlp1(i,j,'X3',Tmp3Name(i,j))
        end do
        end do
c
c
        return
        end
