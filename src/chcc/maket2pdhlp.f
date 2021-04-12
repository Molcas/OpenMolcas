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
        subroutine makeT2pdHlp (T2p,Tau,aGrp,aSGrp,
     c                          dimi,dimij,dimapp,dimap,dimabp)
c
c       this routine do:
c       define T2(+)((aa)",ij) = Tau((ab)',i,j) + Tau((ab)',j,i)
c       for the case: aGrp=bGrp and aSGrp=bSGrp
c
c       parameter description:
c       T2p     - array for T2+ (O)
c       Tau     - array for Tau (I)
c       xGrp    - Groups of a',b' (I)
c       xSGrp   - SubGroups of a",b" (I)
c       dimx    - Dimension of i,(i>=j),a",a',(a>=b)' (I)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        integer aGrp,aSGrp
        integer dimi,dimij,dimapp,dimap,dimabp
        real*8 T2p(1:dimapp,1:dimij)
        real*8 Tau(1:dimabp,1:dimi,1:dimi)
c
c
c       help variables
        integer i,j,ij,app,ap,abp,appAdd
c
c
c1      def appAdd
c
        appAdd=0
        if (aSGrp.ne.Grpalow(aGrp)) then
          do i=Grpalow(aGrp),aSGrp-1
          appAdd=appAdd+DimSGrpa(i)
          end do
        end if
c
c2        define T2+(aa",ij)
c
        ij=0
        do i=1,dimi
        do j=1,i
        ij=ij+1
          do app=1,dimapp
          ap=appAdd+app
          abp=ap*(ap+1)/2
            T2p(app,ij)=Tau(abp,i,j)+Tau(abp,j,i)
          end do
        end do
        end do
c
        call mv0sv (dimij*dimapp,dimij*dimapp,
     c              T2p(1,1),0.5d0)
c
        return
c Avoid unused argument warnings
        if (.false.) call Unused_integer(dimap)
        end
