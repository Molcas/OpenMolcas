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
        subroutine InsReqo2v4 (NaGrp,NbeGrp,NaSGrp,NbeSgrp,
     c                   mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe)
c
c       this routine do:
c        Inspect W3 and W4 files requirements of o2v4 procedure
c        on this node. It checks which of the W3 and W4 files
c        are used on this node
c
        use Para_Info, only: MyRank
        implicit none
#include "chcc1.fh"
#include "parcc.fh"
#include "o2v4.fh"
c
        integer NaGrp,NbeGrp,NaSGrp,NbeSgrp
        integer mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe
c
c        help variables
c
        integer aGrp,bGrp,gaGrp,beGrp,aSGrp,bSGrp,gaSGrp,beSGrp
        integer dima,dimb,adda,addb
        integer dimbe,dimga,addbe,addga,addbepp,addgapp
        integer bSGrpUp,gaSGrpUp
        integer LenW3,LenW4
        integer i,j,NSGrp
c
c
c*        Initial Set InqW3,InqW4 = F
        NSGrp=NaGrp*NaSGrp
        do i=1,(NSGrp*(NSGrp+1))/2
          do j=1,NSGrp
          InqW3(i,j)=.False.
          end do
          do j=1,(NSGrp*(NSGrp+1))/2
          InqW4(i,j)=.False.
          end do
        end do
        LenW3=0
        LenW4=0
c
c
c*      cycle over all groups of (a>=b)
        adda=0
        do aGrp=1,naGrp
        dima=DimGrpa(aGrp)
c
c##        test, if on this node atleast one combination with this
c        aGrp is scheduled. Skip if no
        i=0
        do j=1,NaGrp
          i=i+ABID(myRank,aGrp,j)
        end do
        if (i.eq.0) goto 12
c
c
        addb=0
        do bGrp=1,aGrp
        dimb=DimGrpa(bGrp)
c
c##     test, if this a'b' combination is planed to be run on this node
        if (ABID(myRank,aGrp,bGrp).eq.0) goto 11
c
c
cx        cycle over all groups of (be>=ga)
          addbe=0
          do beGrp=1,nbeGrp
          dimbe=DimGrpbe(beGrp)
          addga=0
          do gaGrp=1,beGrp
          dimga=DimGrpbe(gaGrp)
c
c
cxx         cycle over all subgroups of (be>=ga)'
            addbepp=addbe
            do beSGrp=GrpbeLow(beGrp),GrpbeUp(beGrp)
c
            if (beGrp.eq.gaGrp) then
              gaSGrpUp=beSGrp
            else
              gaSGrpUp=GrpbeUp(gaGrp)
            end if
c
            addgapp=addga
            do gaSGrp=GrpbeLow(gaGrp),gaSGrpUp
c
c
cxxx          cycle over all subgroups of (a>=b)'
              do aSGrp=GrpaLow(aGrp),GrpaUp(aGrp)
              if (aGrp.eq.bGrp) then
                bSGrpUp=aSGrp
              else
                bSGrpUp=GrpaUp(bGrp)
              end if
              do bSGrp=GrpaLow(bGrp),bSGrpUp
c
c
c*              clasical reading or generating of (VV|VV) integrals
c
cxxxx                term W1(a",be",b",ga") <- -(b"ga"|a"i) . T1(i,be")
cxxxx1          Get Ww(b",ga",a",i) (Wx used as Aux)
                call InsReaW3 (bSGrp,gaSGrp,aSGrp,LenW3)
c
cxxxx           Upgrade W1(a",be",b",ga") <<- (a",be"|b",ga")
                call InsReaW4 (aSGrp,beSGrp,bSGrp,gaSGrp,LenW4)
c
c
                if (aSGrp.ne.bSGrp) then
cxxxx                  term W2(b",be",a",ga") <- -(a"ga"|b"i) . T1(i,ga")
cxxxx1            Get Ww(a",ga",b",i) (Wx used as Aux)
                  call InsReaW3 (aSGrp,gaSGrp,bSGrp,LenW3)
c
cxxxx                  term W2(b",be",a",ga") <<- -(b"be"|a"i) . T1(i,ga")
cxxxx1            Get Ww(b",be",a",i) (Wx used as Aux)
                  call InsReaW3 (bSGrp,beSGrp,aSGrp,LenW3)
c
cxxxx             Add W2(b",be",a",ga") <<- (b",be"Ia",ga")
c                  Upgrade:W2, destroy:Ww
                  call InsReaW4 (bSGrp,beSGrp,aSGrp,gaSGrp,LenW4)
                end if
c
c
c
cxxx          end cycle over all subgroups of (a>=b)'
              end do
              end do
c
cxx         end cycle over all subgroups of (be>=ga)'
            addgapp=addgapp+DimSGrpbe(gaSGrp)
            end do
            addbepp=addbepp+DimSGrpbe(beSGrp)
            end do
c
cx        end cycle over all groups of (be>=ga)
          addga=addga+dimga
          end do
          addbe=addbe+dimbe
          end do
c
c*      end cycle over all groups of (a>=b)
11        addb=addb+dimb
        end do
12        adda=adda+dima
        end do
c
c
        return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(NbeSgrp)
        call Unused_integer(mdGrpa)
        call Unused_integer(mdGrpbe)
        call Unused_integer(mdSGrpa)
        call Unused_integer(mdSGrpbe)
      end if
        end
