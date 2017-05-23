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
c
        subroutine summary (wrk,wrksize,NvGrp,LunAux,maxdim,E1,E2,E2os)
c
c        thsi routine do:
c        T1
c          - make T1o = T1n/Den
c          - calc E1 = FTau
c        T2:
c          - Join T2 form Tmp2files-(X,Y) and Tmp3files
c          - make T2o in T2files = T2n/Den
c          - make Tau from new amplitudes
c          - calc E2,E2os = WTau
c
c        Scheme of recreation of T2n:
c
c        1) contributions from X and Y matrices
c        T2n(be',ga',u,v) <<-
c         C1                + 1/2 X(be',u,ga',v)
c        C2                + 1/2 X(ga',v,be',u)
c         C3                - 1/2 Y(be',u,ga',v)
c        C4                - 1/2 Y(ga',v,be',u)
c         C5                - 1   Y(ga',u,be',v)
c        C6                - 1   Y(be',v,ga',u)
c
c        2) contributions from T2n1,T2n2 (T2+,T2-)
c        T2n(be',ga',u,v) <<-
c         C7                    T2+(be'ga',uv)
c        C8                    T2-(be'ga',uv)
c
c        Energy contributions:
c        E1  = sum(a,i)       [  2 Fvo(a,i) . T1(a,i)               ]
c        E2  = sum(a,b,i,j)   [ {2(ai|bj) - (aj|bi)} . Tau(a,b,i,j) ]
c        E2os= sum(a,b,i,j)   [   (ai|bj)            . Tau(a,b,i,j) ]
c
c        Memory requirements:
c        V1-3 - {v'ov'o}
c        H1,2 - {v'o}
c
        implicit none
#include "wrk.fh"
#include "chcc1.fh"
#include "o3v3.fh"
#include "o2v4.fh"
#include "chcc_files.fh"
#include "chcc_parcc.fh"
#include "para_info.fh"
c
        integer NvGrp,LunAux,maxdim
        real*8 E1,E2,E2os
c
c        help variables
        integer beGrp,gaGrp,dimbe,dimga,addbe,addga,dim1,dim2
        integer beSGrp,gaSGrp,addbepp,addgapp,dimbepp,dimgapp
        real*8 Ehlp,Eoshlp
        character*6 LunName
c
        integer PossH1,PossH2
        integer PossV1,PossV2,PossV3
        integer PossT
c
c*        Distribute Memory
        PossT=PossFree
        call DistMemSum (NvGrp,maxdim,
     c       PossV1,PossV2,PossV3,
     c       PossH1,PossH2,
     c       PossT)
c
c
c*        Operations on T1 amplitudes
c
c        Divide T1n by denominators
        call T1_div (wrk(PossT1n),wrk(PossOE),no,nv)
c
c        Set T1o <- T1n
        dim1=no*nv
        call mv0u (dim1,wrk(PossT1n),1,wrk(PossT1o),1)
c
c        Calc E1
        call Energy_E1 (wrk(PossT1n),wrk(PossFvo),no,nv,EHlp)
        E1=2.0d0*EHlp
cmp        E1=2.0d0*E1
c
c
c
c*        Operations on T2 amplitudes
c
c        cycle over be'=ga'
c
        E2=0.0d0
        E2os=0.0d0
c
        addbe=0
        do beGrp=1,NvGrp
        dimbe=DimGrpv(beGrp)
c
c         Extract H1(be',I) <- T1(be,I)
          call ExtT1 (wrk(PossH1),wrk(PossT1o),dimbe,addbe)
c
c          vanish V3
          dim1=dimbe*dimbe*no*no
          call mv0zero (dim1,dim1,wrk(PossV3))

c
c1.1          Add contribution C1-6
c
c          read V1(be',u,ga',v) <- X(be',u,ga',v)
c               V2(be',u,ga',v) <- Y(be',u,ga',v)
          LunName=Tmp2Name(beGrp,beGrp)
          dim1=dimbe*dimbe*no*no
          if (XYyes(beGrp,beGrp).eq.1) then
            call GetX (wrk(PossV1),dim1,LunAux,LunName,1,0)
            call GetX (wrk(PossV2),dim1,LunAux,LunName,0,1)
          else if (Xyes(beGrp,beGrp).eq.1) then
            call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
            call mv0zero (dim1,dim1,wrk(PossV2))
          else
            call mv0zero (dim1,dim1,wrk(PossV1))
            call mv0zero (dim1,dim1,wrk(PossV2))
          end if

c
c          V3(be',ga',u,v) <-
c            C1                + 1/2 V1(be',u,ga',v)
c           C2                + 1/2 V1(ga',v,be',u)
c            C3                - 1/2 V2(be',u,ga',v)
c           C4                - 1/2 V2(ga',v,be',u)
c            C5                - 1   V2(ga',u,be',v)
c           C6                - 1   V2(be',v,ga',u)
          call MkT_CAlld (wrk(PossV3),wrk(PossV1),wrk(PossV2),
     c                    dimbe,no)
c1.2          Add contribution C7,C8
c
c         cycle over all subgroups of (be">=ga")
c@c
c        goto 9911
c@c
          addbepp=0
          do beSGrp=GrpbeLow(beGrp),GrpbeUp(beGrp)
          dimbepp=DimSGrpbe(beSGrp)
          addgapp=0
          do gaSGrp=GrpbeLow(beGrp),beSGrp
          dimgapp=DimSGrpbe(gaSGrp)
c
          if (beSGrp.eq.gaSGrp) then
c            read V1(bega",uv+) <- T2+(bega",uv+)
c                 V2(bega",uv-) <- T2-(bega",uv-)
            LunName=Tmp3Name(beSGrp,beSGrp)
            dim1=dimbepp*(dimbepp+1)*no*(no+1)/4
            dim2=dimbepp*(dimbepp-1)*no*(no-1)/4
            if (T2o2v4yes(beSGrp,beSGrp).eq.1) then
              call GetX (wrk(PossV1),dim1,LunAux,LunName,1,0)
              call GetX (wrk(PossV2),dim2,LunAux,LunName,0,1)
            else
              call mv0zero (dim1,dim1,wrk(PossV1))
              call mv0zero (dim2,dim2,wrk(PossV2))
            end if

c            V3(be',ga',u,v) <<- T2+(bega",uv)
c                               T2-(bega",uv)
            call MkT_C78d (wrk(PossV3),wrk(PossV1),wrk(PossV2),
     c                      dimbe,dimbepp,addbepp,no)
          else
            LunName=Tmp3Name(beSGrp,gaSGrp)
c            read V1(be",ga",uv+) <- T2+(be",ga",uv+)
c                 V2(be",ga",uv-) <- T2-(be",ga",uv-)
            dim1=dimbepp*dimgapp*no*(no+1)/2
            dim2=dimbepp*dimgapp*no*(no-1)/2
            if (T2o2v4yes(beSGrp,gaSGrp).eq.1) then
              call GetX (wrk(PossV1),dim1,LunAux,LunName,1,0)
              call GetX (wrk(PossV2),dim2,LunAux,LunName,0,1)
            else
              call mv0zero (dim1,dim1,wrk(PossV1))
              call mv0zero (dim2,dim2,wrk(PossV2))
            end if
c            V3(be',ga',u,v) <<- T2+(be",ga",uv)
c                               T2-(be",ga",uv)
            call MkT_C78od (wrk(PossV3),wrk(PossV1),wrk(PossV2),
     c                   dimbe,dimbe,dimbepp,dimgapp,addbepp,addgapp,no)
          end if
c
          addgapp=addgapp+dimgapp
          end do
          addbepp=addbepp+dimbepp
          end do
c@c
c9911        continue
c@c
c
c
c1.3          T2 - V3(be',ga',u,v) are now reconstructed (for be'>=ga')
c          Divide by denominators and completed to be'<ga' also
           call T2d_div (wrk(PossV3),wrk(PossOE),
     c                  dimbe,dimbe,addbe,addbe,no,nv)
c
c          Prepair reduced Set V2(be'>=ga',u,v) <- V3(be',ga',u,v)
          call MkT_red (wrk(PossV2),wrk(PossV3),dimbe,no)
c
#ifdef _MOLCAS_MPP_
c##          Synchronizacny bod:
c          Allreduce V2
           dim1=dimbe*(dimbe+1)*no*no/2
          call gadgop (wrk(PossV2),dim1,'+')
c
c##          For parallel runs only:
c          Inverse expansion  V3(be',ga',u,v) <- V2(be'>=ga',u,v)
           if (nProcs.gt.1) then
             call MkT_exp (wrk(PossV2),wrk(PossV3),dimbe,no)
           end if
#endif
c@@
c        call Chck_T2n (wrk(PossV3),dimbe,addbe,dimbe,addbe,1)
c@@
c
c          Save into corresponding T2file
          LunName=T2Name(beGrp,beGrp)
          dim1=dimbe*(dimbe+1)*no*no/2
          call SaveX (wrk(PossV2),dim1,LunAux,LunName,1,1)
c
c          Make Tau in V3(be',ga',u,v)
c          (note, that T1 are already new ones)
c          @@ s tym 2x H1 to nieje celkom koser
               call MkTau_chcc (wrk(PossV3),wrk(PossH1),wrk(PossH1),
     c                     dimbe,dimbe,no,1.0d0,1.0d0)
c
c          read V1(be',I,ga',J) <- I2(be'I|ga'J)
           LunName=I2Name(beGrp,beGrp)
          dim1=dimbe*dimbe*no*no
          call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
c          Calc E2 = sum(a,b,i,j) (2V1(ai|bj)-V1(aj|bi)) . V3(a,b,i,j)
          call Energy_E2d (wrk(PossV1),wrk(PossV3),Ehlp,Eoshlp,dimbe,no)

          E2=E2+1.0d0*Ehlp
          E2os=E2os+1.0d0*Eoshlp
c
        addbe=addbe+dimbe
        end do
c
c
c        cycle over be>ga
c        reserved files:
c          V3(be',ga',u,v) = T2n(be',ga',u,v)
c        free for use: V1,V2
c
        addbe=DimGrpv(1)
        do beGrp=2,NvGrp
        dimbe=DimGrpv(beGrp)
c
c         Extract H1(be',I) <- T1(be,I)
          call ExtT1 (wrk(PossH1),wrk(PossT1o),dimbe,addbe)
c
c
          addga=0
          do gaGrp=1,beGrp-1
          dimga=DimGrpv(gaGrp)
c
c           Extract H2(ga',I) <- T1(ga,I)
            call ExtT1 (wrk(PossH2),wrk(PossT1o),dimga,addga)
c
c            vanish V3
            dim1=dimbe*dimga*no*no
            call mv0zero (dim1,dim1,wrk(PossV3))
c
c            read V1(be',u,ga',v) <- X(be',u,ga',v)
c                 V2(be',u,ga',v) <- Y(be',u,ga',v)
            LunName=Tmp2Name(beGrp,gaGrp)
            dim1=dimbe*dimga*no*no
            if (XYyes(beGrp,gaGrp).eq.1) then
              call GetX (wrk(PossV1),dim1,LunAux,LunName,1,0)
              call GetX (wrk(PossV2),dim1,LunAux,LunName,0,1)
            else if (Xyes(beGrp,gaGrp).eq.1) then
              call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
              call mv0zero (dim1,dim1,wrk(PossV2))
            else
              call mv0zero (dim1,dim1,wrk(PossV1))
              call mv0zero (dim1,dim1,wrk(PossV2))
            end if
c
c            Add contribution C1,C3,C6
c            V3(be',ga',u,v) < - + 1/2 V1(be',u,ga',v)
c                                - 1/2 V2(be',u,ga',v)
c                               - 1   V2(be',v,ga',u)
            call MkT_C136od (wrk(PossV3),wrk(PossV1),wrk(PossV2),
     c                       dimbe,dimga,no)
c
c            read V1(ga',u,be',v) <- X(ga',u,be',v)
c                 V2(ga',u,be',v) <- Y(ga',u,be',v)
            LunName=Tmp2Name(gaGrp,beGrp)
            dim1=dimbe*dimga*no*no
            if (XYyes(gaGrp,beGrp).eq.1) then
              call GetX (wrk(PossV1),dim1,LunAux,LunName,1,0)
              call GetX (wrk(PossV2),dim1,LunAux,LunName,0,1)
            else if (Xyes(gaGrp,beGrp).eq.1) then
              call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
              call mv0zero (dim1,dim1,wrk(PossV2))
            else
              call mv0zero (dim1,dim1,wrk(PossV1))
              call mv0zero (dim1,dim1,wrk(PossV2))
            end if
c
c            Add contribution C2,C4,C5
c            V3(be',ga',u,v) <<- + 1/2 V1(ga',v,be',u)
c                               - 1/2 V2(ga',v,be',u)
c                                - 1   V2(ga',u,be',v)
            call MkT_C245od (wrk(PossV3),wrk(PossV1),wrk(PossV2),
     c                       dimbe,dimga,no)
c
c
c2.2            Add contribution C7,C8
c
c           cycle over all subgroups of (be">=ga")
            addbepp=0
            do beSGrp=GrpbeLow(beGrp),GrpbeUp(beGrp)
            dimbepp=DimSGrpbe(beSGrp)
            addgapp=0
            do gaSGrp=GrpbeLow(gaGrp),GrpbeUp(gaGrp)
            dimgapp=DimSGrpbe(gaSGrp)
c
c              read V1(be",ga",uv+) <- T2+(be",ga",uv+)
c                   V2(be",ga",uv-) <- T2-(be",ga",uv-)
              LunName=Tmp3Name(beSGrp,gaSGrp)
              dim1=dimbepp*dimgapp*no*(no+1)/2
              dim2=dimbepp*dimgapp*no*(no-1)/2
              if (T2o2v4yes(beSGrp,gaSGrp).eq.1) then
                call GetX (wrk(PossV1),dim1,LunAux,LunName,1,0)
                call GetX (wrk(PossV2),dim2,LunAux,LunName,0,1)
              else
                call mv0zero (dim1,dim1,wrk(PossV1))
                call mv0zero (dim2,dim2,wrk(PossV2))
              end if
c              V3(be',ga',u,v) <<- T2+(be",ga",uv)
c                                 T2-(be",ga",uv)
              call MkT_C78od (wrk(PossV3),wrk(PossV1),wrk(PossV2),
     c                   dimbe,dimga,dimbepp,dimgapp,addbepp,addgapp,no)
c
            addgapp=addgapp+dimgapp
            end do
            addbepp=addbepp+dimbepp
              end do
c
c
c
c            T2 - V3(be',ga',u,v) are now completely reconstructed
c            Divide by denominators
             call T2od_div (wrk(PossV3),wrk(PossOE),
     c                     dimbe,dimga,addbe,addga,no,nv)
c
c
#ifdef _MOLCAS_MPP_
c##            Synchronizacny bod:
c            Allreduce V3
            dim1=dimbe*dimga*no*no
            call gadgop (wrk(PossV3),dim1,'+')
#endif
c@@
c        call Chck_T2n (wrk(PossV3),dimbe,addbe,dimga,addga,0)
c@@
c
c
c            Save V3 into corresponding T2file
            LunName=T2Name(beGrp,gaGrp)
            dim1=dimbe*dimga*no*no
            call SaveX (wrk(PossV3),dim1,LunAux,LunName,1,1)
c
c            Map V1(ga',be',v,u) <- V3(be',ga',u,v)
            call Map4_2143 (wrk(PossV3),wrk(PossV1),dimbe,dimga,no,no)
c
c            Save V1 into corresponding T2file
            LunName=T2Name(gaGrp,beGrp)
            dim1=dimbe*dimga*no*no
            call SaveX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
c
c            Make Tau in V3(be',ga',u,v)
c            (note, that T1 are already new ones)
               call MkTau_chcc (wrk(PossV3),wrk(PossH1),wrk(PossH2),
     c                       dimbe,dimga,no,1.0d0,1.0d0)
c
c            read V1(be',I,ga',J) <- I2(be'I|ga'J)
            LunName=I2Name(beGrp,gaGrp)
            dim1=dimbe*dimga*no*no
            call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
c            Calc E2 = sum(a,b,i,j) (2V1(ai|bj)-V1(aj|bi)) . V3(a,b,i,j)
            call Energy_E2od (wrk(PossV1),wrk(PossV3),Ehlp,Eoshlp,
     c                        dimbe,dimga,no)

            E2=E2+2.0d0*Ehlp
            E2os=E2os+2.0d0*Eoshlp
c
          addga=addga+dimga
          end do
c
        addbe=addbe+dimbe
        end do
c
c@@
         return
c
c        Calc Energy from Chck
         call Chck_energ
c
c        upgrade checkeroo T1c and T2c
c
c        Upgrade T1
         call UpG_T1 (wrk(PossT1o))
c
c        Upgrade T2
        addbe=0
        do beGrp=1,NvGrp
        dimbe=DimGrpv(beGrp)
c
          addga=0
          do gaGrp=1,beGrp
          dimga=DimGrpv(gaGrp)
c
          if (beGrp.eq.gaGrp) then
            LunName=T2Name(beGrp,beGrp)
            dim1=dimbe*(dimbe+1)*no*no/2
            call GetX (wrk(PossV3),dim1,LunAux,LunName,1,1)
            call UpG_T2d (wrk(PossV3),dimbe,addbe)
          else
            LunName=T2Name(beGrp,gaGrp)
            dim1=dimbe*dimga*no*no
            call GetX (wrk(PossV3),dim1,LunAux,LunName,1,1)
            call UpG_T2od (wrk(PossV3),dimbe,addbe,dimga,addga)
          end if
c
          addga=addga+dimga
          end do
c
        addbe=addbe+dimbe
        end do
c
        return
        end
c
c        -----------------------------
c
        subroutine Chck_T2n (T2,dimbe,addbe,dimga,addga,key)
c
c        chek T2n bez menovatelov, nediagonalne
c
        implicit none
        integer dimbe,addbe,dimga,addga,key
#include "chcc1.fh"

         real*8 T2(1:dimbe,1:dimga,1:no,1:no)
c
        integer bstart,bup
        integer u,v,be,ga,bad,ntot
        integer i,j,a,b
        real*8 s,s1
c
        bad=0
        ntot=0
c
        do v=1,no
        do u=1,no
        do ga=addga+1,addga+dimga
        if (key.eq.1) then
             bstart=ga
          bup=addbe+dimbe
        else
          bstart=addbe+1
          bup=addbe+dimbe
        end if
        do be=bstart,bup
c
        ntot=ntot+1
        s=0.0d0
c
c1
          s1=Q21(be,u,ga,v)
          s=s+s1
c
c2
          s1=0.0d0
          do j=1,no
          do i=1,no
          s1=s1+Ac(i,j,u,v)*(T2c(be,ga,i,j)+T1c(be,i)*T1c(ga,j))
          end do
          end do
          s=s+s1
c
c3
          s1=0.0d0
          do b=1,nv
          do a=1,nv
           s1=s1+Bc(a,b,be,ga)*(T2c(a,b,u,v)+T1c(a,u)*T1c(b,v))
          end do
          end do
          s=s+s1
c
c4
          s1=0.0d0
          do a=1,nv
          s1=s1+Gvvc(be,a)*T2c(a,ga,u,v)
          s1=s1+Gvvc(ga,a)*T2c(a,be,v,u)
          end do
          s=s+s1
c
c5
          s1=0.0d0
          do i=1,no
          s1=s1+Gooc(i,u)*T2c(be,ga,i,v)
          s1=s1+Gooc(i,v)*T2c(ga,be,i,u)
          end do
          s=s-s1
c
c6
          s1=0.0d0
          do a=1,nv
          s1=s1+Q3(ga,a,be,u)*T1c(a,v)
          s1=s1+Q3(be,a,ga,v)*T1c(a,u)
          end do
          do i=1,no
          do a=1,nv
          s1=s1-(Q22(ga,a,i,u)*T1c(be,i))*T1c(a,v)
          s1=s1-(Q22(be,a,i,v)*T1c(ga,i))*T1c(a,u)
          end do
          end do
          s=s+s1
c
c7
          s1=0.0d0
          do i=1,no
          s1=s1+Q1(be,u,i,v)*T1c(ga,i)
          s1=s1+Q1(ga,v,i,u)*T1c(be,i)
          end do
          do i=1,no
          do a=1,nv
          s1=s1+(Q21(be,u,a,i)*T1c(a,v))*T1c(ga,i)
          s1=s1+(Q21(ga,v,a,i)*T1c(a,u))*T1c(be,i)
          end do
          end do
          s=s-s1
c
c8
          s1=0.0d0
          do i=1,no
          do a=1,nv
            s1=s1+(2.0d0*Jc(be,i,u,a)-Kc(i,be,u,a))*
     c          (2.0d0*T2c(a,ga,i,v)-T2c(ga,a,i,v))
            s1=s1+(2.0d0*Jc(ga,i,v,a)-Kc(i,ga,v,a))*
     c          (2.0d0*T2c(a,be,i,u)-T2c(be,a,i,u))

          end do
          end do
          s=s+0.5d0*s1
c
c9
          s1=0.0d0
          do i=1,no
          do a=1,nv
          s1=s1+Kc(i,be,u,a)*T2c(ga,a,i,v)
          s1=s1+Kc(i,ga,v,a)*T2c(be,a,i,u)
          end do
          end do
          s=s-0.5d0*s1
c
c10
          s1=0.0d0
          do i=1,no
          do a=1,nv
          s1=s1+Kc(i,ga,u,a)*T2c(be,a,i,v)
          s1=s1+Kc(i,be,v,a)*T2c(ga,a,i,u)
          end do
          end do
          s=s-s1
c
          s=s/(Oeo(u)+Oeo(v)-Oev(be)-Oev(ga))
c
             if  (abs(T2(be-addbe,ga-addga,u,v)-s).gt.1.0d-10) then
c          write (6,*) 'Bad', abs(T2(be-addbe,ga-addga,u,v)
          bad=bad+1
          end if
c
        end do
        end do
        end do
        end do
c
        if (key.eq.1) then
          write (6,*) ' Final test T2 dia',bad,ntot
        else
          write (6,*) ' Final test T2 off',bad,ntot
        end if
c
        return
        end
