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
        subroutine summary (wrk,wrksize,NvGrp,LunAux,maxdim,E1,E2,E2os)
!
!        thsi routine do:
!        T1
!          - make T1o = T1n/Den
!          - calc E1 = FTau
!        T2:
!          - Join T2 form Tmp2files-(X,Y) and Tmp3files
!          - make T2o in T2files = T2n/Den
!          - make Tau from new amplitudes
!          - calc E2,E2os = WTau
!
!        Scheme of recreation of T2n:
!
!        1) contributions from X and Y matrices
!        T2n(be',ga',u,v) <<-
!         C1                + 1/2 X(be',u,ga',v)
!        C2                + 1/2 X(ga',v,be',u)
!         C3                - 1/2 Y(be',u,ga',v)
!        C4                - 1/2 Y(ga',v,be',u)
!         C5                - 1   Y(ga',u,be',v)
!        C6                - 1   Y(be',v,ga',u)
!
!        2) contributions from T2n1,T2n2 (T2+,T2-)
!        T2n(be',ga',u,v) <<-
!         C7                    T2+(be'ga',uv)
!        C8                    T2-(be'ga',uv)
!
!        Energy contributions:
!        E1  = sum(a,i)       [  2 Fvo(a,i) . T1(a,i)               ]
!        E2  = sum(a,b,i,j)   [ {2(ai|bj) - (aj|bi)} . Tau(a,b,i,j) ]
!        E2os= sum(a,b,i,j)   [   (ai|bj)            . Tau(a,b,i,j) ]
!
!        Memory requirements:
!        V1-3 - {v'ov'o}
!        H1,2 - {v'o}
!
#ifdef _MOLCAS_MPP_
        use Para_Info, only: nProcs
#endif
        implicit none
#include "wrk.fh"
#include "chcc1.fh"
#include "o3v3.fh"
#include "o2v4.fh"
#include "chcc_files.fh"
#include "parcc.fh"
!
        integer NvGrp,LunAux,maxdim
        real*8 E1,E2,E2os
!
!        help variables
        integer beGrp,gaGrp,dimbe,dimga,addbe,addga,dim1,dim2
        integer beSGrp,gaSGrp,addbepp,addgapp,dimbepp,dimgapp
        real*8 Ehlp,Eoshlp
        character*6 LunName
!
        integer PossH1,PossH2
        integer PossV1,PossV2,PossV3
        integer PossT
!
!*        Distribute Memory
        PossT=PossFree
        call DistMemSum (NvGrp,maxdim,                                  &
     &       PossV1,PossV2,PossV3,                                      &
     &       PossH1,PossH2,                                             &
     &       PossT)
!
!
!*        Operations on T1 amplitudes
!
!        Divide T1n by denominators
        call T1_div (wrk(PossT1n),wrk(PossOE),no,nv)
!
!        Set T1o <- T1n
        dim1=no*nv
        call mv0u (dim1,wrk(PossT1n),1,wrk(PossT1o),1)
!
!        Calc E1
        call Energy_E1 (wrk(PossT1n),wrk(PossFvo),no,nv,EHlp)
        E1=2.0d0*EHlp
!mp        E1=2.0d0*E1
!
!
!
!*        Operations on T2 amplitudes
!
!        cycle over be'=ga'
!
        E2=0.0d0
        E2os=0.0d0
!
        addbe=0
        do beGrp=1,NvGrp
        dimbe=DimGrpv(beGrp)
!
!         Extract H1(be',I) <- T1(be,I)
          call ExtT1 (wrk(PossH1),wrk(PossT1o),dimbe,addbe)
!
!          vanish V3
          dim1=dimbe*dimbe*no*no
          call mv0zero (dim1,dim1,wrk(PossV3))

!
!1.1          Add contribution C1-6
!
!          read V1(be',u,ga',v) <- X(be',u,ga',v)
!               V2(be',u,ga',v) <- Y(be',u,ga',v)
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

!
!          V3(be',ga',u,v) <-
!            C1                + 1/2 V1(be',u,ga',v)
!           C2                + 1/2 V1(ga',v,be',u)
!            C3                - 1/2 V2(be',u,ga',v)
!           C4                - 1/2 V2(ga',v,be',u)
!            C5                - 1   V2(ga',u,be',v)
!           C6                - 1   V2(be',v,ga',u)
          call MkT_CAlld (wrk(PossV3),wrk(PossV1),wrk(PossV2),          &
     &                    dimbe,no)
!1.2          Add contribution C7,C8
!
!         cycle over all subgroups of (be">=ga")
!@c
!        goto 9911
!@c
          addbepp=0
          do beSGrp=GrpbeLow(beGrp),GrpbeUp(beGrp)
          dimbepp=DimSGrpbe(beSGrp)
          addgapp=0
          do gaSGrp=GrpbeLow(beGrp),beSGrp
          dimgapp=DimSGrpbe(gaSGrp)
!
          if (beSGrp.eq.gaSGrp) then
!            read V1(bega",uv+) <- T2+(bega",uv+)
!                 V2(bega",uv-) <- T2-(bega",uv-)
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

!            V3(be',ga',u,v) <<- T2+(bega",uv)
!                               T2-(bega",uv)
            call MkT_C78d (wrk(PossV3),wrk(PossV1),wrk(PossV2),         &
     &                      dimbe,dimbepp,addbepp,no)
          else
            LunName=Tmp3Name(beSGrp,gaSGrp)
!            read V1(be",ga",uv+) <- T2+(be",ga",uv+)
!                 V2(be",ga",uv-) <- T2-(be",ga",uv-)
            dim1=dimbepp*dimgapp*no*(no+1)/2
            dim2=dimbepp*dimgapp*no*(no-1)/2
            if (T2o2v4yes(beSGrp,gaSGrp).eq.1) then
              call GetX (wrk(PossV1),dim1,LunAux,LunName,1,0)
              call GetX (wrk(PossV2),dim2,LunAux,LunName,0,1)
            else
              call mv0zero (dim1,dim1,wrk(PossV1))
              call mv0zero (dim2,dim2,wrk(PossV2))
            end if
!            V3(be',ga',u,v) <<- T2+(be",ga",uv)
!                               T2-(be",ga",uv)
            call MkT_C78od (wrk(PossV3),wrk(PossV1),wrk(PossV2),        &
     &                   dimbe,dimbe,dimbepp,dimgapp,addbepp,addgapp,no)
          end if
!
          addgapp=addgapp+dimgapp
          end do
          addbepp=addbepp+dimbepp
          end do
!@c
!9911        continue
!@c
!
!
!1.3          T2 - V3(be',ga',u,v) are now reconstructed (for be'>=ga')
!          Divide by denominators and completed to be'<ga' also
           call T2d_div (wrk(PossV3),wrk(PossOE),                       &
     &                  dimbe,dimbe,addbe,addbe,no,nv)
!
!          Prepair reduced Set V2(be'>=ga',u,v) <- V3(be',ga',u,v)
          call MkT_red (wrk(PossV2),wrk(PossV3),dimbe,no)
!
#ifdef _MOLCAS_MPP_
!##          Synchronizacny bod:
!          Allreduce V2
           dim1=dimbe*(dimbe+1)*no*no/2
          call gadgop (wrk(PossV2),dim1,'+')
!
!##          For parallel runs only:
!          Inverse expansion  V3(be',ga',u,v) <- V2(be'>=ga',u,v)
           if (nProcs.gt.1) then
             call MkT_exp (wrk(PossV2),wrk(PossV3),dimbe,no)
           end if
#endif
!@@
!        call Chck_T2n (wrk(PossV3),dimbe,addbe,dimbe,addbe,1)
!@@
!
!          Save into corresponding T2file
          LunName=T2Name(beGrp,beGrp)
          dim1=dimbe*(dimbe+1)*no*no/2
          call SaveX (wrk(PossV2),dim1,LunAux,LunName,1,1)
!
!          Make Tau in V3(be',ga',u,v)
!          (note, that T1 are already new ones)
!          @@ s tym 2x H1 to nieje celkom koser
               call MkTau_chcc (wrk(PossV3),wrk(PossH1),wrk(PossH1),    &
     &                     dimbe,dimbe,no,1.0d0,1.0d0)
!
!          read V1(be',I,ga',J) <- I2(be'I|ga'J)
           LunName=I2Name(beGrp,beGrp)
          dim1=dimbe*dimbe*no*no
          call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
!
!          Calc E2 = sum(a,b,i,j) (2V1(ai|bj)-V1(aj|bi)) . V3(a,b,i,j)
          call Energy_E2d (wrk(PossV1),wrk(PossV3),Ehlp,Eoshlp,dimbe,no)

          E2=E2+1.0d0*Ehlp
          E2os=E2os+1.0d0*Eoshlp
!
        addbe=addbe+dimbe
        end do
!
!
!        cycle over be>ga
!        reserved files:
!          V3(be',ga',u,v) = T2n(be',ga',u,v)
!        free for use: V1,V2
!
        addbe=DimGrpv(1)
        do beGrp=2,NvGrp
        dimbe=DimGrpv(beGrp)
!
!         Extract H1(be',I) <- T1(be,I)
          call ExtT1 (wrk(PossH1),wrk(PossT1o),dimbe,addbe)
!
!
          addga=0
          do gaGrp=1,beGrp-1
          dimga=DimGrpv(gaGrp)
!
!           Extract H2(ga',I) <- T1(ga,I)
            call ExtT1 (wrk(PossH2),wrk(PossT1o),dimga,addga)
!
!            vanish V3
            dim1=dimbe*dimga*no*no
            call mv0zero (dim1,dim1,wrk(PossV3))
!
!            read V1(be',u,ga',v) <- X(be',u,ga',v)
!                 V2(be',u,ga',v) <- Y(be',u,ga',v)
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
!
!            Add contribution C1,C3,C6
!            V3(be',ga',u,v) < - + 1/2 V1(be',u,ga',v)
!                                - 1/2 V2(be',u,ga',v)
!                               - 1   V2(be',v,ga',u)
            call MkT_C136od (wrk(PossV3),wrk(PossV1),wrk(PossV2),       &
     &                       dimbe,dimga,no)
!
!            read V1(ga',u,be',v) <- X(ga',u,be',v)
!                 V2(ga',u,be',v) <- Y(ga',u,be',v)
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
!
!            Add contribution C2,C4,C5
!            V3(be',ga',u,v) <<- + 1/2 V1(ga',v,be',u)
!                               - 1/2 V2(ga',v,be',u)
!                                - 1   V2(ga',u,be',v)
            call MkT_C245od (wrk(PossV3),wrk(PossV1),wrk(PossV2),       &
     &                       dimbe,dimga,no)
!
!
!2.2            Add contribution C7,C8
!
!           cycle over all subgroups of (be">=ga")
            addbepp=0
            do beSGrp=GrpbeLow(beGrp),GrpbeUp(beGrp)
            dimbepp=DimSGrpbe(beSGrp)
            addgapp=0
            do gaSGrp=GrpbeLow(gaGrp),GrpbeUp(gaGrp)
            dimgapp=DimSGrpbe(gaSGrp)
!
!              read V1(be",ga",uv+) <- T2+(be",ga",uv+)
!                   V2(be",ga",uv-) <- T2-(be",ga",uv-)
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
!              V3(be',ga',u,v) <<- T2+(be",ga",uv)
!                                 T2-(be",ga",uv)
              call MkT_C78od (wrk(PossV3),wrk(PossV1),wrk(PossV2),      &
     &                   dimbe,dimga,dimbepp,dimgapp,addbepp,addgapp,no)
!
            addgapp=addgapp+dimgapp
            end do
            addbepp=addbepp+dimbepp
              end do
!
!
!
!            T2 - V3(be',ga',u,v) are now completely reconstructed
!            Divide by denominators
             call T2od_div (wrk(PossV3),wrk(PossOE),                    &
     &                     dimbe,dimga,addbe,addga,no,nv)
!
!
#ifdef _MOLCAS_MPP_
!##            Synchronizacny bod:
!            Allreduce V3
            dim1=dimbe*dimga*no*no
            call gadgop (wrk(PossV3),dim1,'+')
#endif
!@@
!        call Chck_T2n (wrk(PossV3),dimbe,addbe,dimga,addga,0)
!@@
!
!
!            Save V3 into corresponding T2file
            LunName=T2Name(beGrp,gaGrp)
            dim1=dimbe*dimga*no*no
            call SaveX (wrk(PossV3),dim1,LunAux,LunName,1,1)
!
!            Map V1(ga',be',v,u) <- V3(be',ga',u,v)
            call Map4_2143 (wrk(PossV3),wrk(PossV1),dimbe,dimga,no,no)
!
!            Save V1 into corresponding T2file
            LunName=T2Name(gaGrp,beGrp)
            dim1=dimbe*dimga*no*no
            call SaveX (wrk(PossV1),dim1,LunAux,LunName,1,1)
!
!
!            Make Tau in V3(be',ga',u,v)
!            (note, that T1 are already new ones)
               call MkTau_chcc (wrk(PossV3),wrk(PossH1),wrk(PossH2),    &
     &                       dimbe,dimga,no,1.0d0,1.0d0)
!
!            read V1(be',I,ga',J) <- I2(be'I|ga'J)
            LunName=I2Name(beGrp,gaGrp)
            dim1=dimbe*dimga*no*no
            call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
!
!            Calc E2 = sum(a,b,i,j) (2V1(ai|bj)-V1(aj|bi)) . V3(a,b,i,j)
            call Energy_E2od (wrk(PossV1),wrk(PossV3),Ehlp,Eoshlp,      &
     &                        dimbe,dimga,no)

            E2=E2+2.0d0*Ehlp
            E2os=E2os+2.0d0*Eoshlp
!
          addga=addga+dimga
          end do
!
        addbe=addbe+dimbe
        end do
!
!@@
         return
!
!        Calc Energy from Chck
         call Chck_energ
!
!        upgrade checkeroo T1c and T2c
!
!        Upgrade T1
         call UpG_T1 (wrk(PossT1o))
!
!        Upgrade T2
        addbe=0
        do beGrp=1,NvGrp
        dimbe=DimGrpv(beGrp)
!
          addga=0
          do gaGrp=1,beGrp
          dimga=DimGrpv(gaGrp)
!
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
!
          addga=addga+dimga
          end do
!
        addbe=addbe+dimbe
        end do
!
        return
        end
