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
        subroutine o3v3t2 (wrk,wrksize,NvGrp,maxdim,LunAux)
c
c       This routine do:
c
c       --- Part III, generation of X, Y intermediates
c
c        T2(be,ga,u,v) <-
c        T21  +             [ (be,u|ga,v)                 ]
c        T22  + sum(i,j)    [ Ta(be,ga,ij) . A(i,j,u,v)   ]
c          T21 implemented as X0:
c          T22 implemented as X4:
c
c       X(be,u,ga,v) <<-
c          X0 (T21)     1/2 (be,u|ga,v)
c         X1         + Q(be,u,i,a) . E(ga,a,i,v)
c          X2 (T24)   + Gvv(a,be)   . t2(ga,a,v,u)
c          X3 (T25)   - Goo(i,u)    . t2(ga,be,v,i)
c          X4 (T22)   + sum(i>=j)   [ Ta(be,ga,ij) . A(ij,u,v)  ]
c                    - 1/2 sum(i)  [ Ta(be,ga,ii) . A(ii,u,v)  ]
c          Xe (T2ex)  - sum(i,j)    [ Aex(i,j,u,v) . T1(be,i) .T1(ga,j) ]
c                       Xe - o4v2 term is factorized to o4v + o3v2 terms
c
c       Y(be,u,ga,v) = K(be,u,i,a) . t2(i,a,ga,v)
c
c       where E(ga,a,i,v) = 2 t2(ga,a,v,i) - t2(ga,a,i,v)
c
c       T1(be,u) << -
c       T11  +            [  F(be,u)                   ]
c       T12  - sum(a,i)   [ 2F(a,i) . t(be,i) . t(a,u) ]
c       T13  + sum(a)     [  Hvv(a,be) . t(a,u)        ]
c       T14  - sum(i)     [  Hoo(i,u)  . t(be,i)       ]
c       T15  + sum(a,i)   [  Hoo(a,i)  .
c           (2t2(be,a,u,i) - t2(be,a,i,u) + t(a,u).t(be,i)) ]
c
c       Extra cost - multiple reading:
c                    1) N_be . T;
c                    2) N_ga . Q
c                    3) N_ga . K
c
c
c1      intermediates Q and K will be temporarry stored in files
c       QFil, KFil as follows
c       do beGrp
c       do aGrp
c         Q(K) (be',u,i,a')
c       end do
c       end do
c
c2      Structure of files, where selected group of (pq|rs) are
c       stored (V'O|OO) - I1 ; (V'O|V'O) - I2 ; (V'V'|OO) - I3
c
c       (A'I|JK)  I1inxx xx - Group of A'
c
c       (A'I|B'J) I2xxyy xx - Group of A'
c                        yy - Group of B'
c
c       (A'B'|IJ) I3xxyy xx - Group of A'
c                        yy - Group of B'
c
c
c3      Structure of Cholesky vector files
c
c       L1(m,I ,A')  L1vcxx xx - Group of A'
c
c       L2(m,A',B')  L2xxyy xx - Group of A'
c                           yy - Group of B'
c
c4      Structure of Amplitude file
c       t2(A',B',IJ)  T2xxyy xx - Group of A'
c                            yy - Group of B'
c
c
        use Para_Info, only: MyRank
        implicit none
#include "chcc1.fh"
#include "chcc_parcc.fh"
#include "o3v3.fh"
#include "chcc_files.fh"
#include "wrk.fh"
#include "chcc_casy.fh"
c
        integer NvGrp,LunAux

c       help variables
        integer dim1,dim2,dima,dimbe,dimga
        integer aGrp,beGrp,gaGrp,adda,addbe,addga
        integer PossV1,PossV2,PossV3,PossV4
        integer PossH1,PossH2,PossH3,PossH4
        integer PossK,PossQ,PossX,PossY
c
        integer PossT,maxdim

c
        character*6 LunName
c
c       --- introduction part ---
c
cx      Distribute memory
c       @@ yet distribution is the same for all three o3v3 drivers,
c       v buducnosti nosnost usposobit vsetky 3 vlastne
c
        PossT=PossFree
        call DistMemo3v3t2 (NvGrp,maxdim,
     c       PossV1,PossV2,PossV3,PossV4,
     c       PossH1,PossH2,PossH3,PossH4,
     c       PossK,PossQ,
     c       PossT)
c
        if (printkey.ge.10) then
        write (6,*) ' Last Value :',PossT,wrksize
        end if
        if (PossT.gt.wrksize) then
cmp!           write (6,*) ' Nieje dobre, Dr. Ch.  Kokotopuloss',
        write (6,*) ' Not Enough memory in o3v3t2 step!',
     & 'Increase large and/or small segmentation ',
     c                  (1.0d0*PossT)/(1.0d0*wrksize)
           call abend()
        end if
c
c
c##G     Vanish XYyes
c
        do dim1=1,NvGrp
        do dim2=1,NvGrp
          XYyes(dim1,dim2)=0
          if (printkey.ge.10) then
          write (6,*) 'Xyes',dim1,dim2,Xyes(dim1,dim2)
          end if
        end do
        end do
c
c
c
        PossX=PossQ
        PossY=PossK
c
        addbe=0
        do beGrp=1,NvGrp
        dimbe=DimGrpv(beGrp)
c
c##     test, if this beGrp is planed to be run on this node
        dim1=0
        do dim2=1,NvGrp
          dim1=dim1+BeAID(myRank,beGrp,dim2)
        end do
        dim1=dim1+BetaID(myRank,beGrp)
        if (dim1.eq.0) goto 11
        if (printkey.gt.1) then ! nie som si isty ....
        write (6,*) ' o3v3 T2 - ID, beGrp',myRank,beGrp
        end if
cmp
        Call CWTime(TCpu,TWall)
        if (printkey.gt.1) then
        write (6,*)
        write (6,'(A,f18.1)') ' Cpu last call [s] = ',
     & TCpu-TCpu_l
        write (6,'(A,f18.1)') 'Wall last call [s] = ',
     & TWall-TWall_l
        write (6,*)
        write (6,'(A,f18.1)') 'Total Cpu  [s] = ',
     & TCpu
        write (6,'(A,f18.1)') 'Total Wall [s] = ',
     & TWall-TWall0
        write (6,'(A,f18.2)') 'TCpu/TWall [%] = ',
     & 100.0d0*TCpu/(TWall-TWall0)
        write (6,*)
        end if
        TCpu_l=TCpu
        TWall_l=TWall
cmp
c
cT1G          vanish H4(be',u)
          dim1=dimbe*no
          call mv0zero (dim1,dim1,wrk(PossH4))
c
cG             Extract H1(be',I) <- T1(be,I)
          call ExtT1 (wrk(PossH1),wrk(PossT1o),dimbe,addbe)
c
          addga=0
          do gaGrp=1,NvGrp
          dimga=DimGrpv(gaGrp)
c
c           initialize X(be',u,ga',v) as X0, calc in step II
c           vanish Y(be',u,ga',v)
            LunName=Tmp2Name(beGrp,gaGrp)
            dim1=no*no*DimGrpv(beGrp)*DimGrpv(gaGrp)
            if (Xyes(beGrp,gaGrp).eq.1) then
              call GetX (wrk(PossX),dim1,LunAux,LunName,1,1)
              Xyes(beGrp,gaGrp)=2
            else
              call mv0zero (dim1,dim1,wrk(PossX))
            end if
            call mv0zero (dim1,dim1,wrk(PossY))
c
cG               Extract H2(ga',I) <- T1(ga,I)
            call ExtT1 (wrk(PossH2),wrk(PossT1o),dimga,addga)
c
c
            adda=0
            do aGrp=1,NvGrp
            dima=DimGrpv(aGrp)
c
c           Test, if this Be' A' combination is to be run on this node
            if (BeAID(myRank,beGrp,aGrp).eq.0) goto 12
            if (printkey.ge.10) then
            write (6,*) ' o3v3 T2 - ID,be,a,ga',myRank,beGrp,aGrp,gaGrp
            end if
c
c
cXY1.1.1      read V1(ga',a',p,q) <- t2(ga',a',p,q)
              LunName=T2Name(gaGrp,aGrp)
              if (aGrp.eq.gaGrp) then
                dim1=dimga*(dimga+1)*no*no/2
                call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
cXY1.1.2        Expand V1(be',b',p,q) <- V2(be'b',p,q)
                dim1=dima*(dima+1)/2
                call ExpT2 (wrk(possV2),wrk(PossV1),dimga,dim1,no)
              else
                dim1=dimga*dima*no*no
                call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
              end if
c
              if (gaGrp.eq.beGrp) then
c@x              tento if nebol prevereny dostatocne
c              term T15 only for ga'=be'
c              use     : V1(be',a',I,J) = T2(be',a',I,J)
c                        H1(be',I) = T1(be,I) also H2
c              destroy : V2,V3 (also V4 is free)
c
cT15.1                Extract V2(a',I) <- T1(a,i)
                call ExtT1 (wrk(PossV2),wrk(PossT1o),dima,adda)
c
cT15.2                Make V3(be',u,a',i) <- 2t2(be,a,u,i)-t2(be,a,i,u)+t(a,u).t(be,i)
                call MkT_T15 (wrk(PossV3),wrk(PossV1),wrk(PossH1),
     c                        wrk(possV2),dimbe,dima,no)
c
cT15.3                Extract V2(a',i) <- Hvo(a',i) rutina identicka s ExtT1
                  call ExtT1 (wrk(PossV2),wrk(PossHvo),dima,adda)
c
cT15.4f         Calc H4(be',u) <<- V3(be',u,a',i) . V2(a',i)
                call mv0v1a3u (dimbe*no,dima*no,dima*no,dimbe*no,
     c                         dimbe*no,dima*no,1,1,
     c                         wrk(PossV3),wrk(PossV2),wrk(PossH4))
c
              end if
c
cX2.1              Map V3(a',u,ga',v) <- V1(ga',a',v,u)
               call Map4_3142 (wrk(PossV1),wrk(PossV3),dimga,dima,no,no)
c
cX2.2              Extract H3(a',be') <- Gvv(be,a)
              call ExH_X2 (wrk(PossGvv),wrk(PossH3),
     c                     dima,dimbe,nv,adda,addbe)
c
cX2.3              Calc V4(be',u,ga',v) <- H3(T)(a',be') . V3(a',u,ga',v)
              dim1=dimbe*dimga*no*no
                 call mv0zero (dim1,dim1,wrk(PossV4))
              dim1=dimga*no*no
              call mc0c1at3b (dima,dimbe,dima,dim1,dimbe,dim1,
     c                        dimbe,dima,dim1,
     c                        wrk(PossH3),wrk(PossV3),wrk(PossV4))
c
cX2.4f              Add X(be',u,ga',v) <<- +V4(be',u,ga',v)
cc            pozor na faktor, cele X je s vahou 0.5, sem teda asi 2
              dim1=dimbe*dimga*no*no
                 call mv0v1u (dim1,wrk(PossV4),1,wrk(PossX),1,2.0d0)
c
               if (aGrp.eq.beGrp) then
c                now a=be, i.e in V1 there is t2(ga',be',I,J)
c
cX3.1                Calc V2(ga',be',v,u) <- V1(ga',be',v,i) . Goo(i,u)
                dim1=dimga*dimbe*no*no
                call mv0zero (dim1,dim1,wrk(PossV2))
                dim1=dimga*dimbe*no
                call mc0c1a3b (dim1,no,no,no,dim1,no,
     c                         dim1,no,no,
     c                         wrk(PossV1),wrk(PossGoo),wrk(PossV2))
c
cX3.2                map V3((be',u,ga',v) <- V2(ga',be',v,u)
                call Map4_3142 (wrk(possV2),wrk(PossV3),
     c                          dimga,dimbe,no,no)
c
cX3.3f                Add X(be',u,ga',v) <<- -V3(be',u,ga',v)
cc              pozor na faktor, cele X je s vahou 0.5, sem teda asi -2
                dim1=dimbe*dimga*no*no
                call mv0v1u (dim1,wrk(PossV3),1,wrk(PossX),1,-2.0d0)
c
c
cX4.1                Map V2(be',ga',J,I) <- V1(ga',be',I,J)
                 Call Map4_2143 (wrk(PossV1),wrk(PossV2),
     c                          dimga,dimbe,no,no)
c
cX4.2                Make Tau in V2(be',ga',I,J)
                call MkTau_chcc (wrk(PossV2),wrk(PossH1),wrk(PossH2),
     c                           dimbe,dimga,no,1.0d0,1.0d0)
c
cX4.3                Extract V3(be',ga',ij) <- V2(be',ga',I,J)
                call ExV_X41 (wrk(PossV3),wrk(PossV2),dimbe*dimga,no)
c
cX4.5                V4(be',ga',u,v) <- V3(be',ga',ij) . A(ij,u,v)
                dim1=dimga*dimbe*no*no
                call mv0zero (dim1,dim1,wrk(PossV4))
                dim1=dimbe*dimga
                dim2=no*(no+1)/2
                 call mc0c1a3b (dim1,dim2,dim2,no*no,dim1,no*no,
     c                         dim1,dim2,no*no,
     c                       wrk(PossV3),wrk(PossA),wrk(PossV4))
c
cX4.6                Map V3(ij,V,U) <- A(ij,U,V)
                dim1=no*(no+1)/2
                call Map3_132 (wrk(PossA),wrk(PossV3),dim1,no,no)
c
cX4.7                Set A(ij,V,U) <- V3(ij,V,U)
                dim1=no*no*no*(no+1)/2
                call mv0u (dim1,wrk(PossV3),1,wrk(PossA),1)
c
cX4.8           Make V3(be',ga',ij) <- V2(be',ga',J,I) for i>=j
                call ExV_X43 (wrk(PossV3),wrk(PossV2),dimbe*dimga,no)
c
cX4.9                V4(be',ga',u,v) << - V3(be',ga',ij) . A(ij,u,v)
                dim1=dimbe*dimga
                dim2=no*(no+1)/2
                 call mc0c1a3b (dim1,dim2,dim2,no*no,dim1,no*no,
     c                         dim1,dim2,no*no,
     c                       wrk(PossV3),wrk(PossA),wrk(PossV4))
c
cX4.10          Map V3(ij,U,V) <- A(ij,V,U)
                dim1=no*(no+1)/2
                call Map3_132 (wrk(PossA),wrk(PossV3),dim1,no,no)
c
cX4.11          Set A(ij,U,V) <- V3(ij,U,V)
                dim1=no*no*no*(no+1)/2
                call mv0u (dim1,wrk(PossV3),1,wrk(PossA),1)
c
cX4.12          Extract H3(i,u,v) <- A(ii,u,v)
                call ExA_X4 (wrk(PossA),wrk(PossH3),no)
c
c
cX4.13          Extract V3(be',ga',i) <- V2(be',ga',I,J)
                call ExV_X42 (wrk(PossV3),wrk(PossV2),dimbe*dimga,no)
c
cX4.14          V4(be',ga',u,v) <<- - V3(be',ga',i) . H3(i,u,v)
                dim1=dimbe*dimga
                 call mc0c2a3b (dim1,no,no,no*no,dim1,no*no,
     c                         dim1,no,no*no,
     c                         wrk(PossV3),wrk(PossH3),wrk(PossV4))
c
cX4.15          Map V3(be',u,ga',v) <- V4(be',ga',u,v)
                 call Map4_1324 (wrk(PossV4),wrk(PossV3),
     c                          dimbe,dimga,no,no)
c
cX4.16f         Add X(be',u,ga',v) <<- V3(be',u,ga',v)
cc              pozor na faktor, cele X je s vahou 0.5, sem teda asi 1
                dim1=dimbe*dimga*no*no
                   call mv0v1u (dim1,wrk(PossV3),1,wrk(PossX),1,1.0d0)
c
cXe                Term Xe (V2,V3 and V4 can be used)
c                only in cholesky based approach
c
                if (intkey.eq.0) then
cXe.1                Ext V2(i,u,v,j) <- Aex(ij,uv)
                Call Ext_Aex (wrk(PossAex),wrk(PossV2),no)
c
cXe.2           V3(be',u,v,j) <- H1(=T1)(be',i) . V2(i,u,v,j)
                dim1=dimbe*no*no*no
                call mv0zero (dim1,dim1,wrk(PossV3))
                dim1=no*no*no
                call mc0c1a3b (dimbe,no,no,dim1,dimbe,dim1,
     c                         dimbe,no,dim1,
     c                         wrk(PossH1),wrk(PossV2),wrk(PossV3))
c
cXe.3                Map V4(j,ga') <- H2(ga',j)
                call Map2_21 (wrk(PossH2),wrk(PossV4),dimga,no)
c
cXe.4           V2(be',u,v,ga') <- V3(be',u,v,j) . V4(j,ga')
                dim1=dimbe*dimga*no*no
                call mv0zero (dim1,dim1,wrk(PossV2))
                dim1=dimbe*no*no
                call mc0c1a3b (dim1,no,no,dimga,dim1,dimga,
     c                         dim1,no,dimga,
     c                         wrk(PossV3),wrk(PossV4),wrk(PossV2))
c
cXe.5           Map V3(be',u,ga',v) <- V2(be,u,v,ga')
                call Map3_132(wrk(PossV2),wrk(PossV3),dimbe*no,no,dimga)
c
cXe.6f          Add X(be',u,ga',v) <<- - V3(be',u,ga',v)
cc              pozor na faktor, cele X je s vahou 0.5, sem teda asi -1
                dim1=dimbe*dimga*no*no
                   call mv0v1u (dim1,wrk(PossV3),1,wrk(PossX),1,-1.0d0)
                end if
c
               end if
c
c
cX1.2         map V2(a',i,ga',v) <- V1(ga',a',i,v)
              call Map4_3124 (wrk(PossV1),wrk(PossV2),dimga,dima,no,no)
c
cX1.3         read V1(be',u,i,a') <- Q(be',u,i,a')
              LunName=Tmp1Name(beGrp,aGrp)
              dim1=dimbe*dima*no*no
              call GetX (wrk(PossV1),dim1,LunAux,LunName,1,0)
cX1.4              Map V3(be',u,a',i) <- V1(be',u,i,a')
              dim1=dimbe*no
              call Map3_132 (wrk(PossV1),wrk(PossV3),dim1,no,dima)
c
cY1.2         read V1(be',u,i,a') <- K(be',u,i,a')
              dim1=dimbe*dima*no*no
              call GetX (wrk(PossV1),dim1,LunAux,LunName,0,1)
cY1.3              Map V4(be',u,a',i) <- V1(be',u,i,a')
              dim1=dimbe*no
              call Map3_132 (wrk(PossV1),wrk(PossV4),dim1,no,dima)
c
cY1.4f        Y(be',u,ga',v) <<- V4(be',u,a',i) . V2(a',i,ga',v)
              call mc0c1a3b (dimbe*no,dima*no,dima*no,dimga*no,
     c                       dimbe*no,dimga*no,
     c                       dimbe*no,dima*no,dimga*no,
     c                       wrk(PossV4),wrk(PossV2),wrk(PossY))
c
cX1.5         Make E: V1(a',i,ga',v) = 2 V2(a',v,ga',i) - V2(a',i,ga',v)
              call MkE_Y3 (wrk(PossV1),wrk(PossV2),dima,dimga,no)
c
cX1.6f        X(be',u,ga',v) <<- V3(be',u,a',i) . V1(a',i,ga',v)
              call mc0c1a3b (dimbe*no,dima*no,dima*no,dimga*no,
     c                       dimbe*no,dimga*no,
     c                       dimbe*no,dima*no,dimga*no,
     c                       wrk(PossV3),wrk(PossV1),wrk(PossX))
c
12             adda=adda+dima
            end do
c
cX0.1          read V1(be',u,ga',j) <- (be',u|ga',v)
c##          this contribution need to be taken only once for given be'ga'
c         thus will be executed on that node, where be',a'=1 is
c          calculated
          if (BeAID(myRank,beGrp,1).eq.1) then
            LunName=I2Name(beGrp,gaGrp)
            dim1=dimbe*dimga*no*no
            call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
          else
            dim1=dimbe*dimga*no*no
            call mv0zero (dim1,dim1,wrk(PossV1))
          end if
c
cX0.2f          Add X(be',u,ga',v) <<- 1/2 V1(be',u,ga',v)
cc            pozor na faktor, cele X je s vahou 0.5, sem teda asi 1
          dim1=dimbe*dimga*no*no
          call mv0v1u (dim1,wrk(PossV1),1,wrk(PossX),1,1.0d0)
c
c         save X,Y
          LunName=Tmp2Name(beGrp,gaGrp)
          dim1=dimbe*dimga*no*no
          XYyes(beGrp,gaGrp)=1
          call SaveX (wrk(PossX),dim1,LunAux,LunName,1,0)
          call SaveX (wrk(PossY),dim1,LunAux,LunName,0,1)
c
c@@
c        call Chck_Y (wrk(PossY),dimbe,addbe,dimga,addga)
c        call Chck_X (wrk(PossX),dimbe,addbe,dimga,addga)
c@@
          addga=addga+dimga
          end do
c
c
c##       These contributions must be taken only once per Beta'
          if (BetaID(myRank,beGrp).eq.1) then
c
cT12.1      calc V1(i,u) <<- Fvo(T)(a,i) . T1o(a,u)
            dim1=no*no
            call mv0zero (dim1,dim1,wrk(PossV1))
            call mc0c1at3b (nv,no,nv,no,no,no,
     c                   no,nv,no,
     c                   wrk(PossFvo),wrk(PossT1o),wrk(PossV1))
c
cT12.2      calc V2(be',u) <<- T1o(be',i) . V1(i,u)
            dim1=no*dimbe
            call mv0zero (dim1,dim1,wrk(PossV2))
            call mc0c1a3b (dimbe,no,no,no,dimbe,no,
     c                   dimbe,no,no,
     c                   wrk(PossT1o),wrk(PossV1),wrk(PossV2))
c
cT12.3f     Add H4(be',u) <<- -2.0 V2(be',u)
            call mv0v1u (dim1,wrk(PossV2),1,wrk(PossH4),1,-2.0d0)
c
c
cT13.1            Extract V1(a,be') <- Hvv(a,be)
              call ExH_T13 (wrk(PossV1),wrk(PossHvv),dimbe,addbe,nv)
c
cT13.2f     Add H4(be',u) <<- V1(T)(a,be') . T1(a,i)
            call mc0c1at3b (nv,dimbe,nv,no,dimbe,no,
     c                     dimbe,nv,no,
     c                     wrk(PossV1),wrk(PossT1o),wrk(PossH4))
c
c
cT14.1      calc V1(be',u) <- H1(be',i) . Hoo(i,u)
            dim1=no*dimbe
            call mv0zero (dim1,dim1,wrk(PossV1))
            call mc0c1a3b (dimbe,no,no,no,dimbe,no,
     c                     dimbe,no,no,
     c                     wrk(PossH1),wrk(PossHoo),wrk(PossV1))
c
cT14.2f     Add H4(be',u) <<- - V1(be',u)
            dim1=no*dimbe
            call mv0v1u (dim1,wrk(PossV1),1,wrk(PossH4),1,-1.0d0)
c
          end if
c
c
cT1G          Add T1n(be,u) <<- H4(be',u)
              call AdT_T17 (wrk(PossT1n),wrk(PossH4),
     c                  nv,dimbe,no,addbe,1.0d0)
c
c
11        addbe=addbe+dimbe
        end do
c
c
c##        correct Xyes (those, modified to 2 reset to 1)
        do dim1=1,NvGrp
        do dim2=1,NvGrp
          if (Xyes(dim1,dim2).eq.2) then
            Xyes(dim1,dim2)=1
          end if
        end do
        end do
c
c
#ifdef _MOLCAS_MPP_
c##        Synchronizacny bod:
c1        Allreduce T1n
c2        Allreduce X,Y (yet not needed)
        dim1=nv*no
        call gadgop (wrk(PossT1n),dim1,'+')
#endif
c
c
cT11.1        Add T1n(be,u) <<- Fvo(be,u)
        dim1=nv*no
           call mv0v1u (dim1,wrk(PossFvo),1,wrk(PossT1n),1,1.0d0)
c
c@@
c        call Chck_T1 (wrk(PossT1n),0)
c@@
c
c
        return
        end
c
c        ------------------------------
c
        subroutine Chck_T1 (T1,key)
c
c        check T1
c
        implicit none
#include "chcc1.fh"
        real*8 T1(1:nv,1:no)
c
c        help var
        integer be,u,i,j,a,b,bad
        integer key
        real*8 s
c
        bad=0
c
        do u=1,no
        do be=1,nv
c
          s=0.0d0
c
          do a=1,nv
          s=s+Hvvc(be,a)*T1c(a,u)
          end do
c
          do i=1,no
          s=s-Hooc(i,u)*T1c(be,i)
          end do
c
          do i=1,no
          do a=1,nv
             s=s+Hvoc(a,i)*
     c      (2.0d0*T2c(a,be,i,u)-T2c(a,be,u,i)+T1c(a,u)*T1c(be,i))
          end do
          end do
c
          do i=1,no
          do a=1,nv
           s=s+(2.0d0*Q21(a,i,be,u)-Q22(a,be,i,u))*T1c(a,i)
          end do
          end do
c
          do i=1,no
          do b=1,nv
          do a=1,nv
             s=s+(2.0d0*Q3(b,be,a,i)-Q3(a,be,b,i))*
     c          (T2c(a,b,i,u)+T1c(a,i)*T1c(b,u))
          end do
          end do
          end do
c
          do i=1,no
          do j=1,no
          do a=1,nv
             s=s-(2.0d0*Q1(a,i,u,j)-Q1(a,j,u,i))*
     c          (T2c(a,be,i,j)+T1c(a,i)*T1c(be,j))
          end do
          end do
          end do
c
c            s=s/(Oeo(u)-Oev(be))
c
          if (abs(T1(be,u)-s).gt.1.0d-10) then
            bad=bad+1
            if (key.eq.1) then
              T1(be,u)=s
            end if
          end if
c
        end do
        end do
c
         write (6,*) ' T1 test :',bad
c
        return
        end
c
