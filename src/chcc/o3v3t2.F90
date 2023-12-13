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

subroutine o3v3t2(wrk,wrksize,NvGrp,maxdim,LunAux)
! This routine does:
!
! --- Part III, generation of X, Y intermediates
!
!  T2(be,ga,u,v) <-
!  T21  +             [ (be,u|ga,v)                 ]
!  T22  + sum(i,j)    [ Ta(be,ga,ij) . A(i,j,u,v)   ]
!  T21 implemented as X0:
!  T22 implemented as X4:
!
! X(be,u,ga,v) <<-
!    X0 (T21)     1/2 (be,u|ga,v)
!    X1         + Q(be,u,i,a) . E(ga,a,i,v)
!    X2 (T24)   + Gvv(a,be)   . t2(ga,a,v,u)
!    X3 (T25)   - Goo(i,u)    . t2(ga,be,v,i)
!    X4 (T22)   + sum(i>=j)   [ Ta(be,ga,ij) . A(ij,u,v)  ]
!               - 1/2 sum(i)  [ Ta(be,ga,ii) . A(ii,u,v)  ]
!    Xe (T2ex)  - sum(i,j)    [ Aex(i,j,u,v) . T1(be,i) .T1(ga,j) ]
!                 Xe - o4v2 term is factorized to o4v + o3v2 terms
!
! Y(be,u,ga,v) = K(be,u,i,a) . t2(i,a,ga,v)
!
! where E(ga,a,i,v) = 2 t2(ga,a,v,i) - t2(ga,a,i,v)
!
! T1(be,u) << -
! T11  +            [  F(be,u)                   ]
! T12  - sum(a,i)   [ 2F(a,i) . t(be,i) . t(a,u) ]
! T13  + sum(a)     [  Hvv(a,be) . t(a,u)        ]
! T14  - sum(i)     [  Hoo(i,u)  . t(be,i)       ]
! T15  + sum(a,i)   [  Hoo(a,i)  .
!     (2t2(be,a,u,i) - t2(be,a,i,u) + t(a,u).t(be,i)) ]
!
! Extra cost - multiple reading:
!              1) N_be . T;
!              2) N_ga . Q
!              3) N_ga . K
!
!1 intermediates Q and K will be temporarily stored in files
!  QFil, KFil as follows
!  do beGrp
!    do aGrp
!      Q(K) (be',u,i,a')
!    end do
!  end do
!
!2 Structure of files, where selected group of (pq|rs) are
!  stored (V'O|OO) - I1 ; (V'O|V'O) - I2 ; (V'V'|OO) - I3
!
!  (A'I|JK)  I1inxx xx - Group of A'
!
!  (A'I|B'J) I2xxyy xx - Group of A'
!                   yy - Group of B'
!
!  (A'B'|IJ) I3xxyy xx - Group of A'
!                   yy - Group of B'
!
!3 Structure of Cholesky vector files
!
!  L1(m,I ,A')  L1vcxx xx - Group of A'
!
!  L2(m,A',B')  L2xxyy xx - Group of A'
!                      yy - Group of B'
!
!4 Structure of Amplitude file
!  t2(A',B',IJ)  T2xxyy xx - Group of A'
!                       yy - Group of B'

use Index_Functions, only: nTri_Elem
use chcc_global, only: BeAID, BetaID, DimGrpv, I2Name, intkey, no, nv, PosA, PosAex, PosFree, PosFvo, PosGoo, PosGvv, PosHoo, &
                       PosHvo, PosHvv, PosT1n, PosT1o, printkey, T2Name, TCpu, TCpu_l, Tmp1Name, Tmp2Name, TWall, TWall_l, TWall0, &
                       Xyes, XYyes
use Para_Info, only: MyRank
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: wrksize, NvGrp, maxdim, LunAux
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp) :: adda, addbe, addga, aGrp, beGrp, dim_1, dim_2, dima, dimbe, dimga, gaGrp, PosH1, PosH2, PosH3, PosH4, PosK, &
                     PosQ, PosT, PosV1, PosV2, PosV3, PosV4, PosX, PosY
character(len=6) :: LunName

! --- introduction part ---
!
! Distribute memory
! @@ yet distribution is the same for all three o3v3 drivers,
! v buducnosti nosnost usposobit vsetky 3 vlastne

PosT = PosFree
call DistMemo3v3t2(maxdim,PosV1,PosV2,PosV3,PosV4,PosH1,PosH2,PosH3,PosH4,PosK,PosQ,PosT)

if (printkey >= 10) write(u6,*) ' Last Value :',PosT,wrksize
if (PosT > wrksize) then
  !mp write(u6,*) ' Nieje dobre, Dr. Ch.  Kokotopuloss',
  write(u6,*) ' Not Enough memory in o3v3t2 step! Increase large and/or small segmentation ', &
              real(PosT,kind=wp)/real(wrksize,kind=wp)
  call abend()
end if

!G Vanish XYyes

XYyes(1:NvGrp,1:NvGrp) = 0

PosX = PosQ
PosY = PosK

addbe = 0
do beGrp=1,NvGrp
  dimbe = DimGrpv(beGrp)

  !## test, if this beGrp is planed to be run on this node
  dim_1 = sum(BeAID(myRank,beGrp,1:NvGrp))
  dim_1 = dim_1+BetaID(myRank,beGrp)
  if (dim_1 /= 0) then
    if (printkey > 1) write(u6,*) ' o3v3 T2 - ID, beGrp',myRank,beGrp ! nie som si isty ....
    !mp
    call CWTime(TCpu,TWall)
    if (printkey > 1) then
      write(u6,*)
      write(u6,'(A,f18.1)') ' Cpu last call [s] = ',TCpu-TCpu_l
      write(u6,'(A,f18.1)') 'Wall last call [s] = ',TWall-TWall_l
      write(u6,*)
      write(u6,'(A,f18.1)') 'Total Cpu  [s] = ',TCpu
      write(u6,'(A,f18.1)') 'Total Wall [s] = ',TWall-TWall0
      write(u6,'(A,f18.2)') 'TCpu/TWall [%] = ',100.0_wp*TCpu/(TWall-TWall0)
      write(u6,*)
    end if
    TCpu_l = TCpu
    TWall_l = TWall
    !mp

    !T1G vanish H4(be',u)
    dim_1 = dimbe*no
    wrk(PosH4:PosH4+dim_1-1) = Zero

    !G Extract H1(be',I) <- T1(be,I)
    call ExtT1(wrk(PosH1),wrk(PosT1o),dimbe,addbe)

    addga = 0
    do gaGrp=1,NvGrp
      dimga = DimGrpv(gaGrp)

      ! initialize X(be',u,ga',v) as X0, calc in step II
      ! vanish Y(be',u,ga',v)
      LunName = Tmp2Name(beGrp,gaGrp)
      dim_1 = no*no*DimGrpv(beGrp)*DimGrpv(gaGrp)
      if (Xyes(beGrp,gaGrp) == 1) then
        call GetX(wrk(PosX),dim_1,LunAux,LunName,1,1)
        Xyes(beGrp,gaGrp) = 2
      else
        wrk(PosX:PosX+dim_1-1) = Zero
      end if
      wrk(PosY:PosY+dim_1-1) = Zero

      !G Extract H2(ga',I) <- T1(ga,I)
      call ExtT1(wrk(PosH2),wrk(PosT1o),dimga,addga)

      adda = 0
      do aGrp=1,NvGrp
        dima = DimGrpv(aGrp)

        ! Test, if this Be' A' combination is to be run on this node
        if (BeAID(myRank,beGrp,aGrp) /= 0) then
           if (printkey >= 10) write(u6,*) ' o3v3 T2 - ID,be,a,ga',myRank,beGrp,aGrp,gaGrp

           !XY1.1.1 read V1(ga',a',p,q) <- t2(ga',a',p,q)
           LunName = T2Name(gaGrp,aGrp)
           if (aGrp == gaGrp) then
             dim_1 = nTri_Elem(dimga)*no*no
             call GetX(wrk(PosV2),dim_1,LunAux,LunName,1,1)
             !XY1.1.2 Expand V1(be',b',p,q) <- V2(be'b',p,q)
             dim_1 = nTri_Elem(dima)
             call ExpT2(wrk(posV2),wrk(PosV1),dimga,dim_1,no)
           else
             dim_1 = dimga*dima*no*no
             call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)
           end if

           if (gaGrp == beGrp) then
             !@ tento if nebol prevereny dostatocne
             ! term T15 only for ga'=be'
             ! use     : V1(be',a',I,J) = T2(be',a',I,J)
             !           H1(be',I) = T1(be,I) also H2
             ! destroy : V2,V3 (also V4 is free)

             !T15.1 Extract V2(a',I) <- T1(a,i)
             call ExtT1(wrk(PosV2),wrk(PosT1o),dima,adda)

             !T15.2 Make V3(be',u,a',i) <- 2t2(be,a,u,i)-t2(be,a,i,u)+t(a,u).t(be,i)
             call MkT_T15(wrk(PosV3),wrk(PosV1),wrk(PosH1),wrk(posV2),dimbe,dima,no)

             !T15.3 Extract V2(a',i) <- Hvo(a',i) rutina identicka s ExtT1
             call ExtT1(wrk(PosV2),wrk(PosHvo),dima,adda)

             !T15.4f Calc H4(be',u) <<- V3(be',u,a',i) . V2(a',i)
             call mv0v1a3u(dimbe*no,dima*no,dima*no,dimbe*no,dimbe*no,dima*no,1,1,wrk(PosV3),wrk(PosV2),wrk(PosH4))

           end if

           !X2.1 Map V3(a',u,ga',v) <- V1(ga',a',v,u)
           call Map4_3142(wrk(PosV1),wrk(PosV3),dimga,dima,no,no)

           !X2.2 Extract H3(a',be') <- Gvv(be,a)
           call ExH_X2(wrk(PosGvv),wrk(PosH3),dima,dimbe,nv,adda,addbe)

           !X2.3 Calc V4(be',u,ga',v) <- H3(T)(a',be') . V3(a',u,ga',v)
           dim_1 = dimbe*dimga*no*no
           wrk(PosV4:PosV4+dim_1-1) = Zero
           dim_1 = dimga*no*no
           call mc0c1at3b(dima,dimbe,dima,dim_1,dimbe,dim_1,dimbe,dima,dim_1,wrk(PosH3),wrk(PosV3),wrk(PosV4))

           !X2.4f Add X(be',u,ga',v) <<- +V4(be',u,ga',v)
           ! pozor na faktor, cele X je s vahou 0.5, sem teda asi 2
           dim_1 = dimbe*dimga*no*no
           wrk(PosX:PosX+dim_1-1) = wrk(PosX:PosX+dim_1-1)+Two*wrk(PosV4:PosV4+dim_1-1)

           if (aGrp == beGrp) then
             ! now a=be, i.e in V1 there is t2(ga',be',I,J)

             !X3.1 Calc V2(ga',be',v,u) <- V1(ga',be',v,i) . Goo(i,u)
             dim_1 = dimga*dimbe*no*no
             wrk(PosV2:PosV2+dim_1-1) = Zero
             dim_1 = dimga*dimbe*no
             call mc0c1a3b(dim_1,no,no,no,dim_1,no,dim_1,no,no,wrk(PosV1),wrk(PosGoo),wrk(PosV2))

             !X3.2 map V3((be',u,ga',v) <- V2(ga',be',v,u)
             call Map4_3142(wrk(posV2),wrk(PosV3),dimga,dimbe,no,no)

             !X3.3f Add X(be',u,ga',v) <<- -V3(be',u,ga',v)
             ! pozor na faktor, cele X je s vahou 0.5, sem teda asi -2
             dim_1 = dimbe*dimga*no*no
             wrk(PosX:PosX+dim_1-1) = wrk(PosX:PosX+dim_1-1)-Two*wrk(PosV3:PosV3+dim_1-1)

             !X4.1 Map V2(be',ga',J,I) <- V1(ga',be',I,J)
             call Map4_2143(wrk(PosV1),wrk(PosV2),dimga,dimbe,no,no)

             !X4.2 Make Tau in V2(be',ga',I,J)
             call MkTau_chcc(wrk(PosV2),wrk(PosH1),wrk(PosH2),dimbe,dimga,no,One,One)

             !X4.3 Extract V3(be',ga',ij) <- V2(be',ga',I,J)
             call ExV_X41(wrk(PosV3),wrk(PosV2),dimbe*dimga,no)

             !X4.5 V4(be',ga',u,v) <- V3(be',ga',ij) . A(ij,u,v)
             dim_1 = dimga*dimbe*no*no
             wrk(PosV4:PosV4+dim_1-1) = Zero
             dim_1 = dimbe*dimga
             dim_2 = nTri_Elem(no)
             call mc0c1a3b(dim_1,dim_2,dim_2,no*no,dim_1,no*no,dim_1,dim_2,no*no,wrk(PosV3),wrk(PosA),wrk(PosV4))

             !X4.6 Map V3(ij,V,U) <- A(ij,U,V)
             dim_1 = nTri_Elem(no)
             call Map3_132(wrk(PosA),wrk(PosV3),dim_1,no,no)

             !X4.7 Set A(ij,V,U) <- V3(ij,V,U)
             dim_1 = no*no*nTri_Elem(no)
             wrk(PosA:PosA+dim_1-1) = wrk(PosV3:PosV3+dim_1-1)

             !X4.8 Make V3(be',ga',ij) <- V2(be',ga',J,I) for i>=j
             call ExV_X43(wrk(PosV3),wrk(PosV2),dimbe*dimga,no)

             !X4.9 V4(be',ga',u,v) << - V3(be',ga',ij) . A(ij,u,v)
             dim_1 = dimbe*dimga
             dim_2 = nTri_Elem(no)
             call mc0c1a3b(dim_1,dim_2,dim_2,no*no,dim_1,no*no,dim_1,dim_2,no*no,wrk(PosV3),wrk(PosA),wrk(PosV4))

             !X4.10 Map V3(ij,U,V) <- A(ij,V,U)
             dim_1 = nTri_Elem(no)
             call Map3_132(wrk(PosA),wrk(PosV3),dim_1,no,no)

             !X4.11 Set A(ij,U,V) <- V3(ij,U,V)
             dim_1 = no*no*nTri_Elem(no)
             wrk(PosA:PosA+dim_1-1) = wrk(PosV3:PosV3+dim_1-1)

             !X4.12 Extract H3(i,u,v) <- A(ii,u,v)
             call ExA_X4(wrk(PosA),wrk(PosH3),no)

             !X4.13 Extract V3(be',ga',i) <- V2(be',ga',I,J)
             call ExV_X42(wrk(PosV3),wrk(PosV2),dimbe*dimga,no)

             !X4.14 V4(be',ga',u,v) <<- - V3(be',ga',i) . H3(i,u,v)
             dim_1 = dimbe*dimga
             call mc0c2a3b(dim_1,no,no,no*no,dim_1,no*no,dim_1,no,no*no,wrk(PosV3),wrk(PosH3),wrk(PosV4))

             !X4.15 Map V3(be',u,ga',v) <- V4(be',ga',u,v)
             call Map4_1324(wrk(PosV4),wrk(PosV3),dimbe,dimga,no,no)

             !X4.16f Add X(be',u,ga',v) <<- V3(be',u,ga',v)
             ! pozor na faktor, cele X je s vahou 0.5, sem teda asi 1
             dim_1 = dimbe*dimga*no*no
             wrk(PosX:PosX+dim_1-1) = wrk(PosX:PosX+dim_1-1)+wrk(PosV3:PosV3+dim_1-1)

             !Xe Term Xe (V2,V3 and V4 can be used)
             ! only in cholesky based approach

             if (intkey == 0) then
               !Xe.1 Ext V2(i,u,v,j) <- Aex(ij,uv)
               call Ext_Aex(wrk(PosAex),wrk(PosV2),no)

               !Xe.2 V3(be',u,v,j) <- H1(=T1)(be',i) . V2(i,u,v,j)
               dim_1 = dimbe*no*no*no
               wrk(PosV3:PosV3+dim_1-1) = Zero
               dim_1 = no*no*no
               call mc0c1a3b(dimbe,no,no,dim_1,dimbe,dim_1,dimbe,no,dim_1,wrk(PosH1),wrk(PosV2),wrk(PosV3))

               !Xe.3 Map V4(j,ga') <- H2(ga',j)
               call Map2_21(wrk(PosH2),wrk(PosV4),dimga,no)

               !Xe.4 V2(be',u,v,ga') <- V3(be',u,v,j) . V4(j,ga')
               dim_1 = dimbe*dimga*no*no
               wrk(PosV2:PosV2+dim_1-1) = Zero
               dim_1 = dimbe*no*no
               call mc0c1a3b(dim_1,no,no,dimga,dim_1,dimga,dim_1,no,dimga,wrk(PosV3),wrk(PosV4),wrk(PosV2))

               !Xe.5 Map V3(be',u,ga',v) <- V2(be,u,v,ga')
               call Map3_132(wrk(PosV2),wrk(PosV3),dimbe*no,no,dimga)

               !Xe.6f Add X(be',u,ga',v) <<- - V3(be',u,ga',v)
               ! pozor na faktor, cele X je s vahou 0.5, sem teda asi -1
               dim_1 = dimbe*dimga*no*no
               wrk(PosX:PosX+dim_1-1) = wrk(PosX:PosX+dim_1-1)-wrk(PosV3:PosV3+dim_1-1)
             end if

           end if

           !X1.2 map V2(a',i,ga',v) <- V1(ga',a',i,v)
           call Map4_3124(wrk(PosV1),wrk(PosV2),dimga,dima,no,no)

           !X1.3 read V1(be',u,i,a') <- Q(be',u,i,a')
           LunName = Tmp1Name(beGrp,aGrp)
           dim_1 = dimbe*dima*no*no
           call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,0)
           !X1.4 Map V3(be',u,a',i) <- V1(be',u,i,a')
           dim_1 = dimbe*no
           call Map3_132(wrk(PosV1),wrk(PosV3),dim_1,no,dima)

           !Y1.2 read V1(be',u,i,a') <- K(be',u,i,a')
           dim_1 = dimbe*dima*no*no
           call GetX(wrk(PosV1),dim_1,LunAux,LunName,0,1)
           !Y1.3 Map V4(be',u,a',i) <- V1(be',u,i,a')
           dim_1 = dimbe*no
           call Map3_132(wrk(PosV1),wrk(PosV4),dim_1,no,dima)

           !Y1.4f Y(be',u,ga',v) <<- V4(be',u,a',i) . V2(a',i,ga',v)
           call mc0c1a3b(dimbe*no,dima*no,dima*no,dimga*no,dimbe*no,dimga*no,dimbe*no,dima*no,dimga*no,wrk(PosV4),wrk(PosV2), &
                         wrk(PosY))

           !X1.5 Make E: V1(a',i,ga',v) = 2 V2(a',v,ga',i) - V2(a',i,ga',v)
           call MkE_Y3(wrk(PosV1),wrk(PosV2),dima,dimga,no)

           !X1.6f X(be',u,ga',v) <<- V3(be',u,a',i) . V1(a',i,ga',v)
           call mc0c1a3b(dimbe*no,dima*no,dima*no,dimga*no,dimbe*no,dimga*no,dimbe*no,dima*no,dimga*no,wrk(PosV3),wrk(PosV1), &
                         wrk(PosX))

        end if
        adda = adda+dima
      end do

      !X0.1 read V1(be',u,ga',j) <- (be',u|ga',v)
      !## this contribution need to be taken only once for given be'ga'
      !   thus will be executed on that node, where be',a'=1 is calculated
      if (BeAID(myRank,beGrp,1) == 1) then
        LunName = I2Name(beGrp,gaGrp)
        dim_1 = dimbe*dimga*no*no
        call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)
      else
        dim_1 = dimbe*dimga*no*no
        wrk(PosV1:PosV1+dim_1-1) = Zero
      end if

      !X0.2f Add X(be',u,ga',v) <<- 1/2 V1(be',u,ga',v)
      ! pozor na faktor, cele X je s vahou 0.5, sem teda asi 1
      dim_1 = dimbe*dimga*no*no
      wrk(PosX:PosX+dim_1-1) = wrk(PosX:PosX+dim_1-1)+wrk(PosV1:PosV1+dim_1-1)

      ! save X,Y
      LunName = Tmp2Name(beGrp,gaGrp)
      dim_1 = dimbe*dimga*no*no
      XYyes(beGrp,gaGrp) = 1
      call SaveX(wrk(PosX),dim_1,LunAux,LunName,1,0)
      call SaveX(wrk(PosY),dim_1,LunAux,LunName,0,1)

      !@@
      !call Chck_Y(wrk(PosY),dimbe,addbe,dimga,addga)
      !call Chck_X(wrk(PosX),dimbe,addbe,dimga,addga)
      !@@
      addga = addga+dimga
    end do

    !## These contributions must be taken only once per Beta'
    if (BetaID(myRank,beGrp) == 1) then

      !T12.1 calc V1(i,u) <<- Fvo(T)(a,i) . T1o(a,u)
      dim_1 = no*no
      wrk(PosV1:PosV1+dim_1-1) = Zero
      call mc0c1at3b(nv,no,nv,no,no,no,no,nv,no,wrk(PosFvo),wrk(PosT1o),wrk(PosV1))

      !T12.2 calc V2(be',u) <<- T1o(be',i) . V1(i,u)
      dim_1 = no*dimbe
      wrk(PosV2:PosV2+dim_1-1) = Zero
      call mc0c1a3b(dimbe,no,no,no,dimbe,no,dimbe,no,no,wrk(PosT1o),wrk(PosV1),wrk(PosV2))

      !T12.3f Add H4(be',u) <<- -2.0 V2(be',u)
      wrk(PosX:PosX+dim_1-1) = wrk(PosX:PosX+dim_1-1)-Two*wrk(PosV2:PosV2+dim_1-1)

      !T13.1 Extract V1(a,be') <- Hvv(a,be)
      call ExH_T13(wrk(PosV1),wrk(PosHvv),dimbe,addbe,nv)

      !T13.2f Add H4(be',u) <<- V1(T)(a,be') . T1(a,i)
      call mc0c1at3b(nv,dimbe,nv,no,dimbe,no,dimbe,nv,no,wrk(PosV1),wrk(PosT1o),wrk(PosH4))

      !T14.1 calc V1(be',u) <- H1(be',i) . Hoo(i,u)
      dim_1 = no*dimbe
      wrk(PosV1:PosV1+dim_1-1) = Zero
      call mc0c1a3b(dimbe,no,no,no,dimbe,no,dimbe,no,no,wrk(PosH1),wrk(PosHoo),wrk(PosV1))

      !T14.2f Add H4(be',u) <<- - V1(be',u)
      dim_1 = no*dimbe
      wrk(PosH4:PosH4+dim_1-1) = wrk(PosH4:PosH4+dim_1-1)-wrk(PosV1:PosV1+dim_1-1)

    end if

    !T1G Add T1n(be,u) <<- H4(be',u)
    call AdT_T17(wrk(PosT1n),wrk(PosH4),nv,dimbe,no,addbe,One)

  end if
  addbe = addbe+dimbe
end do

!## correct Xyes (those, modified to 2 reset to 1)
do dim_1=1,NvGrp
  do dim_2=1,NvGrp
    if (Xyes(dim_1,dim_2) == 2) Xyes(dim_1,dim_2) = 1
  end do
end do

#ifdef _MOLCAS_MPP_
!## Synchronizacny bod:
!1 Allreduce T1n
!2 Allreduce X,Y (yet not needed)
dim_1 = nv*no
call gadgop(wrk(PosT1n),dim_1,'+')
#endif

!T11.1 Add T1n(be,u) <<- Fvo(be,u)
dim_1 = nv*no
wrk(PosT1n:PosT1n+dim_1-1) = wrk(PosT1n:PosT1n+dim_1-1)+wrk(PosFvo:PosFvo+dim_1-1)

!@@
!call Chck_T1(wrk(PosT1n),0)
!@@

return

end subroutine o3v3t2
