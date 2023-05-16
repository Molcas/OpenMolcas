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

subroutine o3v3chol(wrk,wrksize,NvGrp,maxdim,LunAux)
! This routine does:
!
! --- Part II - generation of Q3,K3,Gvv,Goo,T26-9 in X form
!
! T(be,u)      <<-
! T162 - sum(b,i)   [ (iu|be,b) . t(b,i) ]
! T17  + sum(a,b,i) [ 2(a,i|b,be) - (b,i|a,be) ] . Ta(a,b,i,u)
!    = + sum(a,b,i)  (a,i|b,be) . [ 2 Ta(a,b,i,u) - Ta(a,b,u,i) ]
!
! N.B. T17 implementede in expression:
! T(a,u)      <<-
! T17  + sum(be,b,i)  (b,i|be,a) . [ 2 Ta(b,be,i,u) - Ta(b,be,u,i) ]
!
! Q(be,u,i,a) <-
!        Q3     +   [ 2(b,be|a,i) - (b,i|a,be) ] . t(b,u )
!
! K(be,u,i,a) <-
!        K3     +   [  (b,i |a,be)             ] . t(b ,u)
!
! Gvv(be,a)    <-
! Gvv1  +           [ Hvv(a,be) ]
! Gvv2  - sum(i)    [ Fvo(a,i) . t1(be,i) ]
! Gvv3  + sum(b,i)  [ -(b,be|a,i) +2(b,i|a,be) ] . t(b,i )
!
! Goo(i,u) <-
! Goo1  +  <-       [ Hoo(i,u) ]
! Goo2  + sum(a)    [ Fvo(a,i) . t1(a,u) ]
! Goo3  + sum(a,j)  [ (2(aj|iu) - (ai|ju)) . t1(a,j) ]
!
! terms Goo2,3 implemented as:
! Goo2  + sum(be)   [ Fvo(be,i) . t1(be,u) ]
! Goo3  + sum(be,j) [ (2(be,j|iu) - (be,i|ju)) . t1(be,j) ]
!
! T2(be,ga,u,v) <<-
! T26  + sum(a)     [ (be,u|ga,a) . t1(a,v)         ]
! T27  - sum(a,i)   [ (iu|ga,a) . t1(be,i) . t1(a,v)]
! T28  - sum(i)     [ (be,u|iv)          . t1(ga,i) ]
! T29  - sum(a,i)   [ (be,u|ia). t1(a,v) . t1(ga,i) ]
!
! N.B. T26-T29 implementede in expression:
! T2(a,be,u,v)  <-
! T26  + sum(b)     [           (u,a|be,b) . t1(b,v)  ]
! T27  - sum(b,j)   [ t1(a,j) . (u,j|be,b) . t1(b,v)  ]
! T28  - sum(j)     [           (u,a|j,v)  . t1(be,j) ]
! T29  - sum(b,j)   [ t1(b,v) . (u,a|j,b)  . t1(be,j) ]
! T2(a,be,u,v) is stored as X(a,u,be,v) in Tmp2Name(a,be)
!
! N.B. Process V3(j,v,u,a') <- V1(m,j,v) . L1(T)(m,u,a') run
! for Na*Nbe times, ie. Nbe . o3mv competee with Nbe . I+O (o3v)
!
! Extra cost - multiple reading:
!              1) N_be . T;
!              2) N_ga . Q
!              3) N_ga . K
!
!1 intermediates Q and K will be temporarily stored in files
!  Tmp1Name(be,a)
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
use chcc_global, only: BeAID, BetaID, DimGrpv, I1Name, L0Name, L1Name, L2Name, nc, no, nv, PosFree, PosFvo, PosGoo, PosGvv, &
                       PosHoo, PosHvv, PosT1n, PosT1o, printkey, T2Name, TCpu, TCpu_l, Tmp1Name, Tmp2Name, TWall, TWall_l, TWall0, &
                       Xyes
use Para_Info, only: MyRank
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: wrksize, NvGrp, maxdim, LunAux
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp) :: adda, addb, addbe, aGrp, beGrp, bGrp, dim_1, dim_2, dima, dimb, dimbe, PosH1, PosH2, PosH3, PosH4, PosK, &
                     PosM1, PosM2, PosM3, PosM4, PosM5, PosQ, PosT, PosV1, PosV2, PosV3, PosV4
character(len=6) :: LunName

! --- introduction part ---

! Distribute memory
! @@ vyladit distmem, mozno netreba vsetko, toto je len prefackany
!   DistMemo3v3

PosT = PosFree
if (printkey >= 10) write(u6,*) 'PosFree',PosFree,NvGrp,maxdim
call DistMemo3v3chol(maxdim,PosV1,PosV2,PosV3,PosV4,PosH1,PosH2,PosH3,PosH4,PosM1,PosM2,PosM3,PosM4,PosM5,PosK,PosQ,PosT)

if (printkey >= 10) write(u6,*) ' Last Value :',PosT,wrksize
if (PosT > wrksize) then
  !mp write(u6,*) ' Nieje dobre, Dr. Ch.  Kokotopuloss',
  write(u6,*) ' Not Enough memory in o3v3chol step! Increase large and/or small segmentation ', &
              real(PosT,kind=wp)/real(wrksize,kind=wp)
  call abend()
end if

!@@
!call Chck_mkJ()
!call Chck_mkK()
!@@

! read V1(m,ij) <- L1(m,ij)
LunName = L0Name
dim_1 = nc*nTri_Elem(no)
call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)

! Expand M4(m,i,j) <- V1(m,ij)
dim_1 = nTri_Elem(no)
call Exp1(wrk(PosV1),wrk(PosM4),nc,dim_1,no)

!G Vanish Gvv,Goo
dim_1 = nv*nv
wrk(PosGvv:PosGvv+dim_1-1) = Zero
dim_1 = no*no
wrk(PosGoo:PosGoo+dim_1-1) = Zero

!##G Vanish Xyes

Xyes(1:NvGrp,1:NvGrp) = 0

addbe = 0

do beGrp=1,NvGrp
  dimbe = DimGrpv(beGrp)

  !## test, if this beGrp is planed to be run on this node
  dim_1 = sum(BeAID(myRank,beGrp,1:NvGrp))
  dim_1 = dim_1+BetaID(myRank,beGrp)
  if (dim_1 /= 0) then
    if (printkey > 1) then
      write(u6,*) ' Chol PID, beGrp',myRank,beGrp
      !mp
      call CWTime(TCpu,TWall)
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

    ! Extract H1(be',I) <- T1(be,I)
    call ExtT1(wrk(PosH1),wrk(PosT1o),dimbe,addbe)

    ! vanish M1(m,i,u),M2(m,be',u),M3(m,be',u)
    dim_1 = nc*no*no
    wrk(PosM1:PosM1+dim_1-1) = Zero
    dim_1 = nc*no*dimbe
    wrk(PosM2:PosM2+dim_1-1) = Zero
    wrk(PosM3:PosM3+dim_1-1) = Zero

    !T1G vanish H4(be',u)
    dim_1 = dimbe*no
    wrk(PosH4:PosH4+dim_1-1) = Zero

    if (BetaID(myRank,beGrp) == 1) then
      !## contributions Goo2,3  must be taken only once per Beta'

      !Goo2.1 Extract H2(be',i) <- Fvo(be,i) (vlastnu rutinu netreba)
      call ExtT1(wrk(PosH2),wrk(PosFvo),dimbe,addbe)

      !Goo2.2f Goo(i,u) <<- H2(T)(be',i) . H1(be',u)
      call mc0c1at3b(dimbe,no,dimbe,no,no,no,no,dimbe,no,wrk(PosH2),wrk(PosH1),wrk(PosGoo))

      !Goo3.1 read V1(be',I,JK) <- I1(be',I,JK)
      dim_1 = dimbe*no*nTri_Elem(no)
      LunName = I1Name(beGrp)
      call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)

      !Goo3.2 Make V2(be',j,i,u) <- 2 V1(be',j|iu) - V1(be',i|ju)
      call MkV_Goo3(wrk(PosV1),wrk(PosV2),dimbe,no)

      !Goo3.3f Goo(i,u) <<- V2(T)(be',j,i,u) . H1(be,j)
      dim_1 = dimbe*no
      !@x dorobit rutinu  mv0v1at3u
      !call mv0v1at3u (dim_1,no*no,dim_1,no*no,no*no,dim_1,1,1,wrk(PosV2),wrk(PosH1),wrk(PosGoo))

      ! zatial odflaknute takto
      call mc0c1at3b(dim_1,no*no,dim_1,1,no*no,1,no*no,dim_1,1,wrk(PosV2),wrk(PosH1),wrk(PosGoo))

      !@x
    end if

    addb = 0

    do bGrp=1,NvGrp
      dimb = DimGrpv(bGrp)

      !G Extract H2(b',I) <- T1(b,I)
      call ExtT1(wrk(PosH2),wrk(PosT1o),dimb,addb)

      !QK3.1 read V3(m,i,b') <- L1(m,i,b')
      LunName = L1Name(bGrp)
      dim_1 = nc*no*dimb
      call GetX(wrk(PosV3),dim_1,LunAux,LunName,1,1)

      !QK3.3 M1(m,i,u)   <<-  V3(m,i  ,b') . H2(b',u)
      call mc0c1a3b(nc*no,dimb,dimb,no,nc*no,no,nc*no,dimb,no,wrk(PosV3),wrk(PosH2),wrk(PosM1))

      !Q3.4 read V1(m,be',b')<- [V2(m,beb') or V2(m,b'be)]<- L2(m,be',b')
      !     in cases be <= b through interstep V2
      if (beGrp > bGrp) then
        LunName = L2Name(beGrp,bGrp)
        dim_1 = nc*dimb*dimbe
        call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)
      else if (beGrp == bGrp) then
        LunName = L2Name(beGrp,bGrp)
        dim_1 = nc*dimb*(dimbe+1)/2
        call GetX(wrk(PosV2),dim_1,LunAux,LunName,1,1)
        dim_1 = dimb*(dimbe+1)/2
        call Exp1(wrk(PosV2),wrk(PosV1),nc,dim_1,dimb)
      else
        LunName = L2Name(bGrp,beGrp)
        dim_1 = nc*dimb*dimbe
        call GetX(wrk(PosV2),dim_1,LunAux,LunName,1,1)
        call Map3_132(wrk(PosV2),wrk(PosV1),nc,dimb,dimbe)
      end if

      !Q3.5 M2(m,be',u) <<-  V1(m,be',b') . H2(b',u)
      call mc0c1a3b(nc*dimbe,dimb,dimb,no,nc*dimbe,no,nc*dimbe,dimb,no,wrk(PosV1),wrk(PosH2),wrk(PosM2))

      ! read V2(be'b',I,J) <- T2(be',b,I,J)
      LunName = T2Name(beGrp,bGrp)
      if (bGrp == beGrp) then
        dim_1 = nTri_Elem(dimbe)*no*no
        call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)
        ! Expand V2(be',b',I,J) <- V1(be'b',I,J)
        dim_1 = nTri_Elem(dimb)
        call ExpT2(wrk(posV1),wrk(PosV2),dimbe,dim_1,no)
      else
        dim_1 = dimbe*dimb*no*no
        call GetX(wrk(PosV2),dim_1,LunAux,LunName,1,1)
      end if

      ! Make Tau(be',b',I,J) in V2(be',b',I,J)
      call MkTau_chcc(wrk(PosV2),wrk(PosH1),wrk(PosH2),dimbe,dimb,no,One,One)

      ! Make V1(i,b',be',u) <- 2 V2(be',b',u,i) - V2(be',b',i,u)
      call MkT_T17(wrk(PosV1),wrk(PosV2),dimb,dimbe,no)

      ! M3(m,be',u) <<- V3(m,i,b') . V1(i,b',be',u)
      dim_1 = dimb*no
      dim_2 = dimbe*no
      call mc0c1a3b(nc,dim_1,dim_1,dim_2,nc,dim_2,nc,dim_1,dim_2,wrk(PosV3),wrk(PosV1),wrk(PosM3))

      addb = addb+dimb
    end do

    if (BetaID(myRank,beGrp) == 1) then
      !## contribution T16  must be taken only once per Beta'

      ! Map V1(be',m,i) <- M2(m,be',i)
      call Map3_213(wrk(PosM2),wrk(PosV1),nc,dimbe,no)

      !T162.f H4(be',u) <<- - V1(be',m,i) . M4(m,i,u)
      call mc0c2a3b(dimbe,nc*no,nc*no,no,dimbe,no,dimbe,nc*no,no,wrk(PosV1),wrk(PosM4),wrk(PosH4))

    end if

    !T27.1 V4(be',v,u,j ) <- M2(T)(m,be',v) . M4(m,u,j)
    dim_1 = dimbe*no*no*no
    wrk(PosV4:PosV4+dim_1-1) = Zero
    call mc0c1at3b(nc,dimbe*no,nc,no*no,dimbe*no,no*no,dimbe*no,nc,no*no,wrk(PosM2),wrk(PosM4),wrk(PosV4))

    adda = 0

    do aGrp=1,NvGrp
      dima = DimGrpv(aGrp)

      ! Test, if this Be' A' combination is to be run on this node
      if (BeAID(myRank,beGrp,aGrp) /= 0) then
        if (printkey >= 10) write(u6,*) ' Chol PID, be,a',myRank,beGrp,aGrp

        ! Extract H2(a',I) <- T1(a,I)
        call ExtT1(wrk(PosH2),wrk(PosT1o),dima,adda)

        ! read Q(be',u,i,a'),K(be',u,i,a')
        LunName = Tmp1Name(beGrp,aGrp)
        dim_1 = dimbe*dima*no*no
        call GetX(wrk(PosQ),dim_1,LunAux,LunName,1,0)
        call GetX(wrk(PosK),dim_1,LunAux,LunName,0,1)

        ! read V1(m,be',a')<- [V2(m,bea') or V2(m,a',be')] <- L2(m,be',a')
        ! in cases be <= b through interstep V2
        if (beGrp > aGrp) then
          LunName = L2Name(beGrp,aGrp)
          dim_1 = nc*dima*dimbe
          call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)
        else if (beGrp == aGrp) then
          LunName = L2Name(beGrp,aGrp)
          dim_1 = nc*dima*(dimbe+1)/2
          call GetX(wrk(PosV2),dim_1,LunAux,LunName,1,1)
          dim_1 = dima*(dimbe+1)/2
          call Exp1(wrk(PosV2),wrk(PosV1),nc,dim_1,dimbe)
        else
          LunName = L2Name(aGrp,beGrp)
          dim_1 = nc*dima*dimbe
          call GetX(wrk(PosV2),dim_1,LunAux,LunName,1,1)
          call Map3_132(wrk(PosV2),wrk(PosV1),nc,dima,dimbe)
        end if

        ! V3(a',u) <- V1(T)(m,be',a')(L2) . M3(m,be',u)
        wrk(PosV3:PosV3+dima*no-1) = Zero
        call mc0c1at3b(nc*dimbe,dima,nc*dimbe,no,dima,no,dima,nc*dimbe,no,wrk(PosV1),wrk(PosM3),wrk(PosV3))

        !T17.f Add T1n(a',u) <<- V3(a',u)
        call AdT_T17(wrk(PosT1n),wrk(PosV3),nv,dima,no,adda,One)

        ! V3(be',a',i,u) <-  V1(m,be',a')(L2) . M1(m,i,u)
        dim_1 = dima*dimbe*no*no
        wrk(PosV3:PosV3+dim_1-1) = Zero
        call mc0c1at3b(nc,dimbe*dima,nc,no*no,dimbe*dima,no*no,dima*dimbe,nc,no*no,wrk(PosV1),wrk(PosM1),wrk(PosV3))

        ! map V2(be',u,i,a') <-  V3(be',a',i,u)
        call Map4_1432(wrk(PosV3),wrk(PosV2),dimbe,dima,no,no)

        ! Q(be',u,i,a')  <<- -V2(be',u,i,a')
        wrk(PosQ:PosQ+dim_1-1) = wrk(PosQ:PosQ+dim_1-1)-wrk(PosV2:PosV2+dim_1-1)

        ! Gvv(be',a') <<- 2 sum(i) [V2(be',i,i,a')
        call AdV_G2(wrk(PosGvv),wrk(PosV2),nv,dimbe,dima,no,addbe,adda,Two)

        ! K(be',u,i,a')  <<-  V2(be',u,i,a')
        wrk(PosK:PosK+dim_1-1) = wrk(PosK:PosK+dim_1-1)+wrk(PosV2:PosV2+dim_1-1)

        ! read M5(m,i,a') <- L1(m,i,a')
        LunName = L1Name(aGrp)
        dim_1 = dima*nc*no
        call GetX(wrk(PosM5),dim_1,LunAux,LunName,1,1)

        ! V2(be',u,i,a') <-  M2(m,be',u ) . M5(m,i,a')(L1)
        dim_1 = dima*dimbe*no*no
        wrk(PosV2:PosV2+dim_1-1) = Zero
        call mc0c1at3b(nc,dimbe*no,nc,no*dima,dimbe*no,dima*no,dimbe*no,nc,dima*no,wrk(PosM2),wrk(PosM5),wrk(PosV2))

        ! Q(be',u,i,a') <<- 2 V2(be',u,i,a')
        wrk(PosQ:PosQ+dim_1-1) = wrk(PosQ:PosQ+dim_1-1)+Two*wrk(PosV2:PosV2+dim_1-1)

        ! write Q and K submatrix to corresponding files
        LunName = Tmp1Name(beGrp,aGrp)
        dim_1 = dimbe*dima*no*no
        call SaveX(wrk(PosQ),dim_1,LunAux,LunName,1,0)
        call SaveX(wrk(PosK),dim_1,LunAux,LunName,0,1)
        !@@
        !call Chck_K(wrk(PosK),dimbe,addbe,dima,adda)
        !call Chck_Q(wrk(PosQ),dimbe,addbe,dima,adda)
        !@@

        ! Gvv(be',a') <<- - sum(i) [V2(be,i,i,a')
        call AdV_G2(wrk(PosGvv),wrk(PosV2),nv,dimbe,dima,no,addbe,adda,-One)

        ! Q array (already released) will be served for X

        !T26.1f Map Q(a',U,be',V) <-  2 . V2(be',V,U,a')
        !       i.e. here a' is a first index, that means X (Q) is saved
        !       in different order as K,Q
        call Map4_3421(wrk(PosV2),wrk(PosQ),dimbe,no,no,dima)
        dim_1 = dima*no*no*dimbe
        wrk(PosQ:PosQ+dim_1-1) = Two*wrk(PosQ:PosQ+dim_1-1)

        !T27.2 Map H3(j,a') <- H2(a',j)(T1)
        call Map2_21(wrk(PosH2),wrk(PosH3),dima,no)

        !T27.3f V2(be',v,u,a') <- - V4(be',v,u,j) . H3(j,a')
        dim_1 = dimbe*no*no*dima
        wrk(PosV2:PosV2+dim_1-1) = Zero
        dim_1 = dimbe*no*no
        call mc0c2a3b(dim_1,no,no,dima,dim_1,dima,dim_1,no,dima,wrk(PosV4),wrk(PosH3),wrk(PosV2))

        !T289.1 V1(m,j,v) <- M4(m,j,v) (L1)
        dim_1 = nc*no*no
        wrk(PosV1:PosV1+dim_1-1) = wrk(PosM4:PosM4+dim_1-1)

        !T289.2 V1(m,j,v) <<- M1(m,j,v)
        dim_1 = nc*no*no
        wrk(PosV1:PosV1+dim_1-1) = wrk(PosV1:PosV1+dim_1-1)+wrk(PosM1:PosM1+dim_1-1)

        !T289.3 V3(j,v,u,a') <- V1(T) (m,j,v) . M5(m,u,a')(L1)
        dim_1 = no*no*no*dima
        wrk(PosV3:PosV3+dim_1-1) = Zero
        call mc0c1at3b(nc,no*no,nc,no*dima,no*no,no*dima,no*no,nc,no*dima,wrk(PosV1),wrk(PosM5),wrk(PosV3))

        !T289.4f V2(be',v,u,a') <<- - H1(be',j)(T1) . V3(j,v,u,a')
        dim_1 = no*no*dima
        call mc0c2a3b(dimbe,no,no,dim_1,dimbe,dim_1,dimbe,no,dim_1,wrk(PosH1),wrk(PosV3),wrk(PosV2))

        !T2G Map V1(a',U,be',V) <- V2(be',V,U,a')
        call Map4_3421(wrk(PosV2),wrk(PosV1),dimbe,no,no,dima)

        !T2G Add Q(a',u,be',v)(X) <<- V1(a',U,be',V) (Faktor 2 koli X)
        dim_1 = dimbe*dima*no*no
        wrk(PosQ:PosQ+dim_1-1) = wrk(PosQ:PosQ+dim_1-1)+Two*wrk(PosV1:PosV1+dim_1-1)

        !T2G write X(a',u,be',v) = Q(a',u,be',v)
        LunName = Tmp2Name(aGrp,beGrp)
        dim_1 = dimbe*dima*no*no
        Xyes(aGrp,beGrp) = 1
        call SaveX(wrk(PosQ),dim_1,LunAux,LunName,1,1)

      end if
      adda = adda+dima
    end do

    !T1G Add T1n(be,u) <<- H4(be',u)
    call AdT_T17(wrk(PosT1n),wrk(PosH4),nv,dimbe,no,addbe,One)

  end if
  addbe = addbe+dimbe
end do

#ifdef _MOLCAS_MPP_
!## Synchronizacny bod:
! Allreduce Goo,Gvv
dim_1 = no*no
call gadgop(wrk(PosGoo),dim_1,'+')
dim_1 = nv*nv
call gadgop(wrk(PosGvv),dim_1,'+')
#endif

!Goo1.1f Add Goo(i,u) <<- Hoo(i,u)
dim_1 = no*no
wrk(PosGoo:PosGoo+dim_1-1) = wrk(PosGoo:PosGoo+dim_1-1)+wrk(PosHoo:PosHoo+dim_1-1)

!Gvv1.1 Map V1(be,a) <- Hvv(a,be)
call Map2_21(wrk(PosHvv),wrk(PosV1),nv,nv)

!Gvv1.2f Add Gvv(be,a) <<- V1(be,a)
dim_1 = nv*nv
wrk(PosGvv:PosGvv+dim_1-1) = wrk(PosGvv:PosGvv+dim_1-1)+wrk(PosV1:PosV1+dim_1-1)

!Gvv2.1 Map V2(i,a) <- Fvo(a,i)
call Map2_21(wrk(PosFvo),wrk(PosV2),nv,no)

!Gvv2.2 Calc V1(be,a) <- t1(be,i) . V2(i,a)
dim_1 = nv*nv
wrk(PosV1:PosV1+dim_1-1) = Zero
call mc0c1a3b(nv,no,no,nv,nv,nv,nv,no,nv,wrk(PosT1o),wrk(PosV2),wrk(PosV1))

!Gvv2.3f Add Gvv(be,a) <<- - V1(be,a)
dim_1 = nv*nv
wrk(PosGvv:PosGvv+dim_1-1) = wrk(PosGvv:PosGvv+dim_1-1)-wrk(PosV1:PosV1+dim_1-1)

!@@
!call Chck_Goo(wrk(PosGoo))
!call Chck_Gvv(wrk(PosGvv))
!@@

return

end subroutine o3v3chol
