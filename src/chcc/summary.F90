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

subroutine summary(wrk,wrksize,NvGrp,LunAux,maxdim,E1,E2,E2os)
! this routine does:
! T1
!   - make T1o = T1n/Den
!   - calc E1 = FTau
! T2:
!   - Join T2 form Tmp2files-(X,Y) and Tmp3files
!   - make T2o in T2files = T2n/Den
!   - make Tau from new amplitudes
!   - calc E2,E2os = WTau
!
! Scheme of recreation of T2n:
!
! 1) contributions from X and Y matrices
! T2n(be',ga',u,v) <<-
! C1                + 1/2 X(be',u,ga',v)
! C2                + 1/2 X(ga',v,be',u)
! C3                - 1/2 Y(be',u,ga',v)
! C4                - 1/2 Y(ga',v,be',u)
! C5                - 1   Y(ga',u,be',v)
! C6                - 1   Y(be',v,ga',u)
!
! 2) contributions from T2n1,T2n2 (T2+,T2-)
! T2n(be',ga',u,v) <<-
! C7                    T2+(be'ga',uv)
! C8                    T2-(be'ga',uv)
!
! Energy contributions:
! E1  = sum(a,i)       [  2 Fvo(a,i) . T1(a,i)               ]
! E2  = sum(a,b,i,j)   [ {2(ai|bj) - (aj|bi)} . Tau(a,b,i,j) ]
! E2os= sum(a,b,i,j)   [   (ai|bj)            . Tau(a,b,i,j) ]
!
! Memory requirements:
! V1-3 - {v'ov'o}
! H1,2 - {v'o}

use Index_Functions, only: nTri_Elem
use chcc_global, only: DimGrpv, DimSGrpbe, GrpbeLow, GrpbeUp, I2Name, no, nv, PosFree, PosFvo, PosOE, PosT1n, PosT1o, T2Name, &
                       T2o2v4yes, Tmp2Name, Tmp3Name, Xyes, XYyes
#ifdef _MOLCAS_MPP_
use Para_Info, only: nProcs
#endif
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, NvGrp, LunAux, maxdim
real(kind=wp), intent(inout) :: wrk(wrksize)
real(kind=wp), intent(out) :: E1, E2, E2os
integer(kind=iwp) :: addbe, addbepp, addga, addgapp, beGrp, beSGrp, dim_1, dim_2, dimbe, dimbepp, dimga, dimgapp, gaGrp, gaSGrp, &
                     PosH1, PosH2, PosT, PosV1, PosV2, PosV3
real(kind=wp) :: Ehlp, Eoshlp
character(len=6) :: LunName

! Distribute Memory
PosT = PosFree
call DistMemSum(maxdim,PosV1,PosV2,PosV3,PosH1,PosH2,PosT)

! Operations on T1 amplitudes

! Divide T1n by denominators
call T1_div(wrk(PosT1n),wrk(PosOE),no,nv)

! Set T1o <- T1n
dim_1 = no*nv
wrk(PosT1o:PosT1o+dim_1-1) = wrk(PosT1n:PosT1n+dim_1-1)

! Calc E1
call Energy_E1(wrk(PosT1n),wrk(PosFvo),no,nv,EHlp)
E1 = Two*EHlp
!mp E1 = Two*E1

! Operations on T2 amplitudes

! cycle over be'=ga'

E2 = Zero
E2os = Zero

addbe = 0
do beGrp=1,NvGrp
  dimbe = DimGrpv(beGrp)

  ! Extract H1(be',I) <- T1(be,I)
  call ExtT1(wrk(PosH1),wrk(PosT1o),dimbe,addbe)

  ! vanish V3
  dim_1 = dimbe*dimbe*no*no
  wrk(PosV3:PosV3+dim_1-1) = Zero

  !1.1 Add contribution C1-6

  ! read V1(be',u,ga',v) <- X(be',u,ga',v)
  !      V2(be',u,ga',v) <- Y(be',u,ga',v)
  LunName = Tmp2Name(beGrp,beGrp)
  dim_1 = dimbe*dimbe*no*no
  if (XYyes(beGrp,beGrp) == 1) then
    call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,0)
    call GetX(wrk(PosV2),dim_1,LunAux,LunName,0,1)
  else if (Xyes(beGrp,beGrp) == 1) then
    call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)
    wrk(PosV2:PosV2+dim_1-1) = Zero
  else
    wrk(PosV1:PosV1+dim_1-1) = Zero
    wrk(PosV2:PosV2+dim_1-1) = Zero
  end if

  ! V3(be',ga',u,v) <-
  ! C1                + 1/2 V1(be',u,ga',v)
  ! C2                + 1/2 V1(ga',v,be',u)
  ! C3                - 1/2 V2(be',u,ga',v)
  ! C4                - 1/2 V2(ga',v,be',u)
  ! C5                - 1   V2(ga',u,be',v)
  ! C6                - 1   V2(be',v,ga',u)
  call MkT_CAlld(wrk(PosV3),wrk(PosV1),wrk(PosV2),dimbe,no)
  !1.2 Add contribution C7,C8

  ! cycle over all subgroups of (be">=ga")
  !@
  !if (.false.) then
  !@
  addbepp = 0
  do beSGrp=GrpbeLow(beGrp),GrpbeUp(beGrp)
    dimbepp = DimSGrpbe(beSGrp)
    addgapp = 0
    do gaSGrp=GrpbeLow(beGrp),beSGrp
      dimgapp = DimSGrpbe(gaSGrp)

      if (beSGrp == gaSGrp) then
        ! read V1(bega",uv+) <- T2+(bega",uv+)
        !      V2(bega",uv-) <- T2-(bega",uv-)
        LunName = Tmp3Name(beSGrp,beSGrp)
        dim_1 = nTri_Elem(dimbepp)*nTri_Elem(no)
        dim_2 = nTri_Elem(dimbepp-1)*nTri_Elem(no-1)
        if (T2o2v4yes(beSGrp,beSGrp) == 1) then
          call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,0)
          call GetX(wrk(PosV2),dim_2,LunAux,LunName,0,1)
        else
          wrk(PosV1:PosV1+dim_1-1) = Zero
          wrk(PosV2:PosV2+dim_2-1) = Zero
        end if

        ! V3(be',ga',u,v) <<- T2+(bega",uv)
        !                     T2-(bega",uv)
        call MkT_C78d(wrk(PosV3),wrk(PosV1),wrk(PosV2),dimbe,dimbepp,addbepp,no)
      else
        LunName = Tmp3Name(beSGrp,gaSGrp)
        ! read V1(be",ga",uv+) <- T2+(be",ga",uv+)
        !      V2(be",ga",uv-) <- T2-(be",ga",uv-)
        dim_1 = dimbepp*dimgapp*nTri_Elem(no)
        dim_2 = dimbepp*dimgapp*nTri_Elem(no-1)
        if (T2o2v4yes(beSGrp,gaSGrp) == 1) then
          call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,0)
          call GetX(wrk(PosV2),dim_2,LunAux,LunName,0,1)
        else
          wrk(PosV1:PosV1+dim_1-1) = Zero
          wrk(PosV2:PosV2+dim_2-1) = Zero
        end if
        ! V3(be',ga',u,v) <<- T2+(be",ga",uv)
        !                     T2-(be",ga",uv)
        call MkT_C78od(wrk(PosV3),wrk(PosV1),wrk(PosV2),dimbe,dimbe,dimbepp,dimgapp,addbepp,addgapp,no)
      end if

      addgapp = addgapp+dimgapp
    end do
    addbepp = addbepp+dimbepp
  end do
  !@
  !end if
  !@

  !1.3 T2 - V3(be',ga',u,v) are now reconstructed (for be'>=ga')
  ! Divide by denominators and completed to be'<ga' also
  call T2d_div(wrk(PosV3),wrk(PosOE),dimbe,dimbe,addbe,addbe,no,nv)

  ! Prepair reduced Set V2(be'>=ga',u,v) <- V3(be',ga',u,v)
  call MkT_red(wrk(PosV2),wrk(PosV3),dimbe,no)

# ifdef _MOLCAS_MPP_
  !## Synchronizacny bod:
  ! Allreduce V2
  dim_1 = nTri_Elem(dimbe)*no*no
  call gadgop(wrk(PosV2),dim_1,'+')

  !## For parallel runs only:
  ! Inverse expansion  V3(be',ga',u,v) <- V2(be'>=ga',u,v)
  if (nProcs > 1) call MkT_exp(wrk(PosV2),wrk(PosV3),dimbe,no)
# endif
  !@@
  !call Chck_T2n(wrk(PosV3),dimbe,addbe,dimbe,addbe,1)
  !@@

  ! Save into corresponding T2file
  LunName = T2Name(beGrp,beGrp)
  dim_1 = nTri_Elem(dimbe)*no*no
  call SaveX(wrk(PosV2),dim_1,LunAux,LunName,1,1)

  ! Make Tau in V3(be',ga',u,v)
  ! (note, that T1 are already new ones)
  ! @@ s tym 2x H1 to nieje celkom koser
  call MkTau_chcc(wrk(PosV3),wrk(PosH1),wrk(PosH1),dimbe,dimbe,no,One,One)

  ! read V1(be',I,ga',J) <- I2(be'I|ga'J)
  LunName = I2Name(beGrp,beGrp)
  dim_1 = dimbe*dimbe*no*no
  call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)

  ! Calc E2 = sum(a,b,i,j) (2V1(ai|bj)-V1(aj|bi)) . V3(a,b,i,j)
  call Energy_E2d(wrk(PosV1),wrk(PosV3),Ehlp,Eoshlp,dimbe,no)

  E2 = E2+Ehlp
  E2os = E2os+Eoshlp

  addbe = addbe+dimbe
end do

! cycle over be>ga
! reserved files:
!   V3(be',ga',u,v) = T2n(be',ga',u,v)
! free for use: V1,V2

addbe = DimGrpv(1)
do beGrp=2,NvGrp
  dimbe = DimGrpv(beGrp)

  ! Extract H1(be',I) <- T1(be,I)
  call ExtT1(wrk(PosH1),wrk(PosT1o),dimbe,addbe)

  addga = 0
  do gaGrp=1,beGrp-1
    dimga = DimGrpv(gaGrp)

    ! Extract H2(ga',I) <- T1(ga,I)
    call ExtT1(wrk(PosH2),wrk(PosT1o),dimga,addga)

    ! vanish V3
    dim_1 = dimbe*dimga*no*no
    wrk(PosV3:PosV3+dim_1-1) = Zero

    ! read V1(be',u,ga',v) <- X(be',u,ga',v)
    !      V2(be',u,ga',v) <- Y(be',u,ga',v)
    LunName = Tmp2Name(beGrp,gaGrp)
    dim_1 = dimbe*dimga*no*no
    if (XYyes(beGrp,gaGrp) == 1) then
      call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,0)
      call GetX(wrk(PosV2),dim_1,LunAux,LunName,0,1)
    else if (Xyes(beGrp,gaGrp) == 1) then
      call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)
      wrk(PosV2:PosV2+dim_1-1) = Zero
    else
      wrk(PosV1:PosV1+dim_1-1) = Zero
      wrk(PosV2:PosV2+dim_1-1) = Zero
    end if

    ! Add contribution C1,C3,C6
    ! V3(be',ga',u,v) < - + 1/2 V1(be',u,ga',v)
    !                     - 1/2 V2(be',u,ga',v)
    !                     - 1   V2(be',v,ga',u)
    call MkT_C136od(wrk(PosV3),wrk(PosV1),wrk(PosV2),dimbe,dimga,no)

    ! read V1(ga',u,be',v) <- X(ga',u,be',v)
    !      V2(ga',u,be',v) <- Y(ga',u,be',v)
    LunName = Tmp2Name(gaGrp,beGrp)
    dim_1 = dimbe*dimga*no*no
    if (XYyes(gaGrp,beGrp) == 1) then
      call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,0)
      call GetX(wrk(PosV2),dim_1,LunAux,LunName,0,1)
    else if (Xyes(gaGrp,beGrp) == 1) then
      call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)
      wrk(PosV2:PosV2+dim_1-1) = Zero
    else
      wrk(PosV1:PosV1+dim_1-1) = Zero
      wrk(PosV2:PosV2+dim_1-1) = Zero
    end if

    ! Add contribution C2,C4,C5
    ! V3(be',ga',u,v) <<- + 1/2 V1(ga',v,be',u)
    !                     - 1/2 V2(ga',v,be',u)
    !                     - 1   V2(ga',u,be',v)
    call MkT_C245od(wrk(PosV3),wrk(PosV1),wrk(PosV2),dimbe,dimga,no)

    !2.2 Add contribution C7,C8

    ! cycle over all subgroups of (be">=ga")
    addbepp = 0
    do beSGrp=GrpbeLow(beGrp),GrpbeUp(beGrp)
      dimbepp = DimSGrpbe(beSGrp)
      addgapp = 0
      do gaSGrp=GrpbeLow(gaGrp),GrpbeUp(gaGrp)
        dimgapp = DimSGrpbe(gaSGrp)

        ! read V1(be",ga",uv+) <- T2+(be",ga",uv+)
        !      V2(be",ga",uv-) <- T2-(be",ga",uv-)
        LunName = Tmp3Name(beSGrp,gaSGrp)
        dim_1 = dimbepp*dimgapp*nTri_Elem(no)
        dim_2 = dimbepp*dimgapp*nTri_Elem(no-1)
        if (T2o2v4yes(beSGrp,gaSGrp) == 1) then
          call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,0)
          call GetX(wrk(PosV2),dim_2,LunAux,LunName,0,1)
        else
          wrk(PosV1:PosV1+dim_1-1) = Zero
          wrk(PosV2:PosV2+dim_2-1) = Zero
        end if
        ! V3(be',ga',u,v) <<- T2+(be",ga",uv)
        !                     T2-(be",ga",uv)
        call MkT_C78od(wrk(PosV3),wrk(PosV1),wrk(PosV2),dimbe,dimga,dimbepp,dimgapp,addbepp,addgapp,no)

        addgapp = addgapp+dimgapp
      end do
      addbepp = addbepp+dimbepp
    end do

    ! T2 - V3(be',ga',u,v) are now completely reconstructed
    !  Divide by denominators
    call T2od_div(wrk(PosV3),wrk(PosOE),dimbe,dimga,addbe,addga,no,nv)

#   ifdef _MOLCAS_MPP_
    !## Synchronizacny bod:
    ! Allreduce V3
    dim_1 = dimbe*dimga*no*no
    call gadgop(wrk(PosV3),dim_1,'+')
#   endif
    !@@
    !call Chck_T2n(wrk(PosV3),dimbe,addbe,dimga,addga,0)
    !@@

    ! Save V3 into corresponding T2file
    LunName = T2Name(beGrp,gaGrp)
    dim_1 = dimbe*dimga*no*no
    call SaveX(wrk(PosV3),dim_1,LunAux,LunName,1,1)

    ! Map V1(ga',be',v,u) <- V3(be',ga',u,v)
    call Map4_2143(wrk(PosV3),wrk(PosV1),dimbe,dimga,no,no)

    ! Save V1 into corresponding T2file
    LunName = T2Name(gaGrp,beGrp)
    dim_1 = dimbe*dimga*no*no
    call SaveX(wrk(PosV1),dim_1,LunAux,LunName,1,1)

    ! Make Tau in V3(be',ga',u,v)
    ! (note, that T1 are already new ones)
    call MkTau_chcc(wrk(PosV3),wrk(PosH1),wrk(PosH2),dimbe,dimga,no,One,One)

    ! read V1(be',I,ga',J) <- I2(be'I|ga'J)
    LunName = I2Name(beGrp,gaGrp)
    dim_1 = dimbe*dimga*no*no
    call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)

    ! Calc E2 = sum(a,b,i,j) (2V1(ai|bj)-V1(aj|bi)) . V3(a,b,i,j)
    call Energy_E2od(wrk(PosV1),wrk(PosV3),Ehlp,Eoshlp,dimbe,dimga,no)

    E2 = E2+Two*Ehlp
    E2os = E2os+Two*Eoshlp

    addga = addga+dimga
  end do

  addbe = addbe+dimbe
end do

!@@
return

! Calc Energy from Chck
call Chck_energ()

! upgrade checkeroo T1c and T2c

! Upgrade T1
call UpG_T1(wrk(PosT1o))

! Upgrade T2
addbe = 0
do beGrp=1,NvGrp
  dimbe = DimGrpv(beGrp)

  addga = 0
  do gaGrp=1,beGrp
    dimga = DimGrpv(gaGrp)

    if (beGrp == gaGrp) then
      LunName = T2Name(beGrp,beGrp)
      dim_1 = nTri_Elem(dimbe)*no*no
      call GetX(wrk(PosV3),dim_1,LunAux,LunName,1,1)
      call UpG_T2d(wrk(PosV3),dimbe,addbe)
    else
      LunName = T2Name(beGrp,gaGrp)
      dim_1 = dimbe*dimga*no*no
      call GetX(wrk(PosV3),dim_1,LunAux,LunName,1,1)
      call UpG_T2od(wrk(PosV3),dimbe,addbe,dimga,addga)
    end if

    addga = addga+dimga
  end do

  addbe = addbe+dimbe
end do

return

end subroutine summary
