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

subroutine Reord_chcc(wrk,wrksize,NaGrpR,NaSGrpR,NchBlk,LunAux)
! This routine does:
! 1) Read local CD1 file of Cholesky vectors from MC
! 2) Reorder ChV from pq,mloc to mloc,pq
!     Make L0-L2 files (local, dimensioned as ncLoc)
!     (if needed: (JoinLkey=1)
! 3) Make I0-I3 integrals
!    Calc E MP2
!     Make T20 (fist estimation (ai|jb)/Dijab
!     Make L0-L2 files (Global, dimensioned as nc)
!     (if needed: (JoinLkey=2)
!
!  for integral based approach (intkey=1)
!  5.1) W3 (vv|vo) integrals
!  5.2) W4 (vv|vv) integrals
!  5.3) Make L0-L2 files (Global, dimensioned as nc)
!     (if needed: JoinLkey=3)

!1 Structure of files, where selected group of (pq|rs) are
!  stored (V'O|OO) - I1 ; (V'O|V'O) - I2 ; (V'V'|OO) - I3
!
!  (IJ |KL)  I0intg
!
!  (A'I|JK)  I1inxx xx - Group of A'
!
!  (A'I|B'J) I2xxyy xx - Group of A'
!                   yy - Group of B'
!
!  (A'B'|IJ) I3xxyy xx - Group of A'
!                   yy - Group of B'
!
!2 Structure of Cholesky vector files
!
!  L0(m,IJ)    L0vctr  I>=J
!
!  L1(m,I ,A') L1vcxx xx - Group of A'
!@@   kokot som, ze som to takto urobil, prerobit na L(m,a,i) to treba
!
!  L2(m,A'B')  L2xxyy xx - Group of A', A'>=B'
!                     yy - Group of B'
!
!3 Structure of Amplitude file
!  t2(A'B',I,J)  T2xxyy xx - Group of A'
!                       yy - Group of B'
!
! Memory requirements:
!  cholesky:
! V1   - max {ov'ov'; nbas.nbas.m'; ov'm; v'v'm;  oom}
! V2   - max {ov'ov'; v'v'm, ov'm; oom}
! V3   - max {ov'm; oom}
! V4   -      oom
!  integral based:
! V1   - max {ov'ov'; nbas.nbas.m'; ov'm; v'v'm; oom; V"V"V"V"}
! V2   - max {ov'ov'; v'v'm, ov'm; oom}
! V3   - max {ov'm; oom; V'V'M}
! V4   - oom
! M1   - V"V"m
! M2   - max {V"V"M; OV"M)

use Index_Functions, only: nTri_Elem
use chcc_global, only: DimGrpa, DimGrpaR, DimSGrpa, generkey, GrpaLow, GrpaUp, I0Name, I1Name, I2Name, I3Name, intkey, JoinLkey, &
                       L0Name, L1Name, L2Name, nc, nfr, no, nv, PosFoo, PosFree, PosFvo, PosFvv, PosOE, printkey, T2Name
#ifdef _MOLCAS_MPP_
use chcc_global, only: InqW3, InqW4, NChLoc, W34DistKey
use Para_Info, only: MyRank
#endif
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: wrksize, NaGrpR, NaSGrpR, NChBlk, LunAux
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp) :: abGrp, abSGrp, adda, addapp, addb, addbpp, addc, addcpp, addd, adddpp, aGrp, aSGrp, bGrp, bSGrp, bSGrpUp, &
                     cdGrp, cdSGrp, cGrp, ChLow, ChUp, cSGrp, dGrp, dim_1, dima, dimab, dimabpp, dimapp, dimb, dimbpp, dimc, &
                     dimcd, dimcdpp, DimCh(100), dimci, dimcpp, dimd, dimdpp, dimij, dSGrp, dSGrpUp, i, idisk, idum(1), LunChVF, &
                     maxdim, mdGrpa, mdGrpbe, mdSGrpa, mdSGrpbe, NaGrp, NaSGrp, NbeGrp, NbeSgrp, nbs, NCh, ncLoc, nOrbE, Nv4Ints, &
                     PosM1, PosM2, PosT, PosV1, PosV2, PosV3, PosV4, PosX, w3abcy, w3aby, w4abcdy, w4aby
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: NSGrp
#endif
real(kind=wp) :: e2, e2os
logical(kind=iwp) :: Found
character(len=24) :: Label
character(len=10) :: LunName10
character(len=8) :: LunName8
character(len=6) :: LunName
integer(kind=iwp), external :: isfreeunit

! Def parameters
call DefParReord(NaGrpR,maxdim)
NaGrp = NaGrpR
NaSGrp = NaSGrpR
NbeGrp = NaGrpR
NbeSGrp = NaSGrpR
call DefParo2v4(NaGrp,NbeGrp,NaSGrp,NbeSgrp,mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe)
if (printkey >= 10) write(u6,*) ' Maxdim',maxdim,mdSGrpa

! Distribute memory
PosT = PosFree
call DistMemReord(maxdim,mdSGrpa,NchBlk,PosV1,PosV2,PosV3,PosV4,PosM1,PosM2,PosT)
if (printkey >= 10) write(u6,*) ' Last Value :',PosT,wrksize
if (PosT > wrksize) then
  !mp write(u6,*) ' Nieje dobre - Reord_chcc, Dr. Ch. Kokotopuloss',
  write(u6,*) ' Not Enough memory in Reord_chcc step! Increase large and/or small segmentation ', &
              real(PosT,kind=wp)/real(wrksize,kind=wp)
  call abend()
end if

! Get Orbital energies

!nOrbE = nfr+no+nv ! wrong size if ndel/=0
call Get_iArray('nBas',idum,1) ! must read always nBas fr runf
nOrbE = idum(1)
Label = 'OrbE'
call qpg_dArray(Label,Found,nOrbE)
if ((.not. Found) .or. (nOrbE == 0)) call SysAbendMsg('get_orbe','Did not find:',Label)
call Get_dArray(Label,wrk(PosOE),nOrbE)
if (printkey >= 10) then ! toto som si nie isty
  do i=1,nfr+no+nv
    write(u6,*) i,wrk(PosOE+i-1)
  end do
end if

! skip frozen OE
PosOE = PosOE+nfr

! Make Foo,Fvv,Fov
! N.B. Ked bude OE(q) ine ako Fqq,  alebo F nediagonalny,
!      bude treba urobit inak

dim_1 = no*no
wrk(PosFoo:PosFoo+dim_1-1) = Zero
dim_1 = nv*nv
wrk(PosFvv:PosFvv+dim_1-1) = Zero
dim_1 = nv*no
wrk(PosFvo:PosFvo+dim_1-1) = Zero

! Escape, if Reord is not needed

if (generkey == 0) return

! ------- part 1, read  data form _CD file
!         and distributed in the form L(p',q',m') ------

! open _CD file and initialize parameters

!mp <new 21/04/09
LunChVF = 80
LunChVF = isfreeunit(LunChVF)
!mp >
call DaName_mf_wa(LunChVF,'CD1tmp')

#ifdef _MOLCAS_MPP_
ncLoc = NChLoc(myRank)
#else
ncLoc = nc
#endif

Nch = 1
ChLow = 1
ChUp = NChBlk
DimCh(1) = ChUp-ChLow+1
idisk = 1
!mp
if (printkey >= 10) then
  write(u6,*)
  write(u6,*) 'ncLoc ',ncLoc
  write(u6,*)
end if
!mp

!* read the block into V1(p,q,m') from _CD1
do
  dim_1 = (no+nv)*(no+nv)*DimCh(NCh)
  if (printkey >= 10) write(u6,*) 'Read _CD1',DimCh(NCh),dim_1,idisk
  call ddafile(LunChVF,2,wrk(PosV1),dim_1,idisk)

  !1.1 Extract L0(ij,m')

  dimc = DimCh(NCh)
  nbs = no+nv
  dimij = nTri_Elem(no)
  call Ext_L0(wrk(PosV1),wrk(PosV2),no,dimij,dimc,nbs)

  LunName = L0Name
  dim_1 = dimij*dimc
  if (ChLow == 1) then
    call SaveX(wrk(PosV2),dim_1,LunAux,LunName,1,1)
  else
    call SaveX(wrk(PosV2),dim_1,LunAux,LunName,3,1)
  end if

  !1.2 Extract L1(i,a',m')

  adda = 0
  do aGrp=1,NaGrpR
    dima = DimGrpaR(aGrp)

    dimc = DimCh(NCh)
    nbs = no+nv
    call Ext_L1(wrk(PosV1),wrk(PosV2),no,dima,dimc,adda+no,nbs)

    LunName = L1Name(aGrp)
    dim_1 = dima*no*dimc
    if (ChLow == 1) then
      call SaveX(wrk(PosV2),dim_1,LunAux,LunName,1,1)
    else
      call SaveX(wrk(PosV2),dim_1,LunAux,LunName,3,1)
    end if

    adda = adda+dima
  end do

  !1.3 Extract L2(a'b',m')

  adda = 0
  do aGrp=1,NaGrpR
    dima = DimGrpaR(aGrp)

    addb = 0
    do bGrp=1,aGrp
      dimb = DimGrpaR(bGrp)

      dimc = DimCh(NCh)
      nbs = no+nv
      if (aGrp == bGrp) then
        dimab = nTri_Elem(dima)
        call Ext_L2s(wrk(PosV1),wrk(PosV2),dima,dimab,dimc,adda+no,addb+no,nbs)
      else
        dimab = dima*dimb
        call Ext_L2u(wrk(PosV1),wrk(PosV2),dima,dimb,dimc,adda+no,addb+no,nbs)
      end if

      LunName = L2Name(aGrp,bGrp)
      dim_1 = dimab*dimc
      if (ChLow == 1) then
        call SaveX(wrk(PosV2),dim_1,LunAux,LunName,1,1)
      else
        call SaveX(wrk(PosV2),dim_1,LunAux,LunName,3,1)
      end if

      addb = addb+dimb
    end do
    adda = adda+dima
  end do

  !* Upgrade parameters, if needed
  if (ChUp >= ncLoc) exit
  Nch = Nch+1
  ChLow = ChLow+DimCh(Nch-1)
  if ((ChUp+NChBlk) > ncLoc) then
    ChUp = ncLoc
  else
    ChUp = ChUp+NChBlk
  end if
  DimCh(Nch) = ChUp-ChLow+1
end do

! close _CD file

call DaClos(LunChVF)

! ------- part 2, read data in the form L(p',q',m')
!         and redistributed in the form L(m,p'q')   --------

!2.1 Reorder L0(ij,m') to L0(ml,ij)

LunName = L0Name
dimij = nTri_Elem(no)

!2.1.1 get L0(ij,ml) into V1(ij,ml)

if (NCh == 1) then
  ! all in one
  dim_1 = dimij*ncLoc
  call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)

else if (NCh == 2) then
  ! all in two
  PosX = PosV1
  dim_1 = dimij*DimCh(1)
  call GetX(wrk(PosX),dim_1,LunAux,LunName,1,0)
  PosX = PosV1+dim_1
  dim_1 = dimij*DimCh(2)
  call GetX(wrk(PosX),dim_1,LunAux,LunName,0,1)

else
  ! more than 2 records
  PosX = PosV1
  dim_1 = dimij*DimCh(1)
  call GetX(wrk(PosX),dim_1,LunAux,LunName,1,0)
  do i=2,NCh-1
    PosX = PosX+dim_1
    dim_1 = dimij*DimCh(i)
    call GetX(wrk(PosX),dim_1,LunAux,LunName,0,0)
  end do
  PosX = PosX+dim_1
  dim_1 = dimij*DimCh(NCh)
  call GetX(wrk(PosX),dim_1,LunAux,LunName,0,1)

end if

!2.1.2 map V1(ij,ml) -> V2(ml,ij)

call Map2_21(wrk(PosV1),wrk(PosV2),dimij,ncLoc)

!2.1.3 write L0 = V2(ml,ij) back into file

dim_1 = dimij*ncLoc
call SaveX(wrk(PosV2),dim_1,LunAux,LunName,1,1)

!2.2 Reorder L1(i,a',m') to L1(ml,i,a')

do aGrp=1,NaGrpR
  dima = DimGrpaR(aGrp)
  LunName = L1Name(aGrp)

  !2.2.1 get L1(i,a',ml) into V1(i,a',ml)

  if (NCh == 1) then
    ! all in one
    dim_1 = no*dima*ncLoc
    call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)

  else if (NCh == 2) then
    ! all in two
    PosX = PosV1
    dim_1 = no*dima*DimCh(1)
    call GetX(wrk(PosX),dim_1,LunAux,LunName,1,0)
    PosX = PosV1+dim_1
    dim_1 = no*dima*DimCh(2)
    call GetX(wrk(PosX),dim_1,LunAux,LunName,0,1)

  else
    ! more than 2 records
    PosX = PosV1
    dim_1 = no*dima*DimCh(1)
    call GetX(wrk(PosX),dim_1,LunAux,LunName,1,0)
    do i=2,NCh-1
      PosX = PosX+dim_1
      dim_1 = no*dima*DimCh(i)
      call GetX(wrk(PosX),dim_1,LunAux,LunName,0,0)
    end do
    PosX = PosX+dim_1
    dim_1 = no*dima*DimCh(NCh)
    call GetX(wrk(PosX),dim_1,LunAux,LunName,0,1)

  end if

  !2.2.2 map V1(i,a',ml) -> V2(ml,i,a')

  call Map2_21(wrk(PosV1),wrk(PosV2),no*dima,ncLoc)

  !2.2.3 write L1 = V2(ml,i,a') back into file

  dim_1 = no*dima*ncLoc
  call SaveX(wrk(PosV2),dim_1,LunAux,LunName,1,1)

end do

!2.3 Reorder L2(a'b',m') to L2(ml,a'b')

do aGrp=1,NaGrpR
  dima = DimGrpaR(aGrp)
  do bGrp=1,aGrp
    dimb = DimGrpaR(bGrp)
    if (aGrp == bGrp) then
      dimab = nTri_Elem(dima)
    else
      dimab = dima*dimb
    end if
    LunName = L2Name(aGrp,bGrp)

    !2.3.1  get L2(a'b',ml) into V1(a'b',ml)

    if (NCh == 1) then
      ! all in one
      dim_1 = dimab*ncLoc
      call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)

    else if (NCh == 2) then
      ! all in two
      PosX = PosV1
      dim_1 = dimab*DimCh(1)
      call GetX(wrk(PosX),dim_1,LunAux,LunName,1,0)
      PosX = PosV1+dim_1
      dim_1 = dimab*DimCh(2)
      call GetX(wrk(PosX),dim_1,LunAux,LunName,0,1)

    else
      ! more than 2 records
      PosX = PosV1
      dim_1 = dimab*DimCh(1)
      call GetX(wrk(PosX),dim_1,LunAux,LunName,1,0)
      do i=2,NCh-1
        PosX = PosX+dim_1
        dim_1 = dimab*DimCh(i)
        call GetX(wrk(PosX),dim_1,LunAux,LunName,0,0)
      end do
      PosX = PosX+dim_1
      dim_1 = dimab*DimCh(NCh)
      call GetX(wrk(PosX),dim_1,LunAux,LunName,0,1)

    end if

    !2.3.2 map V1(a'b',ml) -> V2(ml,a'b')

    call Map2_21(wrk(PosV1),wrk(PosV2),dimab,ncLoc)

    !2.3.3 write L2 = V2(ml,a'b') back into file

    dim_1 = dimab*ncLoc
    call SaveX(wrk(PosV2),dim_1,LunAux,LunName,1,1)

  end do
end do

!2.4 reconstructing L0-L2 (Global, m=nc) ---
!    from Local L0-L2 (ml=NChLoc(myRank))
!    if needed (JoinLKey=1)

if (JoinLkey == 1) then
  if (ncLoc < nc) then
    call JoinLvec(wrk,wrksize,PosV1,PosV2,NaGrpR,LunAux)
    ncLoc = nc
  end if
end if

! ------- part 3, produce I integrals from Cholesky vectors
!
!3.1 Generate I0(ij,kl)

!3.124.1 read V4(ml,ij) <- L0(ml,ij)
LunName = L0name
dimij = nTri_Elem(no)
call GetX(wrk(PosV4),ncLoc*dimij,LunAux,LunName,1,1)

!3.1.2 V1(ij,kl) = V4(T)(ml,ij) . V4(ml,kl)
dim_1 = dimij*dimij
wrk(PosV1:PosV1+dim_1-1) = Zero
call mc0c1at3b(ncLoc,dimij,ncLoc,dimij,dimij,dimij,dimij,ncLoc,dimij,wrk(PosV4),wrk(PosV4),wrk(PosV1))

#ifdef _MOLCAS_MPP_
!## Synchronizacny bod
!3.1.3 Allreduce V1
if (ncLoc < nc) then
  dim_1 = dimij*dimij
  call gadgop(wrk(PosV1),dim_1,'+')
end if
#endif

!3.1.4 write V1(ij,kl)
LunName = I0name
dim_1 = dimij*dimij
call SaveX(wrk(PosV1),dim_1,LunAux,LunName,1,1)

!3.23 generate I1(a'i,jk), I2(a'i,b'j) :-)

e2 = Zero
e2os = Zero

addb = 0
do bGrp=1,NaGrpR
  dimb = DimGrpaR(bGrp)

  !3.23.2 read V2(ml,i,b') <- L1(ml,i,b')
  LunName = L1name(bGrp)
  dim_1 = ncLoc*no*dimb
  call GetX(wrk(PosV2),dim_1,LunAux,LunName,1,1)

  !3.23.3 Map V3(ml,b',i) <- V2(ml,i,b')
  call Map3_132(wrk(PosV2),wrk(PosV3),ncLoc,no,dimb)

  !3.2.4 calc V1(b',i,kl) <<- V3(T)(ml,b',i) . V4(ml,kl)
  !Bug dim_1 = no*dima*dimij
  dim_1 = no*dimb*dimij
  wrk(PosV1:PosV1+dim_1-1) = Zero
  call mc0c1at3b(ncLoc,dimb*no,ncLoc,dimij,dimb*no,dimij,dimb*no,ncLoc,dimij,wrk(PosV3),wrk(PosV4),wrk(PosV1))

# ifdef _MOLCAS_MPP_
  !## Synchronizacny bod
  !3.2.5 Allreduce V1
  if (ncLoc < nc) then
    dim_1 = no*dimb*dimij
    call gadgop(wrk(PosV1),dim_1,'+')
  end if
# endif

  !3.2.6 write V1(b',i,kl)
  LunName = I1name(bGrp)
  dim_1 = no*dimb*dimij
  call SaveX(wrk(PosV1),dim_1,LunAux,LunName,1,1)

  adda = 0
  do aGrp=1,NaGrpR
    dima = DimGrpaR(aGrp)
    if (printkey >= 10) write(u6,*) aGrp,bGrp,dima,dimb

    !3.3.4 read V1(ml,i,a)
    LunName = L1name(aGrp)
    dim_1 = ncLoc*no*dima
    call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)

    !3.3.5 Map V2(ml,a',i) <- V1(ml,i,a')
    call Map3_132(wrk(PosV1),wrk(PosV2),ncLoc,no,dima)

    !3.3.6 V1(a',i,b',j) <<- V2(T)(ml,a',i) . V3(ml,b',j)
    dim_1 = no*no*dima*dimb
    wrk(PosV1:PosV1+dim_1-1) = Zero
    call mc0c1at3b(ncLoc,dima*no,ncLoc,dimb*no,dima*no,dimb*no,dima*no,ncLoc,dimb*no,wrk(PosV2),wrk(PosV3),wrk(PosV1))

#   ifdef _MOLCAS_MPP_
    !## Synchronizacny bod
    !3.3.7 Allreduce V1
    if (ncLoc < nc) then
      dim_1 = no*no*dima*dimb
      call gadgop(wrk(PosV1),dim_1,'+')
    end if
#   endif

    !3.3.8 write V1(a',i,b',j)
    LunName = I2name(aGrp,bGrp)
    dim_1 = no*dima*no*dimb
    call SaveX(wrk(PosV1),dim_1,LunAux,LunName,1,1)

    !3.3.x cierny vypocet E2
    call CVE2(wrk(PosV1),wrk(PosOE),dima,dimb,adda+no,addb+no,no,e2,e2os)

    !3.3.x Make V2(ab',i,j) <- V1(a',i,b',j)/Dijab
    if (aGrp == bGrp) then
      call MkT20p(wrk(PosV2),wrk(PosV1),wrk(PosOE),dima,adda+no,no)
    else
      call MkT20u(wrk(PosV2),wrk(PosV1),wrk(PosOE),dima,dimb,adda+no,addb+no,no)
    end if

    !3.3.x Save V2 - T2(0) into T2Name
    LunName = T2Name(aGrp,bGrp)
    if (aGrp == bGrp) then
      dim_1 = nTri_Elem(dima)*no*no
    else
      dim_1 = dima*dimb*no*no
    end if
    call SaveX(wrk(PosV2),dim_1,LunAux,LunName,1,1)

    adda = adda+dima
  end do

  addb = addb+dimb
end do

write(6,*)
write(u6,92) ' E2 MP2    :',e2
write(u6,92) ' E2 ss     :',e2-e2os
write(u6,92) ' E2 os     :',e2os
write(u6,*)

!3.4 generate I3(a'b'|ij)

adda = 0
do aGrp=1,NaGrpR
  dima = DimGrpaR(aGrp)

  addb = 0
  do bGrp=1,aGrp
    dimb = DimGrpaR(bGrp)

    !3.4.2 read V2(ml,a'b')
    LunName = L2name(aGrp,bGrp)
    if (aGrp == bGrp) then
      dimab = nTri_Elem(dima)
    else
      dimab = dima*dimb
    end if
    call GetX(wrk(PosV2),ncLoc*dimab,LunAux,LunName,1,1)

    !3.4.3 V1(a'b',ij) <<- V2(T)(ml,a'b') . V4(ml,ij)
    wrk(PosV1:PosV1+dimab*dimij-1) = Zero
    call mc0c1at3b(ncLoc,dimab,ncLoc,dimij,dimab,dimij,dimab,ncLoc,dimij,wrk(PosV2),wrk(PosV4),wrk(PosV1))

#   ifdef _MOLCAS_MPP_
    !## Synchronizacny bod
    !3.4.4 Allreduce V1
    if (ncLoc < nc) then
      dim_1 = dimij*dimab
      call gadgop(wrk(PosV1),dim_1,'+')
    end if
#   endif

    !3.4.5 write V1(a',i,b',j)
    LunName = I3name(aGrp,bGrp)
    dim_1 = dimij*dimab
    call SaveX(wrk(PosV1),dim_1,LunAux,LunName,1,1)

    adda = adda+dima
  end do

  addb = addb+dimb
end do

!3.5 reconstructing L0-L2 (Global, m=nc) ---
!    from Local L0-L2 (ml=NChLoc(myRank))
!    if needed (JoinLKey=2)

if (JoinLkey == 2) then
  if (ncLoc < nc) then
    call JoinLvec(wrk,wrksize,PosV1,PosV2,NaGrpR,LunAux)
    ncLoc = nc
  end if
end if

! in W3 and W4 files are not generated, finish
if (intkey == 0) return

! ------- part 5, produce (vv|vo) and (vv|vv) integrals
!         for integral based approach
!         all terminology in o2v4 language
!
! Extra cost: Mult. read - Na*(Na+1)/2 . L2(m,cd)
!                        - Na*(Na+1)/2 . L1(m,ci)
!
! N.B.
! - Ak je L2 uz spojene (ncLoc=nc) a W34DistKey=1
!   potom sa obskakuju ani nepisu tie kombinacie, ktore na
!   tomto node netreba (vhodne pre velke nProcs aj pri rychlej
!   sieti, pre pomalu siet aj pre mensie nProcs)
! - Ak L2 este nieje spojene (ncLoc<nc) a W34DistKey=1
!   potom sa nic neobskakuje, vsetky W3 a W4 sa pocitaju
!   a allreducuju, ale sa nepisu na disk ak ich na tomto
!   node netreba, cim sa setri miesto
!   (vhodne pre mensie nProcs a rychlu siet)
! - Ak L2 este nieje spojene (ncLoc<nc) a W34DistKey=0
!   potom sa nic neobskakuje, vsetky W3 a W4 sa pocitaju
!   a allreducuju, a pisu sa na disk vsade,
!   cim sa nesetri miesto, ale vsade su vsetky W3,4
!   (zbytocny pripad, pokial sa nechceme hrat s load-ballance
!   a nepotrebujeme mat vsade vsetko)
! - Paralelny pripad (ncLoc=nc) a W34DistKey=0
!   bude pocitat iba kokot
! - Osetrene na paralelny Nprocs=1, aj na skalarny case

#ifdef _MOLCAS_MPP_

if (W34DistKey == 1) then
  ! case: Ditributed W34 files
  !*.1 make a map, which W3 and W4 files are needed on this node
  !    i.e. def InqW3, InqW4
  call Xo2v4ctl(NaGrp,NaSGrp)

else
  ! case: All W34 files on each node
  !*.1 set InqW3,InqW4 - True
  NSGrp = NaGrp*NaSGrp
  InqW3(1:nTri_Elem(NSGrp),1:NSGrp) = .true.
  InqW4(1:nTri_Elem(NSGrp),1:nTri_Elem(NSGrp)) = .true.

end if

#endif

Nv4Ints = 0

! cycle over a'>=b' groups
adda = 0
do aGrp=1,NaGrp
  dima = DimGrpa(aGrp)
  addb = 0
  do bGrp=1,aGrp
    dimb = DimGrpa(bGrp)
    abGrp = nTri_Elem(aGrp-1)+bGrp

    if (printkey >= 10) write(u6,*) ' W3 + W4 ',aGrp,bGrp

    ! test, if at least one file of W3/W4 integrals
    ! needs to be calculated on this node for given a',b'
    ! skip, if there is nothing to calculate for this a'b'
    ! on this node and L is joined (ncLoc=nc)
    call DefW34y(aGrp,bGrp,w3aby,w4aby,NaGrp)
    if ((w3aby+w4aby /= 0) .or. (ncLoc /= nc)) then

      ! read V2(ml,a'b') = L2(ml,a'b')
      LunName = L2name(aGrp,bGrp)
      if (aGrp == bGrp) then
        dimab = nTri_Elem(dima)
      else
        dimab = dima*dimb
      end if
      call GetX(wrk(PosV2),ncLoc*dimab,LunAux,LunName,1,1)

      ! ---- Part 5.1 - Generation of (ab|ci) integrals ----
      !      (Svincaja morda)

      ! skip, if there is no W3 file to calculate for this a'b'
      ! on this node and L is joined (ncLoc=nc)
      if ((w3aby /= 0) .or. (ncLoc /= nc)) then

        ! cycle over c' groups
        addc = 0
        do cGrp=1,NaGrp
          dimc = DimGrpa(cGrp)

          ! test, if at least one file of W3 integrals
          ! needs to be calculated on this node for given a',b',c'
          ! and skip, if there is nothing to calculate for this a',b',c'
          ! on this node and L is joined (ncLoc=nc)
          call DefW3y(aGrp,bGrp,cGrp,w3abcy)
          if ((w3abcy /= 0) .or. (ncLoc /= nc)) then

            ! Get V3(ml,c',i)
            ! Read V1(ml,i,c') <- L1(ml,i,c)
            dimci = dimc*no
            LunName = L1Name(cGrp)
            call GetX(wrk(PosV1),ncLoc*dimci,LunAux,LunName,1,1)
            ! Map V3(ml,c',i) <- V1(ml,i,c')
            call Map3_132(wrk(PosV1),wrk(PosV3),ncLoc,no,dimc)

            ! cycle over a">=b" subgroups
            addapp = 0
            do aSGrp=GrpaLow(aGrp),GrpaUp(aGrp)
              dimapp = DimSGrpa(aSGrp)
              if (aGrp == bGrp) then
                bSGrpUp = aSGrp
              else
                bSGrpUp = GrpaUp(bGrp)
              end if
              addbpp = 0
              do bSGrp=GrpaLow(bGrp),bSGrpUp
                dimbpp = DimSGrpa(bSGrp)
                abSGrp = nTri_Elem(aSGrp-1)+bSGrp

                ! Extract M1(ml,a"b") <- V2(ml,a'b')
                if (aSGrp == bSGrp) then
                  dimabpp = nTri_Elem(dimapp)
                else
                  dimabpp = dimapp*dimbpp
                end if
                call Ext_W4(wrk(PosV2),wrk(PosM1),ncLoc,dima,dimb,dimab,dimapp,dimbpp,dimabpp,addapp,addbpp,aGrp,bGrp,aSGrp,bSGrp)

                ! cycle over c" subgroups
                addcpp = 0
                do cSGrp=GrpaLow(cGrp),GrpaUp(cGrp)
                  dimcpp = DimSGrpa(cSGrp)

#                 ifdef _MOLCAS_MPP_
                  ! skip, if W3(a"b",c"i) is not needed on this node
                  ! and L is joined (ncLoc=nc)
                  if ((InqW3(abSGrp,cSGrp)) .or. (ncLoc /= nc)) then
#                 endif

                    ! Extract M2(ml,c",i) <- V3(ml,c',i)
                    call Ext_W3(wrk(PosV3),wrk(PosM2),ncLoc,no,dimc,dimcpp,addcpp)

                    ! Calc V1(a"b",c"i) <- M1(T)(ml,a"b") . M2(ml,c",i)
                    dim_1 = dimcpp*no
                    wrk(PosV1:PosV1+dimabpp*dim_1-1) = Zero
                    call mc0c1at3b(ncLoc,dimabpp,ncLoc,dim_1,dimabpp,dim_1,dimabpp,ncLoc,dim_1,wrk(PosM1),wrk(PosM2),wrk(PosV1))

#                   ifdef _MOLCAS_MPP_
                    !## Synchronizacny bod
                    ! Allreduce V1
                    if (ncLoc < nc) then
                      dim_1 = dimabpp*dimcpp*no
                      call gadgop(wrk(PosV1),dim_1,'+')
                    end if
                    ! skip, if W3(a"b",c"i) is not needed on this node
                    if (InqW3(abSGrp,cSGrp)) then
#                   endif

                      ! Create Proper LunName8
                      call MkNameV3(aSGrp,bSGrp,cSGrp,'W3',LunName8)

                      ! Write integral block to proper file
                      !open(unit=LunAux,file=LunName8,form='unformatted')
                      call MOLCAS_BinaryOpen_Vanilla(LunAux,LunName8)
                      call wri_chcc(LunAux,dimabpp*dimcpp*no,wrk(PosV1))
                      close(LunAux)
#                 ifdef _MOLCAS_MPP_
                    end if
                  end if
#                 endif

                  ! end cycle over c" subgroups
                  addcpp = addcpp+dimcpp
                end do

                ! end cycle over a">=b" subgroups
                addbpp = addbpp+dimbpp
              end do
              addapp = addapp+dimapp
            end do
          end if

          ! end cycle over c' groups
          addc = addc+dimc
        end do
      end if

      ! ---- Part 5.2 - Generation of (ab|cd) integrals ----
      !     (Kobyljacaja sraka)

      ! skip, if there is no W4 file to calculate for this a'b'
      ! on this node and L is joined (ncLoc=nc)
      if ((w4aby /= 0) .or. (ncLoc /= nc)) then

        ! cycle over c'>=d' groups
        addc = 0
        do cGrp=1,NaGrp
          dimc = DimGrpa(cGrp)
          addd = 0
          do dGrp=1,cGrp
            dimd = DimGrpa(dGrp)
            cdGrp = nTri_Elem(cGrp-1)+dGrp
            if (cGrp <= aGrp) then

              ! test, if at least one file of W4 integrals
              ! needs to be calculated on this node for given a',b',c',d'
              ! and skip, if there is nothing to calculate for this
              ! a',b',c',d' on this node and L is joined (ncLoc=nc)
              call DefW4y(aGrp,bGrp,cGrp,dGrp,w4abcdy)
              if ((w4abcdy /= 0) .or. (ncLoc /= nc)) then

                ! read V3(ml,c'd') = L2(ml,c'd')
                if (abGrp == cdGrp) then
                  dimcd = dimab
                  wrk(PosV3:PosV3+ncLoc*dimcd-1) = wrk(PosV2:PosV2+ncLoc*dimcd-1)
                else
                  LunName = L2name(cGrp,dGrp)
                  if (cGrp == dGrp) then
                    dimcd = nTri_Elem(dimc)
                  else
                    dimcd = dimc*dimd
                  end if
                  call GetX(wrk(PosV3),ncLoc*dimcd,LunAux,LunName,1,1)
                end if

                ! cycle over a">=b" subgroups
                addapp = 0
                do aSGrp=GrpaLow(aGrp),GrpaUp(aGrp)
                  dimapp = DimSGrpa(aSGrp)
                  if (aGrp == bGrp) then
                    bSGrpUp = aSGrp
                  else
                    bSGrpUp = GrpaUp(bGrp)
                  end if
                  addbpp = 0
                  do bSGrp=GrpaLow(bGrp),bSGrpUp
                    dimbpp = DimSGrpa(bSGrp)
                    abSGrp = nTri_Elem(aSGrp-1)+bSGrp

                    ! test, if at least one file of W4 integrals
                    ! needs to be calculated on this node for given a",b",c',d'
                    ! and skip, if there is nothing to calculate for this
                    ! a",b",c',d' on this node and L is joined (ncLoc=nc)
                    call DefW4y2(aSGrp,bSGrp,cGrp,dGrp,w4abcdy)
                    if ((w4abcdy /= 0) .or. (ncLoc /= nc)) then

                      ! Extract M1(ml,a"b") <- V2(ml,a'b')
                      if (aSGrp == bSGrp) then
                        dimabpp = nTri_Elem(dimapp)
                      else
                        dimabpp = dimapp*dimbpp
                      end if
                      call Ext_W4(wrk(PosV2),wrk(PosM1),ncLoc,dima,dimb,dimab,dimapp,dimbpp,dimabpp,addapp,addbpp,aGrp,bGrp,aSGrp, &
                                  bSGrp)

                      ! cycle over c">=d" subgroups
                      addcpp = 0
                      do cSGrp=GrpaLow(cGrp),GrpaUp(cGrp)
                        dimcpp = DimSGrpa(cSGrp)
                        if (cGrp == dGrp) then
                          dSGrpUp = cSGrp
                        else
                          dSGrpUp = GrpaUp(dGrp)
                        end if
                        adddpp = 0
                        do dSGrp=GrpaLow(dGrp),dSGrpUp
                          dimdpp = DimSGrpa(dSGrp)
                          cdSGrp = nTri_Elem(cSGrp-1)+dSGrp
                          if (cdSGrp <= abSGrp) then

#                           ifdef _MOLCAS_MPP_
                            ! skip, if W4(a"b",c"d") is not needed on this node
                            ! and L is joined (ncLoc=nc)
                            if ((InqW4(abSGrp,cdSGrp)) .or. (ncLoc /= nc)) then
#                           endif

                              ! Extract M2(ml,a"b") <- V3(ml,a'b')
                              if (cSGrp == dSGrp) then
                                dimcdpp = nTri_Elem(dimcpp)
                              else
                                dimcdpp = dimcpp*dimdpp
                              end if
                              call Ext_W4(wrk(PosV3),wrk(PosM2),ncLoc,dimc,dimd,dimcd,dimcpp,dimdpp,dimcdpp,addcpp,adddpp,cGrp, &
                                          dGrp,cSGrp,dSGrp)

                              ! Calc V1(a"b",c"d") <- M1(T)(ml,a"b") . M2(ml,c"d")
                              dim_1 = dimabpp*dimcdpp
                              wrk(PosV1:PosV1+dim_1-1) = Zero
                              call mc0c1at3b(ncLoc,dimabpp,ncLoc,dimcdpp,dimabpp,dimcdpp,dimabpp,ncLoc,dimcdpp,wrk(PosM1), &
                                             wrk(PosM2),wrk(PosV1))

#                             ifdef _MOLCAS_MPP_
                              !## Synchronizacny bod
                              ! Allreduce V1
                              if (ncLoc < nc) then
                                dim_1 = dimabpp*dimcdpp
                                call gadgop(wrk(PosV1),dim_1,'+')
                              end if
                              ! skip, if W4(a"b",c"d") not needed on this node
                              if (InqW4(abSGrp,cdSGrp)) then
#                             endif

                                ! Def proper LunName10
                                call MkNameV4(aSGrp,bSGrp,cSGrp,dSGrp,'W4',LunName10)

                                ! Write integral block to proper file
                                !open(unit=LunAux,file=LunName10,form='unformatted')
                                call MOLCAS_BinaryOpen_Vanilla(LunAux,LunName10)
                                call wri_chcc(LunAux,dimabpp*dimcdpp,wrk(PosV1))
                                close(LunAux)
                                Nv4Ints = Nv4Ints+dimabpp*dimcdpp
#                           ifdef _MOLCAS_MPP_
                              end if
                            end if
#                           endif
                          end if

                          ! end cycle over c">=d" subgroups
                          adddpp = adddpp+dimdpp
                        end do
                        addcpp = addcpp+dimcpp
                      end do
                    end if

                    ! end cycle over a">=b" subgroups
                    addbpp = addbpp+dimbpp
                  end do
                  addapp = addapp+dimapp
                end do
              end if
            end if

            ! end cycle over c'>=d' groups
            addd = addd+dimd
          end do
          addc = addc+dimc
        end do
      end if
    end if

    ! end cycle over a'>=b' groups
    addb = addb+dimb
  end do
  adda = adda+dima
end do

write(u6,*) ' V4 ints compressing factor: ',real(nv*nv*nv*nv,kind=wp)/Nv4Ints

!5.3 reconstructing L0-L2 (Global, m=nc) ---
!    if needed (JoinLKey=3)

if (JoinLkey == 3) then
  if (ncLoc < nc) then
    call JoinLvec(wrk,wrksize,PosV1,PosV2,NaGrpR,LunAux)
    ncLoc = nc
  end if
end if

return

92 format(a12,1x,f15.12)

end subroutine Reord_chcc
