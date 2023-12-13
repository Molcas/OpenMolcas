!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Christian Pomelli                                      *
!***********************************************************************

subroutine PolyGen(MaxT,MxSph,Ptype,Pflag,TsAre,TsNum,XEN,YEN,ZEN,REN,TsEff,CV,JTR)
! Polygen: a program to generate spherical polyhedra with triangular
! faces by C.S. Pomelli (cris@ibm550.icqem.cnr.it)
! An equilateral division algorithm is used.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Four, Half, Pi
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: MaxT, MxSph, Pflag
real(kind=wp), intent(in) :: TsAre, XEN, YEN, ZEN, REN
integer(kind=iwp), intent(inout) :: TsNum
integer(kind=iwp), intent(out) :: Ptype, TsEff, JTR(3,MaxT)
real(kind=wp), intent(out) :: CV(3,MxSph)
integer(kind=iwp) :: ednew(90,100), i, ii, j, jj, l, m, n, NDI, NDP, NDT, NE0, NF, NFI, NFP, NFT, NT0, NTPT, NV, NVPT, &
                     oldtr(100,100)
real(kind=wp) :: alpha, beta, cos1, cos2, costheta, dnorm, sintheta, theta, v1(3), v2(3), v3(3)
integer(kind=iwp), allocatable :: TrNew(:,:,:)
integer(kind=iwp), parameter :: Icosa = 1, Pentakis = 2, Tetra = 3
! Data for original polyhedra
! 1) Tetrahedron
integer(kind=iwp), parameter :: NET = 6, NTT = 4, NVT = 4, &
                                teted(2,NET) = reshape([1,2,2,3,1,3,1,4,3,4,2,4],shape(teted)), &
                                tettre(3,NTT) = reshape([1,2,3,1,6,4,3,5,4,2,5,6],shape(tettre)), &
                                tettrv(3,NTT) = reshape([1,2,3,1,2,4,1,3,4,2,3,4],shape(tettrv))
! sqrt(1/3) = 0.577350269
real(kind=wp), parameter :: tetve(3,NVT) = reshape([-0.577350269_wp, 0.577350269_wp, 0.577350269_wp, &
                                                     0.577350269_wp,-0.577350269_wp, 0.577350269_wp, &
                                                    -0.577350269_wp,-0.577350269_wp,-0.577350269_wp, &
                                                     0.577350269_wp, 0.577350269_wp,-0.577350269_wp],shape(tetve))
! 2) Icosaedron
integer(kind=iwp), parameter :: NEI = 30, NTI = 20, NVI = 12, &
                                icoed(2,NEI) = reshape([1,2,2,3,1,3,2,6,1,6,3,4,1,4,4,5,1,5,5,6,3,7,2,7,6,11,2,11,7,11,4,8,3,8,7, &
                                                        8,5,9,4,9,8,9,6,10,5,10,9,10,10,11,8,12,7,12,11,12,9,12,10,12], &
                                                       shape(icoed)), &
                                icotre(3,NTI) = reshape([1,2,3,1,4,5,3,6,7,7,8,9,9,10,5,2,11,12,4,13,14,12,15,14,6,16,17,11,18,17, &
                                                         8,19,20,16,21,20,10,22,23,19,24,23,22,25,13,18,26,27,15,28,27,21,29,26, &
                                                         24,30,29,25,28,30],shape(icotre)), &
                                icotrv(3,NTI) = reshape([1,2,3,1,2,6,1,3,4,1,4,5,1,5,6,2,3,7,2,6,11,2,7,11,3,4,8,3,7,8,4,5,9,4,8, &
                                                         9,5,6,10,5,9,10,6,10,11,7,8,12,7,11,12,8,9,12,9,10,12,10,11,12], &
                                                        shape(icotrv))
! sqrt(4/5) = 0.894427191
! sqrt(1/5) = 0.447213595
! 1/2+sqrt(1/20) = 0.723606798
! 1/2-sqrt(1/20) = 0.276393202
! sqrt(1/2+sqrt(1/20)) = 0.850650808
! sqrt(1/2-sqrt(1/20)) = 0.525731112
real(kind=wp), parameter :: icove(3,NVI) = reshape([ 0.000000000_wp, 0.000000000_wp, 1.000000000_wp, &
                                                     0.276393202_wp, 0.850650808_wp, 0.447213595_wp, &
                                                    -0.723606798_wp, 0.525731112_wp, 0.447213595_wp, &
                                                    -0.723606798_wp,-0.525731112_wp, 0.447213595_wp, &
                                                     0.276393202_wp,-0.850650808_wp, 0.447213595_wp, &
                                                     0.894427191_wp, 0.000000000_wp, 0.447213595_wp, &
                                                    -0.276393202_wp, 0.850650808_wp,-0.447213595_wp, &
                                                    -0.894427191_wp, 0.000000000_wp,-0.447213595_wp, &
                                                    -0.276393202_wp,-0.850650808_wp,-0.447213595_wp, &
                                                     0.723606798_wp,-0.525731112_wp,-0.447213595_wp, &
                                                     0.723606798_wp, 0.525731112_wp,-0.447213595_wp, &
                                                     0.000000000_wp, 0.000000000_wp,-1.000000000_wp],shape(icove))
! 3) Pentakisdodecahedron
integer(kind=iwp), parameter :: NEP = 90, NTP = 60, NVP = 32, &
                                pened(2,NEP) = reshape([1,2,2,3,1,3,2,6,1,6,3,4,1,4,4,5,1,5,5,6,3,8,2,8,6,7,2,7,7,12,2,12,8,12,4, &
                                                        9,3,9,8,13,3,13,9,13,5,10,4,10,9,14,4,14,10,14,6,11,5,11,10,15,5,15,11,15, &
                                                        7,16,6,16,11,16,12,17,7,17,16,17,12,18,8,18,13,18,13,19,9,19,14,19,14,20, &
                                                        10,20,15,20,15,21,11,21,16,21,17,22,12,22,18,22,18,23,13,23,19,23,19,24, &
                                                        14,24,20,24,20,25,15,25,21,25,17,26,16,26,21,26,22,27,17,27,26,27,22,28, &
                                                        18,28,23,28,23,29,19,29,24,29,24,30,20,30,25,30,25,31,21,31,26,31,27,28, &
                                                        28,29,29,30,30,31,27,31,28,32,27,32,31,32,29,32,30,32],shape(pened)), &
                                pentre(3,NTP) = reshape([1,2,3,1,4,5,3,6,7,7,8,9,9,10,5,2,11,12,4,13,14,14,15,16,12,17,16,6,18,19, &
                                                         11,20,21,19,22,21,8,23,24,18,25,26,24,27,26,10,28,29,23,30,31,29,32,31, &
                                                         13,33,34,28,35,34,15,36,37,33,38,37,17,39,40,20,41,40,22,42,43,25,44,43, &
                                                         27,45,46,30,47,46,32,48,49,35,50,49,36,51,52,39,53,52,41,54,55,42,56,55, &
                                                         44,57,58,45,59,58,47,60,61,48,62,61,38,63,64,50,65,64,51,66,67,63,68,67, &
                                                         53,69,70,54,71,70,56,72,73,57,74,73,59,75,76,60,77,76,62,78,79,65,80,79, &
                                                         66,81,69,71,82,72,74,83,75,77,84,78,68,85,80,81,86,87,85,88,87,82,89,86, &
                                                         83,90,89,84,88,90],shape(pentre)), &
                                pentrv(3,NTP) = reshape([1,2,3,1,2,6,1,3,4,1,4,5,1,5,6,2,3,8,2,6,7,2,7,12,2,8,12,3,4,9,3,8,13,3,9, &
                                                         13,4,5,10,4,9,14,4,10,14,5,6,11,5,10,15,5,11,15,6,7,16,6,11,16,7,12,17,7, &
                                                         16,17,8,12,18,8,13,18,9,13,19,9,14,19,10,14,20,10,15,20,11,15,21,11,16, &
                                                         21,12,17,22,12,18,22,13,18,23,13,19,23,14,19,24,14,20,24,15,20,25,15,21, &
                                                         25,16,17,26,16,21,26,17,22,27,17,26,27,18,22,28,18,23,28,19,23,29,19,24, &
                                                         29,20,24,30,20,25,30,21,25,31,21,26,31,22,27,28,23,28,29,24,29,30,25,30, &
                                                         31,26,27,31,27,28,32,27,31,32,28,29,32,29,30,32,30,31,32],shape(pentrv))
! sqrt(4/5) = 0.894427191
! sqrt(1/3) = 0.577350269
! sqrt(1/5) = 0.447213595
! 1/2+sqrt(1/20) = 0.723606798
! 1/2-sqrt(1/20) = 0.276393202
! sqrt(1/2+sqrt(1/20)) = 0.850650808
! sqrt(1/2-sqrt(1/20)) = 0.525731112
! sqrt(1/3+sqrt(4/45)) = 0.794654472
! sqrt(1/3-sqrt(4/45)) = 0.187592474
! sqrt(2/3+sqrt(4/45)) = 0.982246946
! sqrt(2/3-sqrt(4/45)) = 0.607061998
! sqrt(1/6+sqrt(1/180)) = 0.491123473
! sqrt(1/6-sqrt(1/180)) = 0.303530999
! sqrt(1/2+sqrt(5/36)) = 0.934172359
! sqrt(1/2-sqrt(5/36)) = 0.356822090
real(kind=wp), parameter :: penve(3,NVP) = reshape([ 0.000000000_wp, 0.000000000_wp, 1.000000000_wp, &
                                                     0.491123473_wp, 0.356822090_wp, 0.794654472_wp, &
                                                    -0.187592474_wp, 0.577350269_wp, 0.794654472_wp, &
                                                    -0.607061998_wp, 0.000000000_wp, 0.794654472_wp, &
                                                    -0.187592474_wp,-0.577350269_wp, 0.794654472_wp, &
                                                     0.491123473_wp,-0.356822090_wp, 0.794654472_wp, &
                                                     0.894427191_wp, 0.000000000_wp, 0.447213595_wp, &
                                                     0.276393202_wp, 0.850650808_wp, 0.447213595_wp, &
                                                    -0.723606798_wp, 0.525731112_wp, 0.447213595_wp, &
                                                    -0.723606798_wp,-0.525731112_wp, 0.447213595_wp, &
                                                     0.276393202_wp,-0.850650808_wp, 0.447213595_wp, &
                                                     0.794654472_wp, 0.577350269_wp, 0.187592474_wp, &
                                                    -0.303530999_wp, 0.934172359_wp, 0.187592474_wp, &
                                                    -0.982246946_wp, 0.000000000_wp, 0.187592474_wp, &
                                                    -0.303530999_wp,-0.934172359_wp, 0.187592474_wp, &
                                                     0.794654472_wp,-0.577350269_wp, 0.187592474_wp, &
                                                     0.982246946_wp, 0.000000000_wp,-0.187592474_wp, &
                                                     0.303530999_wp, 0.934172359_wp,-0.187592474_wp, &
                                                    -0.794654472_wp, 0.577350269_wp,-0.187592474_wp, &
                                                    -0.794654472_wp,-0.577350269_wp,-0.187592474_wp, &
                                                     0.303530999_wp,-0.934172359_wp,-0.187592474_wp, &
                                                     0.723606798_wp, 0.525731112_wp,-0.447213595_wp, &
                                                    -0.276393202_wp, 0.850650808_wp,-0.447213595_wp, &
                                                    -0.894427191_wp, 0.000000000_wp,-0.447213595_wp, &
                                                    -0.276393202_wp,-0.850650808_wp,-0.447213595_wp, &
                                                     0.723606798_wp,-0.525731112_wp,-0.447213595_wp, &
                                                     0.607061998_wp, 0.000000000_wp,-0.794654472_wp, &
                                                     0.187592474_wp, 0.577350269_wp,-0.794654472_wp, &
                                                    -0.491123473_wp, 0.356822090_wp,-0.794654472_wp, &
                                                    -0.491123473_wp,-0.356822090_wp,-0.794654472_wp, &
                                                     0.187592474_wp,-0.577350269_wp,-0.794654472_wp, &
                                                     0.000000000_wp, 0.000000000_wp,-1.000000000_wp],shape(penve))
integer(kind=iwp) :: edo(max(NET,NEI,NEP),2), treo(max(NTT,NTI,NTP),3), trvo(max(NTT,NTI,NTP),3)

NT0 = 0 ! dummy initialize
NV = 0 ! dummy initialize
NE0 = 0 ! dummy initialize

! A seconda delle due opzioni di funzionamento (Area ottimale o
! numero di tessere ottimale ) vengono stabilite il tipo di poliedro
! (Ptype) e la  frequenza di divisione (NF)
!
! Pflag = 0 Numero di tessere: si genera il poliedro disponibile
! con il numero di tessere piu' prossimo a TsNum
!
! Pflag  =1  Area: si genera il poliedro con il numero di tessere tali
! da avere un area il piu' simile possibile
if (Pflag == 1) then
  TsNum = int(Four*pi*REN**2/TsAre+Half)
end if

! si controlla che TsNum <= MxTs

if (TsNum > MaxT) then
  write(u6,'(a)') 'Too many tesserae in PolyGen'
  call Abend()
end if

! calcola il valore di NF ottimale per ciascuna famiglia di poliedri.

NFI = int(sqrt(TsNum/(NTI*One))+Half)
NFP = int(sqrt(TsNum/(NTP*One))+Half)
NFT = int(sqrt(TsNum/(NTT*One))+Half)

! trova lo scostamento del numero di tessere in input

NDI = abs(TsNum-NTI*NFI**2)
NDP = abs(TsNum-NTP*NFP**2)
NDT = abs(TsNum-NTT*NFT**2)

! sceglie tipo e frequenza di divisione
! a parita' di differenza in numero di tessere sceglie quello con NF
! minore

if ((NDP <= NDI) .and. (NDP <= NDT)) then
  Ptype = Pentakis
  NF = NFP
  TsEff = NTP*NFP**2
else
  if (NDI <= NDT) then
    Ptype = Icosa
    NF = NFI
    TsEff = NTI*NFI**2
  else
    Ptype = Tetra
    NF = NFT
    TsEff = NTT*NFT**2
  end if
end if

! carica gli opportuni valori iniziali

if (Ptype == Icosa) then

  ! icosaedro

  NT0 = NTI
  NV = NVI
  NE0 = NEI
  do i=1,NT0
    trvo(i,:) = icotrv(:,i)
    treo(i,:) = icotre(:,i)
  end do
  do i=1,NV
    CV(:,i) = icove(:,i)
  end do
  do i=1,NE0
    edo(i,:) = icoed(:,i)
  end do
else if (Ptype == Pentakis) then

  ! pentakisdodecaedro

  NT0 = NTP
  NV = NVP
  NE0 = NEP
  do i=1,NT0
    trvo(i,:) = pentrv(:,i)
    treo(i,:) = pentre(:,i)
  end do
  do i=1,NV
    CV(:,i) = penve(:,i)
  end do
  do i=1,NE0
    edo(i,:) = pened(:,i)
  end do
else if (Ptype == Tetra) then

  ! tetraedro

  NT0 = NTT
  NV = NVT
  NE0 = NET
  do i=1,NT0
    trvo(i,:) = tettrv(:,i)
    treo(i,:) = tettre(:,i)
  end do
  do i=1,NV
    CV(:,i) = tetve(:,i)
  end do
  do i=1,NE0
    edo(i,:) = teted(:,i)
  end do
end if
NVPT = NV+1
if (NVPT > 1000) then
  write(u6,*) 'NVPT out of range in polygen',nvpt
  call Abend()
end if

! nuovi vertici posti lungo i vecchi spigoli
! le regole di calcolo derivano da  calcoli di algebra lineare e trigonometria
! i vertici vengono memorizzati in ordine progressivo ed riferiti
! tramite ednew al vertice di appartenenza

do j=1,NE0
  v1(:) = CV(:,edo(j,1))
  v2(:) = CV(:,edo(j,2))
  costheta = (v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3))/(sqrt(v1(1)**2+v1(2)**2+v1(3)**2)*sqrt(v2(1)**2+v2(2)**2+v2(3)**2))
  theta = acos(costheta)
  sintheta = sin(theta)
  do l=1,NF-1
    m = NF-l
    cos1 = cos(theta*l/NF)
    cos2 = cos(theta*(NF-l)/NF)
    alpha = (cos1-costheta*cos2)/sintheta**2
    beta = (cos2-costheta*cos1)/sintheta**2
    v3(:) = alpha*v1(:)+beta*v2(:)
    dnorm = sqrt(v3(1)**2+v3(2)**2+v3(3)**2)
    v3(:) = v3(:)/dnorm
    CV(:,NVPT) = v3(:)
    ednew(j,l+1) = NVPT
    NVPT = NVPT+1
    if (NVPT > 1000) then
      write(u6,*) 'NVPT out of range in polygen',nvpt
      call Abend()
    end if
  end do
end do

! nuovi vertici non posti lungo i vecchi spigoli
! a partire dai vertici in ednew secondo regole analoghe alle
! precedenti, vengono memorizzati a seconda del triangolo in trnew
! trnew(triangolo,fila,n ordine)=nvertice

! Allocate TrNew. Note that we change the order from NT0*NF*NF
! to NF*NF*NT0, (i,j,k) to (k,j,i)

call mma_allocate(TrNew,NF,NF,NT0,label='TrNew')

do j=1,NT0
  ii = treo(j,1)
  jj = treo(j,3)
  do l=3,NF
    v1(:) = CV(:,ednew(ii,l))
    v2(:) = CV(:,ednew(jj,l))
    costheta = (v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3))/(sqrt(v1(1)**2+v1(2)**2+v1(3)**2)*sqrt(v2(1)**2+v2(2)**2+v2(3)**2))
    theta = acos(costheta)
    sintheta = sin(theta)
    do n=1,l-2
      cos1 = cos(theta*n/(l-1))
      cos2 = cos(theta*(l-1-n)/(l-1))
      alpha = (cos1-costheta*cos2)/sintheta**2
      beta = (cos2-costheta*cos1)/sintheta**2
      v3(:) = alpha*v1(:)+beta*v2(:)
      dnorm = sqrt(v3(1)**2+v3(2)**2+v3(3)**2)
      v3(:) = v3(:)/dnorm
      CV(:,NVPT) = v3(:)

      !trnew(j,l,n+1) = NVPT ! Old code

      TrNew(n+1,l,j) = NVPT

      NVPT = NVPT+1
      if (NVPT > 1000) then
        write(u6,*) 'NVPT out of range in polygen',nvpt
        call Abend()
      end if
    end do
  end do
end do
NV = NVPT-1

! ora per ogni triangolo originario vengono posti nella matrice oldtr
! i numeri di ordine dei vertici originali,creati lungo gli spigoli e
! no dei vecchi triangoli secondo lo schema:
!
!         11                    La disposizione e' secondo
!         | \                   la posizione geometrica del
!         |  \                  vertice.
!         21--22
!         | \  |\               i vertici (1,1)-(NF+1,1)-(NF+1,NF+1) sono
!         |  \ | \              quelli originari.
!         31--32--33
!         | \  |  | \           i nuovi triangoli sono:
!         |  \ |  |  \          (i,j)-(i+1,j)-(i+1,j+1)
!         41--42--43--44        i=1,...,NF j=1,...,i
!
!                               (i,j)--(i,j+1),(i+1,j+1)
!                               i=2,...,NF j=1,...,i-1

NTPT = 1
do n=1,NT0

  ! 1 vecchi spigoli

  oldtr(1,1) = trvo(n,1)
  oldtr(NF+1,1) = trvo(n,2)
  oldtr(NF+1,NF+1) = trvo(n,3)

  ! 2 nuovi vertici lungo i vecchi spigoli

  do l=2,NF
    oldtr(l,1) = ednew(treo(n,1),l)
    oldtr(NF+1,l) = ednew(treo(n,2),l)
    oldtr(l,l) = ednew(treo(n,3),l)
  end do

  ! 3 nuovi vertici non lungo i vecchi spigoli

  do l=3,NF
    do m=2,l-1

      !oldtr(l,m) = trnew(n,l,m) ! old code

      oldtr(l,m) = TrNew(m,l,n)

    end do
  end do

  ! ora si creano i nuovi triangoli

  do i=1,NF
    do j=1,i
      JTR(1,NTPT) = oldtr(i,j)
      JTR(2,NTPT) = oldtr(i+1,j)
      JTR(3,NTPT) = oldtr(i+1,j+1)
      NTPT = NTPT+1
    end do
  end do
  do i=2,NF
    do j=1,i-1
      JTR(1,NTPT) = oldtr(i,j)
      JTR(2,NTPT) = oldtr(i,j+1)
      JTR(3,NTPT) = oldtr(i+1,j+1)
      NTPT = NTPT+1
    end do
  end do
end do
call mma_deallocate(TrNew)

! scrittura della JVT

do i=1,NV
  CV(1,i) = CV(1,i)*REN+XEN
  CV(2,i) = CV(2,i)*REN+YEN
  CV(3,i) = CV(3,i)*REN+ZEN
end do

return

end subroutine PolyGen
