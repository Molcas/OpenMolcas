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
!
! Subroutine for calculation of spectroscopic constants for
! diatomic molecules from the raw vibrational term values
! obtained in Vibrot. A least square fit to these data is used
!
! ********** MOLCAS-Release 91 05 01 **********

subroutine Spectc(Req,E0,Atom1,Atom2,nE)

use Vibrot_globals, only: J1A, J2A, Vibwvs, iadrsp
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Three, Half, OneHalf, auTocm, auToeV, Angstrom
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: Req, E0
character(len=2), intent(in) :: Atom1, Atom2
integer(kind=iwp), intent(in) :: nE
real(kind=wp) :: ChkSum, B1, B2, A11, A12, A22, Xj, Xj1, Xj2, Xj3, Xj4, Xd, Det, xnv, xnv2, xnv3, xnv4, dc, diff, Tvjc, Tvjnv, &
                 tmax, bc, X(3), A(3,3), G(3)
real(kind=wp) :: alphae, betae, gammae, dele, be, we, wexe, weye, Re, dE, d0
real(kind=wp), allocatable :: T(:,:), Tvj(:,:), F(:), Fc(:), dF(:), B(:), D(:), Gv(:), GV2(:)
integer(kind=iwp) :: nv, i, j, k, ind, jst, jend, jrng, jind(5)

ChkSum = Zero

call mma_allocate(Tvj,[1,nE],[J1A,J2A],label='Tvj')

do j=J1A,J2A
  call DDafile(Vibwvs,2,Tvj(:,j),nE,iadrsp(j-J1A+1))
end do

! Compute rotaional constants B and D for each vibrational band

if (J2A-J1A < 3) then
  write(u6,*) 'SPECTC Error: J2A-J1A < 3'
  write(u6,*) ' SPECTC requires a range of rotational quanta'
  write(u6,*) ' in order to fit spectroscopic constants.'
  write(u6,*) ' The upper and lower limits of J must differ by'
  write(u6,*) '  at least 3.'
  write(u6,'(1x,a,2i6)') 'J2A,J1A:',J2A,J1A
  call Abend()
end if
call mma_allocate(F,[J1A+1,J2A],label='F')
call mma_allocate(Fc,[J1A+1,J2A],label='Fc')
call mma_allocate(dF,[J1A+1,J2A],label='dF')
call mma_allocate(B,nE,label='B')
call mma_allocate(D,nE,label='D')
do nv=1,nE
  ! Term values for this band
  do j=J1A+1,J2A
    F(j) = Tvj(nv,j)-Tvj(nv,J1A)
    F(j) = F(j)*auTocm
  end do
  B1 = Zero
  B2 = Zero
  A11 = Zero
  A12 = Zero
  A22 = Zero
  do j=J1A+1,J2A
    Xj1 = j*(j+1)-J1A*(J1A+1)
    Xj2 = Xj1**2
    Xj3 = Xj1**3
    Xj4 = Xj2**2
    B1 = B1+F(j)*Xj1
    B2 = B2-F(j)*Xj2
    A11 = A11+Xj2
    A12 = A12+Xj3
    A22 = A22+Xj4
  end do

  Xd = A12/A22
  Det = A11-Xd*A12
  B(nv) = (B1+B2*Xd)/Det
  D(nv) = (B2*A11/A22+B1*Xd)/Det
  do j=J1A+1,J2A
    Xj = j*(j+1)-J1A*(J1A+1)
    Fc(j) = B(nv)*Xj-D(nv)*Xj**2
    dF(j) = F(j)-Fc(j)
  end do
  write(u6,990) nv-1,B(nv),D(nv)
  do j=J1A+1,J2A
    write(u6,991) j,F(j),Fc(j),dF(j)
  end do
end do
write(u6,992)
do nv=1,nE
  write(u6,991) nv,B(nv),D(nv)
end do
call mma_deallocate(F)
call mma_deallocate(Fc)
call mma_deallocate(dF)

! Compute spectroscopic constants de (dele) and betae
! by a least square fit

dele = Zero
betae = Zero
if (nE == 1) then
  de = d(1)
else
  A11 = nE
  B1 = Zero
  B2 = Zero
  A12 = Zero
  A22 = Zero
  do nv=1,nE
    xnv = nv-Half
    B1 = B1+D(nv)
    B2 = B2+D(nv)*xnv
    A12 = A12+xnv
    A22 = A22+xnv**2
  end do
  Det = A11*A22-A12*A12
  dele = (B1*A22-B2*A12)/Det
  betae = (B2*A11-B1*A12)/Det
end if
write(u6,994) dele,betae
do nv=1,nE
  dc = dele+betae*(nv-Half)
  diff = D(nv)-dc
  write(u6,991) nv-1,D(nv),dc,diff
end do

call mma_deallocate(D)

! Spectroscopic constants Be,Alphe, and Gammae from rotational
! constants B(nv)

X(:) = Zero
if (nE == 1) then
  X(1) = B(1)
else if (nE == 2) then
  X(1) = OneHalf*B(1)-Half*B(2)
  X(2) = B(2)-B(1)
else
  G(:) = Zero
  A(:,:) = Zero
  A(1,1) = nE
  do nv=1,nE
    xnv = nv-Half
    xnv2 = xnv**2
    xnv3 = xnv**3
    xnv4 = xnv2**2
    G(1) = G(1)+B(nv)
    G(2) = G(2)+B(nv)*xnv
    G(3) = G(3)+B(nv)*xnv2
    A(1,2) = A(1,2)+xnv
    A(1,3) = A(1,3)+xnv2
    A(2,3) = A(2,3)+xnv3
    A(3,3) = A(3,3)+xnv4
  end do
  A(2,1) = A(1,2)
  A(2,2) = A(1,3)
  A(3,1) = A(1,3)
  A(3,2) = A(2,3)
  call dminv(3,3,A)
  do i=1,3
    X(i) = Zero
    do k=1,3
      X(i) = X(i)+A(i,k)*G(k)
    end do
  end do
end if
be = X(1)
alphae = -X(2)
gammae = X(3)
write(u6,995) be,alphae,gammae
do nv=1,nE
  bc = be-alphae*(nv-Half)+gammae*(nv-Half)**2
  diff = B(nv)-bc
  write(u6,991) nv-1,B(nv),bc,diff
end do

call mma_deallocate(B)

! vibrational constants we,wexe and weye from band origins

X(:) = Zero
if (nE == 1) then
  X(1) = Two*Tvj(1,J1A)
else if (nE == 2) then
  X(1) = Three*Tvj(1,J1A)-Tvj(2,J1A)/Three
  x(2) = Two*Tvj(2,J1A)/Three-Two*Tvj(1,J1A)
else
  G(:) = Zero
  A(:,:) = Zero
  do nv=1,nE
    xnv = nv-Half
    xnv2 = xnv**2
    xnv3 = xnv**3
    G(1) = G(1)+Tvj(nv,J1A)*xnv
    G(2) = G(2)+Tvj(nv,J1A)*xnv2
    G(3) = G(3)+Tvj(nv,J1A)*xnv3
    A(1,1) = A(1,1)+xnv2
    A(1,2) = A(1,2)+xnv3
    A(1,3) = A(1,3)+xnv2**2
    A(2,3) = A(2,3)+xnv2*xnv3
    A(3,3) = A(3,3)+xnv3**2
  end do
  A(2,1) = A(1,2)
  A(2,2) = A(1,3)
  A(3,1) = A(1,3)
  A(3,2) = A(2,3)
  call dminv(3,3,a)
  do i=1,3
    X(i) = Zero
    do k=1,3
      X(i) = X(i)+A(i,k)*G(k)
    end do
  end do
end if
we = X(1)*auTocm
wexe = X(2)*auTocm
weye = X(3)*auTocm
write(u6,996) we,wexe,weye
do nv=1,nE
  xnv = nv-Half
  Tvjc = we*xnv+wexe*xnv**2+weye*xnv**3
  Tvjnv = auTocm*Tvj(nv,J1A)
  diff = Tvjnv-Tvjc
  write(u6,991) nv-1,Tvjnv,Tvjc,diff
end do

! print output of spectroscopic constants

write(u6,1000) Atom1,Atom2
write(u6,1100) J1A,J2A,0,nE-1

Re = Req*Angstrom
dE = -E0*auToeV
d0 = dE-Tvj(1,J1A)*auToeV

write(u6,1200) Re,dE,d0,we,wexe,weye,be,alphae,gammae,dele,betae
ChkSum = ChkSum+we+wexe

! compute term values for check

call mma_allocate(T,[1,nE],[J1A,J2A],label='T')

ind = 0
do j=J1A,J2A
  Xj = j*(j+1)
  do nv=1,nE
    ind = ind+1
    xnv = nv-Half
    T(nv,j) = we*xnv+wexe*xnv**2+weye*xnv**3+xj*(be-alphae*xnv+gammae*xnv**2)-xj**2*(dele+betae*xnv)
  end do
end do

! compute max deviation

tmax = Zero
do j=J1A,J2A
  do nv=1,nE
    Tvj(nv,j) = Tvj(nv,j)*auTocm
    diff = abs(T(nv,j)-Tvj(nv,j))
    if (diff >= tmax) tmax = diff
  end do
end do
ChkSum = ChkSum+tmax
write(u6,1300) tmax

! print output of term values

write(u6,1400)
jst = J1A
do
  jend = jst+size(jind)-1
  if (jend > J2A) jend = J2A
  jrng = jend-jst+1
  do i=1,jrng
    jind(i) = jst+i-1
  end do
  write(u6,1500) jind(1:jrng)
  write(u6,1510)
  do nv=1,nE
    write(u6,1600) nv-1,(Tvj(nv,j),T(nv,j),j=jst,jend)
  end do
  jst = jend+1
  if (jend >= J2A) exit
end do

! observed g-values (Tvj for j=0)

call mma_allocate(Gv,nE,label='Gv')
call mma_allocate(Gv2,nE-1,label='Gv')

do i=1,nE
  Gv(i) = Tvj(i,J1A)
  ChkSum = ChkSum+Gv(i)
end do

! compute delta g(v+1/2) values

do i=1,nE-1
  Gv2(i) = Gv(i+1)-Gv(i)
end do

! print Gv and Gv2

write(u6,1700)
do i=1,nE-1
  write(u6,1800) i-1,Gv(i)
  write(u6,1810) Gv2(i)
end do
write(u6,1800) nE-1,Gv(nE)

call mma_deallocate(Gv)
call mma_deallocate(Gv2)
call mma_deallocate(T)

call mma_deallocate(Tvj)

call Add_Info('VIBROT_SPECTC',[ChkSum],1,2)
!write(u6,*) 'Spectc: ChkSum',ChkSum

return

990 format(/5x,'Rotational constants for vibrational quantum number',i3 &
           /5x,'B=',e13.6,' cm-1     D=',e13.6,' cm-1' &
           /5x,'Observed and computed term values (cm-1)')
991 format(5x,i3,3e20.6)
992 format(/5x,'Rotational constants B(nv) and D(nv) in cm-1')
994 format(/5x,'Spectroscopic constants De=',e13.6,' cm-1  Betae=',e13.6,' cm-1' &
           /5x,'Observed and computed D values')
995 format(/5x,'Spectroscopic constants Be,Alphae and Gammae' &
           /5x,'Be=',e13.6,' cm-1    Alphae=',e13.6,' cm-1    Gammae=',e13.6 &
           /5x,'Observed and computed B values')
996 format(/5x,'Vibrational constants we  =',e13.6,' cm-1' &
           /5x,'                      wexe=',e13.6,' cm-1' &
           /5x,'                      weye=',e13.6,' cm-1' &
           /5x,'Observed and computed band origins')
1000 format(//5x,'Spectroscopic constants for ',2a2)
1100 format(//5x,'Range of J-values used in fit',2i3, &
            /5x,'Range of v-values used in fit',2i3)
1200 format(///5x,'Re(a)',15x,f8.4, &
            /5x,'De(ev)',14x,f8.4, &
            /5x,'D0(ev)',14x,f8.4, &
            /5x,'we(cm-1)',7x,e13.6, &
            /5x,'wexe(cm-1)',5x,e13.6, &
            /5x,'weye(cm-1)',5x,e13.6, &
            /5x,'Be(cm-1)',7x,e13.6, &
            /5x,'Alphae(cm-1)',3x,e13.6, &
            /5x,'Gammae(cm-1)',3x,e13.6, &
            /5x,'Dele(cm-1)',5x,e13.6, &
            /5x,'Betae(cm-1)',4x,e13.6)
1300 format(//5x,'Max deviation in term values is',e10.2,' cm(-1)')
1400 format(//5x,'Term values(observed and computed) in cm(-1)')
1500 format(/1x,'J-value',7x,i3,19x,i3,19x,i3,19x,i3,19x,i3)
1510 format(/1x,'v-value')
1600 format(4x,i2,2x,5(f8.2,2x,f8.2,4x))
1700 format(//5x,'observed G-values in cm(-1)' &
            //5x,'v',6x,'G(v)',5x,'deltaG(v+1/2)')
1800 format(4x,i2,f11.2)
1810 format(17x,f13.2)

end subroutine Spectc
