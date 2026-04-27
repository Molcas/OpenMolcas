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

!ifdef _DEBUPRINT_
subroutine frenkelexc(Frenkeltri,ndim,nst1,nst2)

use Index_Functions, only: nTri_Elem
use frenkel_global_vars, only: excl, iTyp, jTyp, nestla, nestlb, valst
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, auToEV
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ndim, nst1, nst2
real(kind=wp), intent(inout) :: Frenkeltri(nTri_Elem(ndim))
integer(kind=iwp) :: a, b, d1lines, d2lines, I, IO, iPL, J, K, L, LUT1, ntrans, run, tri
real(kind=wp) :: DIPNORM, GS
character(len=13) :: filnam, filnam1, filnam2, filnam3
integer(kind=iwp), allocatable :: I1(:), I2(:), F1(:), F2(:)
real(kind=wp), allocatable :: D1X(:), D1Y(:), D1Z(:), D2X(:), D2Y(:), D2Z(:), Dip1(:,:,:), Dip2(:,:,:), DipFrintermx(:,:), &
                              DipFrintermx2(:,:), DipFrintermy(:,:), DipFrintermy2(:,:), DipFrintermz(:,:), DipFrintermz2(:,:), &
                              DipFrx(:,:), DipFry(:,:), DipFrz(:,:), E_FRENKEL(:), EDIFF(:), EigEn(:), EigEnHa(:), &
                              Freninterm(:,:), Freninterm2(:,:), Frenkeldia(:,:), Frenkelquad(:,:), OSCSTR(:), Rassi1(:), Rassi2(:)
character(len=8), allocatable :: Hamelembra(:), Hamelemket(:)
integer(kind=iwp), external :: iPrintLevel, isFreeUnit

call GETPRINTLEVEL()
iPL = iPrintLevel(-1)

call mma_allocate(Frenkeldia,ndim,ndim)
call mma_allocate(Freninterm,ndim,ndim)
call mma_allocate(Freninterm2,ndim,ndim)
call mma_allocate(Frenkelquad,ndim,ndim)
call mma_allocate(Rassi1,nst1)
call mma_allocate(Rassi2,nst2)
call mma_allocate(E_FRENKEL,ndim)

! read in rassi energies of both monomers to add them
! to the diagonal elements
#ifdef _DEBUGPRINT_
write(u6,*) 'iTyp=',iTyp,'jTyp=',jTyp
#endif

if (iTyp < jTyp) then
  write(u6,*) 'Monomer B was calculated first.'
  write(u6,*) 'Index 1 refers to A.'
else
  write(u6,*) 'Monomer A was calculated first.'
  write(u6,*) 'Index 1 refers to B.'
  iTyp = 2
  jTyp = 1
end if

write(filnam2,'(A,I1)') 'stE',iTyp
LuT1 = isFreeUnit(11)

call molcas_open(LuT1,filnam2)
write(u6,'(A,1X,I1,A,I3.3)') 'Number of states of system',iTyp,'=',nst1
write(u6,'(A,1X,I1,A,I3.3)') 'Number of states of system',jTyp,'=',nst2
do i=1,nst1
  read(LuT1,*) Rassi1(i)
end do
close(LuT1)

write(filnam3,'(A,I1)') 'stE',jTyp
LuT1 = isFreeUnit(11)

call molcas_open(LuT1,filnam3)
do i=1,nst2
  read(LuT1,*) Rassi2(i)
end do
close(LuT1)

! add rassi energies to diag of exciton coupl matrix
do i=1,nst1
  write(u6,*) 'State energy of system 1',Rassi1(i)
end do
do i=1,nst2
  write(u6,*) 'State energy of system 2',Rassi2(i)
end do

tri = 0
do a=1,nst1
  do b=1,nst2
    if ((a /= 1) .or. (b /= 1)) then
      if (EXCL) then
        if ((all(nestla /= a)) .or. (all(nestlb /= b))) cycle
      else
        if ((a <= valst) .and. (b <= valst)) cycle
        if ((a > valst) .and. (b > valst)) cycle
      end if
    end if
    tri = tri+1
    Frenkeltri(nTri_Elem(tri)) = Frenkeltri(nTri_Elem(tri))+Rassi1(a)+Rassi2(b)
#   ifdef _DEBUGPRINT_
    write(u6,*) 'tri counter',tri
    write(u6,*) 'add energies 1 and 2',Rassi1(a),Rassi2(b)
#   endif
  end do
end do

! subtract gs energy
gs = Frenkeltri(1)
do i=1,ndim
  Frenkeltri(nTri_Elem(i)) = Frenkeltri(nTri_Elem(i))-gs
end do
! convert to eV
Frenkeltri(:) = Frenkeltri(:)*auToEV

! make a full Hamilton matrix from the lower triangular
call square(Frenkeltri,Frenkelquad,1,ndim,ndim)
#ifdef _DEBUGPRINT_
write(u6,*) 'Frenkelquad:'
do i=1,ndim
  write(u6,'(1000ES18.8)') (Frenkelquad(j,i),j=1,ndim)
end do
#endif

call mma_allocate(Hamelembra,ndim)
call mma_allocate(Hamelemket,ndim)
tri = 0
do a=1,nst1
  do b=1,nst2
    if ((a /= 1) .or. (b /= 1)) then
      if (EXCL) then
        if ((all(nestla /= a)) .or. (all(nestlb /= b))) cycle
      else
        if ((a <= valst) .and. (b <= valst)) cycle
        if ((a > valst) .and. (b > valst)) cycle
      end if
    end if
    tri = tri+1
    write(Hamelemket(tri),'(A,I3.3,I3.3,A)') '|',a,b,'>'
    write(Hamelembra(tri),'(A,I3.3,I3.3,A)') '<',a,b,'|'
  end do
end do

write(u6,*) 'Frenkel Hamiltonian [eV]:'
write(u6,'(11X)',advance='NO')
write(u6,'(100(A,8X))',advance='YES') (Hamelemket(i),i=1,ndim)

do i=1,ndim
  write(u6,'(A,1X,100(ES16.8))') Hamelembra(i),(Frenkelquad(i,j),j=1,ndim)
end do

write(u6,*) 'Frenkel eigenvectors (column wise):'
write(u6,'(11X)',advance='NO')
write(u6,'(100(I3.3,13X))',advance='YES') (i,i=1,ndim)

call NIdiag_New(Frenkeltri,Frenkeldia,ndim,ndim)

do i=1,ndim
  write(u6,'(A,1X,100(ES16.8))') Hamelembra(i),(Frenkeldia(i,j),j=1,ndim)
end do

write(u6,*) 'Hamiltonian eigenvalues [eV]:'
do i=1,ndim
  E_FRENKEL(i) = Frenkeltri(nTri_Elem(i))
  write(u6,*) Frenkeltri(nTri_ELem(i))
end do

call Add_Info('E_FRENKEL',E_FRENKEL,ndim,6)

call mma_deallocate(E_FRENKEL)
call mma_deallocate(Hamelembra)
call mma_deallocate(Hamelemket)

#ifdef _DEBUGPRINT_
call dgemm_('N','N',ndim,ndim,ndim,One,Frenkelquad,ndim,Frenkeldia,ndim,Zero,Freninterm,ndim)

write(u6,*) 'first trafo H*U'
do i=1,ndim
  write(u6,'(1000ES18.8)') (Freninterm(i,j),j=1,ndim)
end do
call dgemm_('T','N',ndim,ndim,ndim,One,Frenkeldia,ndim,Freninterm,ndim,Zero,Freninterm2,ndim)
write(u6,*) 'second trafo U^(t)HU (should be diagonal)'
do i=1,ndim
  write(u6,'(1000ES18.8)') (Freninterm2(i,j),j=1,ndim)
end do
#endif

ntrans = nTri_Elem(ndim-1)

d1lines = 0
d2lines = 0
! search for number of lines in the dipvec files
write(filnam,'(A,I1)') 'dip_vec',iTyp
LuT1 = isFreeUnit(11)
call molcas_open(LuT1,filnam)
do
  read(LuT1,*,iostat=io)
  if (io /= 0) exit
  d1lines = d1lines+1
end do
close(LuT1)

write(filnam1,'(A,I1)') 'dip_vec',jTyp
LuT1 = isFreeUnit(11)
call molcas_open(LuT1,filnam1)
do
  read(LuT1,*,iostat=io)
  if (io /= 0) exit
  d2lines = d2lines+1
end do
close(LuT1)

#ifdef _DEBUGPRINT_
write(u6,*) 'd1lines ',d1lines
write(u6,*) 'd2lines ',d2lines
write(u6,*) 'ntrans',ntrans
#endif

call mma_allocate(I1,d1lines)
call mma_allocate(F1,d1lines)
call mma_allocate(D1X,d1lines)
call mma_allocate(D1Y,d1lines)
call mma_allocate(D1Z,d1lines)
call mma_allocate(I2,d2lines)
call mma_allocate(F2,d2lines)
call mma_allocate(D2X,d2lines)
call mma_allocate(D2Y,d2lines)
call mma_allocate(D2Z,d2lines)
call mma_allocate(Dip1,nst1,nst1,3)
call mma_allocate(Dip2,nst2,nst2,3)
call mma_allocate(EigEn,ndim)
call mma_allocate(EigEnHa,ndim)
call mma_allocate(EDIFF,ntrans)
call mma_allocate(OSCSTR,ntrans)
call mma_allocate(DipFrx,ndim,ndim)
call mma_allocate(DipFry,ndim,ndim)
call mma_allocate(DipFrz,ndim,ndim)
call mma_allocate(DipFrintermx,ndim,ndim)
call mma_allocate(DipFrintermx2,ndim,ndim)
call mma_allocate(DipFrintermy,ndim,ndim)
call mma_allocate(DipFrintermy2,ndim,ndim)
call mma_allocate(DipFrintermz,ndim,ndim)
call mma_allocate(DipFrintermz2,ndim,ndim)

write(filnam,'(A,I1)') 'dip_vec',iTyp
LuT1 = isFreeUnit(11)
call molcas_open(LuT1,filnam)
do i=1,d1lines
  read(LuT1,222) I1(i),F1(i),D1X(i),D1Y(i),D1Z(i)
end do
close(LuT1)

write(filnam1,'(A,I1)') 'dip_vec',jTyp
LuT1 = isFreeUnit(11)
call molcas_open(LuT1,filnam1)
do i=1,d2lines
  read(LuT1,222) I2(i),F2(i),D2X(i),D2Y(i),D2Z(i)
end do
close(LuT1)

! fill dipole matrices and use symmetry
Dip1(:,:,:) = Zero
do i=1,d1lines
  Dip1(I1(i),F1(i),1) = D1X(i)
  Dip1(I1(i),F1(i),2) = D1Y(i)
  Dip1(I1(i),F1(i),3) = D1Z(i)
  Dip1(F1(i),I1(i),1) = D1X(i)
  Dip1(F1(i),I1(i),2) = D1Y(i)
  Dip1(F1(i),I1(i),3) = D1Z(i)
end do

write(u6,*) 'Transition dipole vectors of system 1:'
write(u6,551) 'from','to','x','y','z'
do i=1,d1lines
  write(u6,443) I1(i),F1(i),(Dip1(I1(i),F1(i),j),j=1,3)
end do

Dip2(:,:,:) = Zero
do i=1,d2lines
  Dip2(I2(i),F2(i),1) = D2X(i)
  Dip2(I2(i),F2(i),2) = D2Y(i)
  Dip2(I2(i),F2(i),3) = D2Z(i)
  Dip2(F2(i),I2(i),1) = D2X(i)
  Dip2(F2(i),I2(i),2) = D2Y(i)
  Dip2(F2(i),I2(i),3) = D2Z(i)
end do

write(u6,*) 'Transition dipole vectors of system 2:'
write(u6,551) 'from','to','x','y','z'
do i=1,d2lines
  write(u6,443) I2(i),F2(i),(Dip2(I2(i),F2(i),j),j=1,3)
end do

write(u6,*) 'relative excitonic state eigenenergies [eV]:'
do i=1,ndim
  EigEn(i) = Frenkeltri(nTri_Elem(i))-Frenkeltri(1)
  write(u6,*) 'EigEn:',i,EigEn(i)
end do

if (iPL >= 3) then
  write(u6,*) 'relative excitonic state eigenenergies in hartree:'
  do i=1,ndim
    EigEnHa(i) = (Frenkeltri(nTri_Elem(i))-Frenkeltri(1))/auToEV
    write(u6,*) 'EigEn:',i,EigEnHa(i)
  end do
end if

! create x
a = 0
b = 0
DipFrx(:,:) = Zero
do i=1,nst1
  do k=1,nst2
    if ((i /= 1) .or. (k /= 1)) then
      if (excl) then
        if ((all(nestla /= i)) .or. (all(nestlb /= k))) cycle
      else
        if ((i <= valst) .and. (k <= valst)) cycle
        if ((i > valst) .and. (k > valst)) cycle
      end if
    end if
    a = a+1
    b = 0
    do j=1,nst1
      do l=1,nst2
        if ((j /= 1) .or. (l /= 1)) then
          if (EXCL) then
            if ((all(nestla /= j)) .or. (all(nestlb /= l))) cycle
          else
            if ((j <= valst) .and. (l <= valst)) cycle
            if ((j > valst) .and. (l > valst)) cycle
          end if
        end if
        b = b+1
        if (k == l) DipFrx(a,b) = Dip1(i,j,1)
        if (i == j) DipFrx(a,b) = Dip2(k,l,1)
      end do
    end do
  end do
end do

if (iPL >= 3) then
  write(u6,*) 'DipFrx:'
  do i=1,ndim
    write(u6,'(100ES18.8)') (DipFrx(i,j),j=1,ndim)
  end do
end if
! transform dipole matrix in dipole basis
call dgemm_('N','N',ndim,ndim,ndim,One,DipFrx,ndim,Frenkeldia,ndim,Zero,DipFrintermx,ndim)

call dgemm_('T','N',ndim,ndim,ndim,One,Frenkeldia,ndim,DipFrintermx,ndim,Zero,DipFrintermx2,ndim)

if (iPL >= 3) then
  write(u6,*) 'trafo U^(t)DU, dipole mtx in exc. basis, X'
  do i=1,ndim
    write(u6,'(100ES18.8)') (DipFrintermx2(i,j),j=1,ndim)
  end do
end if
! create y
a = 0
b = 0
DipFry(:,:) = Zero
do i=1,nst1
  do k=1,nst2
    if ((i /= 1) .or. (k /= 1)) then
      if (EXCL) then
        if ((all(nestla /= i)) .or. (all(nestlb /= k))) cycle
      else
        if ((i <= valst) .and. (k <= valst)) cycle
        if ((i > valst) .and. (k > valst)) cycle
      end if
    end if
    a = a+1
    b = 0
    do j=1,nst1
      do l=1,nst2
        if ((j /= 1) .or. (l /= 1)) then
          if (EXCL) then
            if ((all(nestla /= j)) .or. (all(nestlb /= l))) cycle
          else
            if ((j <= valst) .and. (l <= valst)) cycle
            if ((j > valst) .and. (l > valst)) cycle
          end if
        end if
        b = b+1
        if (k == l) DipFry(a,b) = Dip1(I,J,2)
        if (i == j) DipFry(a,b) = Dip2(K,L,2)
      end do
    end do
  end do
end do

if (iPL >= 3) then
  write(u6,*) 'DipFry:'
  do i=1,ndim
    write(u6,'(100ES18.8)') (DipFry(i,j),j=1,ndim)
  end do
end if

! transform dipole matrix in dipole basis
call dgemm_('N','N',ndim,ndim,ndim,One,DipFry,ndim,Frenkeldia,ndim,Zero,DipFrintermy,ndim)

call dgemm_('T','N',ndim,ndim,ndim,One,Frenkeldia,ndim,DipFrintermy,ndim,Zero,DipFrintermy2,ndim)
if (iPL >= 3) then
  write(u6,*) 'trafo U^(t)DU, dipole mtx in exc. basis, Y'
  do i=1,ndim
    write(u6,'(100ES18.8)') (DipFrintermy2(i,j),j=1,ndim)
  end do
end if

! create z
a = 0
b = 0
DipFrz(:,:) = Zero
do i=1,nst1
  do k=1,nst2
    if ((i /= 1) .or. (k /= 1)) then
      if (excl) then
        if ((all(nestla /= i)) .or. (all(nestlb /= k))) cycle
      else
        if ((i <= valst) .and. (k <= valst)) cycle
        if ((i > valst) .and. (k > valst)) cycle
      end if
    end if
    a = a+1
    b = 0
    do j=1,nst1
      do l=1,nst2
        if ((j /= 1) .or. (l /= 1)) then
          if (excl) then
            if ((all(nestla /= j)) .or. (all(nestlb /= l))) cycle
          else
            if ((j <= valst) .and. (l <= valst)) cycle
            if ((j > valst) .and. (l > valst)) cycle
          end if
        end if
        b = b+1
        if (k == l) then
          if (k > 1) cycle
          DipFrz(a,b) = Dip1(I,J,3)
        end if
        if (i == j) DipFrz(a,b) = Dip2(K,L,3)
      end do
    end do
  end do
end do
if (iPL >= 3) then
  write(u6,*) 'DipFrz:'
  do i=1,ndim
    write(u6,'(100ES18.8)') (DipFrz(i,j),j=1,ndim)
  end do
end if

! transform dipole matrix in dipole basis
call dgemm_('N','N',ndim,ndim,ndim,One,DipFrz,ndim,Frenkeldia,ndim,Zero,DipFrintermz,ndim)

call dgemm_('T','N',ndim,ndim,ndim,One,Frenkeldia,ndim,DipFrintermz,ndim,Zero,DipFrintermz2,ndim)
if (iPL >= 3) then
  write(u6,*) 'trafo U^(t)DU, dipole mtx in exc. basis, Z'
  do i=1,ndim
    write(u6,'(100ES18.8)') (DipFrintermz2(i,j),j=1,ndim)
  end do
end if

! calc and print spectrum
run = 0
write(u6,*) 'Excitonic absorption spectrum'
write(u6,3) 'from','to','excitation energy [eV]','Dx','Dy','Dz','osc.str.'
do i=1,ndim
  do j=i+1,ndim
    run = run+1
    EDIFF(run) = EigEn(j)-EigEn(i)
    dipnorm = DipFrintermx2(i,j)*DipFrintermx2(i,j)+DipFrintermy2(i,j)*DipFrintermy2(i,j)+DipFrintermz2(i,j)*DipFrintermz2(i,j)
    OSCSTR(run) = (Two/Three)*EDIFF(run)*dipnorm
    OSCSTR(run) = OSCSTR(run)/auToEV
    write(u6,2) i,j,EDIFF(run),DipFrintermx2(i,j),DipFrintermy2(i,j),DipFrintermz2(i,j),OSCSTR(run)
  end do
end do

call Add_Info('FRENKEL_OSCSTR',OSCSTR,nTri_Elem(ndim-1),6)

if (iPL >= 3) then
  run = 0
  write(u6,*) 'Excitonic absorption spectrum'
  write(u6,3) 'from','to','excitation energy [Ha]','Dx','Dy','Dz','osc.str.'
  do i=1,ndim
    do j=i+1,ndim
      run = run+1
      EDIFF(run) = (EigEn(j)-EigEn(i))/auToEV
      dipnorm = DipFrintermx2(i,j)*DipFrintermx2(i,j)+DipFrintermy2(i,j)*DipFrintermy2(i,j)+DipFrintermz2(i,j)*DipFrintermz2(i,j)
      OSCSTR(run) = (Two/Three)*EDIFF(run)*dipnorm
      write(u6,2) i,j,EDIFF(run),DipFrintermx2(i,j),DipFrintermy2(i,j),DipFrintermz2(i,j),OSCSTR(run)
    end do
  end do
end if

call mma_deallocate(I1)
call mma_deallocate(F1)
call mma_deallocate(D1X)
call mma_deallocate(D1Y)
call mma_deallocate(D1Z)
call mma_deallocate(I2)
call mma_deallocate(F2)
call mma_deallocate(D2X)
call mma_deallocate(D2Y)
call mma_deallocate(D2Z)
call mma_deallocate(Dip1)
call mma_deallocate(Dip2)
call mma_deallocate(EigEn)
call mma_deallocate(EigEnHa)
call mma_deallocate(EDIFF)
call mma_deallocate(OSCSTR)
call mma_deallocate(Frenkeldia)
call mma_deallocate(Freninterm)
call mma_deallocate(Freninterm2)
call mma_deallocate(Frenkelquad)
call mma_deallocate(Rassi1)
call mma_deallocate(Rassi2)
call mma_deallocate(DipFrx)
call mma_deallocate(DipFry)
call mma_deallocate(DipFrz)
call mma_deallocate(DipFrintermx)
call mma_deallocate(DipFrintermx2)
call mma_deallocate(DipFrintermy)
call mma_deallocate(DipFrintermy2)
call mma_deallocate(DipFrintermz)
call mma_deallocate(DipFrintermz2)

2 format(5X,(1x,I4,3x,I4,3x,F18.8,3x,F18.8,3x,F18.8,3x,F18.8,3x,F18.8))
3 format(6X,A,5x,A,6x,A,6x,A,19x,A,19x,A,19x,A)
222 format(5X,2(1X,I4),5X,3(1X,ES18.8))
443 format(3x,I4,3x,I4,3x,F18.8,3x,F18.8,3x,F18.8)
551 format(5X,A,4x,A,10x,A,20x,A,20x,A)

end subroutine frenkelexc
