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

subroutine G_HIGH_1(iMLTPL,dim,ESOM,GRAD,S_SOM,dipsom,Do_structure_abc,cryst,coord,gtens,maxes,iprint)
! This routine calculates the g-tensor and D-tensor in the basis of the any effective spin,
! (COMING FROM 1 MOLECULAR TERM)
!
! dim ---  the multiplicity of the effective spin

implicit none
integer, parameter :: wp = kind(0.d0)
#include "stdalloc.fh"
integer, intent(in) :: dim, iMLTPL, iprint
real(kind=8), intent(in) :: ESOM(dim), cryst(6), coord(3)
real(kind=8), intent(out) :: gtens(3), maxes(3,3)
complex(kind=8), intent(in) :: s_som(3,dim,dim), dipsom(3,dim,dim)
logical, intent(in) :: Do_structure_abc, GRAD
! local variables:
integer :: I, L, M, J, N, I1, I2, nmax, IsFreeUnit, LuDgrad, rc
real(kind=8) :: ESUM, E0, CHECK_SGN2
real(kind=8) :: knm(12,0:12)
real(kind=8), allocatable :: ELOC(:)
real(kind=8), allocatable :: axes_in_abc(:,:)
complex(kind=8), allocatable :: DIP_O(:,:), DIP_W(:,:), MUX(:,:), MUY(:,:), MUZ(:,:), MUXZ(:,:), MUZX(:,:), HZFS(:,:), &
                                DIP_MOW(:,:), HZFS_MONM(:,:), HZFS_MWNM(:,:), ZOUT(:,:), AMS(:,:,:), AMSSPIN(:,:,:), &
                                DIPSO2(:,:,:), S_SO2(:,:,:), HCF2(:,:,:,:), SP_DIPO(:), SP_DIPW(:)
complex(kind=8) :: B(3,dim,-dim:dim), C(dim,-dim:dim), BNMC(3,dim,0:dim), BNMS(3,dim,0:dim), CNMC(dim,0:dim), CNMS(dim,0:dim), &
                   ES(0:2), FS(0:2), SP_HZFSO, SP_HZFSW, SP_MOW, CHECK_SGN, m_fact, trace
external trace, IsFreeUnit
!-----------------------------------------------------------------------

call mma_allocate(ELOC,dim,'ELOC')
call mma_allocate(axes_in_abc,3,3,'axes_in_abc')

call mma_allocate(DIP_O,dim,dim,'DIP_O')
call mma_allocate(DIP_W,dim,dim,'DIP_W')
call mma_allocate(MUX,dim,dim,'MUX')
call mma_allocate(MUY,dim,dim,'MUY')
call mma_allocate(MUZ,dim,dim,'MUZ')
call mma_allocate(MUXZ,dim,dim,'MUXZ')
call mma_allocate(MUZX,dim,dim,'MUZX')
call mma_allocate(HZFS,dim,dim,'HZFS')
call mma_allocate(DIP_MOW,dim,dim,'DIP_MOW')
call mma_allocate(HZFS_MONM,dim,dim,'HZFS_MONM')
call mma_allocate(HZFS_MWNM,dim,dim,'HZFS_MWNM')
call mma_allocate(ZOUT,dim,dim,'ZOUT')
call mma_allocate(AMS,3,dim,dim,'AMS')
call mma_allocate(AMSSPIN,3,dim,dim,'AMSSPIN')
call mma_allocate(DIPSO2,3,dim,dim,'DIPSO2')
call mma_allocate(S_SO2,3,dim,dim,'S_SO2')
call mma_allocate(HCF2,dim,3,dim,dim,'HCF2')
call mma_allocate(SP_DIPO,3,'SP_DIPO')
call mma_allocate(SP_DIPW,3,'SP_DIPW')

!-----------------------------------------------------------------------
call atens(dipsom,dim,gtens,maxes,2)
! save data for construction of the blocking barriers
!do i=1,3
!  do j=1,3
!    axes(iMLTPL,i,j) = 0.0_wp
!    axes(iMLTPL,i,j) = maxes(i,j)
!  end do
!end do
! compute the magnetic axes in the crystalographic coordinate system, If requested:
if (do_structure_abc) then
  axes_in_abc = 0.0_wp
  if (iprint > 4) then
    write(6,'(A, 6F12.6)') 'cryst = ',(cryst(i),i=1,6)
    write(6,'(A, 3F12.6)') 'coord = ',(coord(i),i=1,3)
  end if

  rc = 0
  call abc_axes(cryst,coord,maxes,axes_in_abc,1,rc)

  write(6,'(19x,32a,3x,a)') '|',('-',i=1,4),'|',('-',i=1,5),' a ',('-',i=1,7),' b ',('-',i=1,7),' c ',('-',i=1,3),'|', &
                            'a , b , c  -- crystallographic axes'
  write(6,'(A,F12.9,A,3F10.6,1x,A,16x,a)') ' gX = ',gtens(1),' | Xm |',(axes_in_abc(j,1),j=1,3),'|','(defined in the input)'
  write(6,'(A,F12.9,A,3F10.6,1x,A)') ' gY = ',gtens(2),' | Ym |',(axes_in_abc(j,2),j=1,3),'|'
  write(6,'(A,F12.9,A,3F10.6,1x,A)') ' gZ = ',gtens(3),' | Zm |',(axes_in_abc(j,3),j=1,3),'|'
  write(6,'(83a)') ('-',i=1,56),'|'
end if ! do_structure_abc
! Compute the matrix elements of the magnetic moment in the coordinate system
! of magnetic axes.  ==> I.e. ROTATE the matrix DipSO to the coordinate system of magnetic axes
!maxes = 0.0_wp
!maxes(1,1) = 1.0_wp
!maxes(2,2) = 1.0_wp
!maxes(3,3) = 1.0_wp
call rotmom2(dipsom,dim,maxes,dipso2)
call rotmom2(s_som,dim,maxes,s_so2)

if (iprint > 2) then
  call prMom('G_HIGH_1:  DIPSO2(l,i,j):',dipso2,dim)
  call prMom('G_HIGH_1:   S_SO2(l,i,j):',s_so2,dim)
end if

call zcopy_(dim*dim,[(0.0_wp,0.0_wp)],0,ZOUT,1)
call zcopy_(3*dim*dim,[(0.0_wp,0.0_wp)],0,AMS,1)
call zcopy_(3*dim*dim,[(0.0_wp,0.0_wp)],0,AMSSPIN,1)
call zcopy_(dim*3*dim*dim,[(0.0_wp,0.0_wp)],0,HCF2,1)

call mu_order(dim,s_so2,dipso2,gtens,1,HCF2,AMS,AMSSPIN,ZOUT,iprint)

check_sgn = (0.0_wp,0.0_wp)
check_sgn2 = 0.0_wp
mux = (0.0_wp,0.0_wp)
muy = (0.0_wp,0.0_wp)
muz = (0.0_wp,0.0_wp)
muxz = (0.0_wp,0.0_wp)
muzx = (0.0_wp,0.0_wp)
do i=1,dim
  do j=1,dim
    mux(i,j) = HCF2(1,1,i,j)
    muy(i,j) = HCF2(1,2,i,j)
    muz(i,j) = HCF2(1,3,i,j)
  end do
end do
call ZGEMM_('N','N',dim,dim,dim,(1.0_wp,0.0_wp),mux,dim,muz,dim,(0.0_wp,0.0_wp),muxz,dim)
call ZGEMM_('N','N',dim,dim,dim,(1.0_wp,0.0_wp),muz,dim,mux,dim,(0.0_wp,0.0_wp),muzx,dim)

if (abs(muy(1,2)) > 1.d-25) then
  check_sgn = (0.0_wp,-1.0_wp)*(muxz(1,2)-muzx(1,2))/muy(1,2)
  check_sgn2 = dble(check_sgn)
else
  write(6,'(A)') 'Is it an Ising Doublet?'
  write(6,'(A)') 'For an Ising Doublet gX=gY=0, therefore, the product gX * gY * gZ is also zero'
end if

if (iprint > 2) then
  write(6,'(5x,A,2F20.14)') 'check_sgn  = ',check_sgn
  write(6,'(5x,A,F20.14)') 'check_sgn2 = ',check_sgn2
end if

write(6,'(A,F11.6)') 'CHECK-SIGN parameter = ',check_sgn2
if (check_sgn2 < 0.0_wp) then
  write(6,'(A,i2,a)') 'The sign of the product gX * gY * gZ for multiplet',iMLTPL,': < 0.'
else if (check_sgn2 > 0.0_wp) then
  write(6,'(A,i2,a,F9.6)') 'The sign of the product gX * gY * gZ for multiplet',iMLTPL,': > 0.'
end if
! Obtain the b3m and c3m coefficients:
B(1:3,1:dim,-dim:dim) = (0.0_wp,0.0_wp)
do N=1,dim-1,2
  do M=0,N
    call zcopy_(dim*dim,[(0.0_wp,0.0_wp)],0,DIP_O,1)
    call zcopy_(dim*dim,[(0.0_wp,0.0_wp)],0,DIP_W,1)

    call Stewens_matrixel(N,M,dim,DIP_O,DIP_W,IPRINT)

    if (iprint > 3) then
      write(6,*)
      write(6,'( 5x,a,i2,a,i3)') 'DIP_O, N = ',N,', M =',m
      write(6,*)
      do i=1,dim
        write(6,'(20(2F10.6,2x))') (DIP_O(i,j),j=1,dim)
      end do

      write(6,*)
      write(6,'( 5x,a,i2,a,i3)') 'DIP_W, N = ',N,', M =',m
      write(6,*)
      do i=1,dim
        write(6,'(20(2F10.6,2x))') (DIP_W(i,j),j=1,dim)
      end do
    end if

    SP_DIPO = (0.0_wp,0.0_wp)
    SP_DIPW = (0.0_wp,0.0_wp)
    SP_MOW = (0.0_wp,0.0_wp)
    SP_MOW = trace(dim,DIP_O,DIP_W)
    do l=1,3
      SP_DIPO(l) = trace(dim,AMS(l,1:dim,1:dim),DIP_O)
      SP_DIPW(l) = trace(dim,AMS(l,1:dim,1:dim),DIP_W)

      B(l,n,-m) = SP_DIPO(l)/SP_MOW
      B(l,n,m) = SP_DIPW(l)/SP_MOW
    end do ! l
  end do !m
end do !n

BNMC(1:3,1:dim,0:dim) = (0.0_wp,0.0_wp)
BNMS(1:3,1:dim,0:dim) = (0.0_wp,0.0_wp)
do n=1,dim-1,2
  do m=0,N
    do l=1,3
      if (M == 0) then
        BNMC(l,n,m) = (0.5_wp,0.0_wp)*(B(l,n,m)+B(l,n,-m))
      else
        m_fact = cmplx((-1)**M,0.0,kind=8)
        BNMC(l,n,m) = B(l,n,m)+m_fact*B(l,n,-m)
        BNMS(l,n,m) = (B(l,n,m)-m_fact*B(l,n,-m))*(0.0_wp,-1.0_wp)
      end if
    end do
  end do
end do !n

write(6,*)
write(6,'(100A)') ('-',i=1,80)
write(6,'(A)') 'DECOMPOSITION OF THE MAGNETIC MOMENT Mu_i IN IRREDUCIBLE TENSOR OPERATORS (ITO):'
write(6,'(100A)') ('-',i=1,80)
write(6,*)
write(6,'(A)') 'The quantization axis is the main magnetic axis of this multiplet (Zm).'
write(6,'(100A)') ('*',i=1,80)
write(6,'(A)') '   Mu_i = SUM_{n,m}: [ B(i,n,m) * O(n,m) +  C(i,n,m) * W(n,m) ]'
write(6,'(A)') 'where:'
write(6,'(A)') '   O(n,m) =  0.5 * ( (-1)**m * Y(n,+m) + Y(n,-m) );'
write(6,'(A)') '   W(n,m) = -0.5 * ( (-1)**m * Y(n,+m) - Y(n,-m) ) * I;   (I = imaginary unit)'
write(6,'(A)') '   n - the rank of the ITO, = 1, 3, 5, ... 2*spin;'
write(6,'(A)') '   m - the component of the ITO, = 0, 1, ... n;'
write(6,'(A)') '   i - the Cartesian projection of the magnetic moment, i = x,y,z;'
write(6,'(A)') 'These operators have been defined in: '
write(6,'(A)') '  L. F. Chibotaru, L.Ungur, J. Chem. Phys., 137, 064112 (2012).'
write(6,'(100A)') ('-',i=1,63),'|'
write(6,'(A)') '  n  |  m  | i |        B(i,n,m)       |        C(i,n,m)       |'
do N=1,dim-1,2
  write(6,'(A)') '-----|-----|---|-----------------------|-----------------------|'
  do M=0,N
    if (M /= 0) write(6,'(A)') '     |-----|---|-----------------------|-----------------------|'
    write(6,'(2(1x,I2,2x,A),1x,A,1x,A,2(ES22.14,1x,A))')N,'|',M,'|','X','|',dble(BNMC(1,N,M)),'|',dble(BNMS(1,N,M)),'|'
    write(6,'(2(1x,I2,2x,A),1x,A,1x,A,2(ES22.14,1x,A))')N,'|',M,'|','Y','|',dble(BNMC(2,N,M)),'|',dble(BNMS(2,N,M)),'|'
    write(6,'(2(1x,I2,2x,A),1x,A,1x,A,2(ES22.14,1x,A))')N,'|',M,'|','Z','|',dble(BNMC(3,N,M)),'|',dble(BNMS(3,N,M)),'|'
  end do
end do
write(6,'(100A)') ('-',i=1,63),'|'
! decomposition of the magnetic moment in Extended Stevens Operators

call Set_knm(knm)

write(6,'(/)')
write(6,'(100A)') ('*',i=1,80)
write(6,'(A)') '   Mu_i = SUM_{k,q} * [ B(i,k,q) * O(k,q) ];'
write(6,'(A)') 'where:'
write(6,'(A)') '   O(k,q) =  Extended Stevens Operators (ESO) as defined in:'
write(6,'(10x,A)') '1. Rudowicz, C.; J.Phys.C: Solid State Phys.,18(1985) 1415-1430.'
write(6,'(10x,A)') '2. Implemented in the "EasySpin" function in MATLAB, www.easyspin.org.'
write(6,'(A    )') '   k - the rank of the ITO, = 1, 3, 5, 7, 9, 11;'
write(6,'(A    )') '   q - the component of the ITO, = -k, -k+1, ... 0, 1, ... k;'
if ((dim-1) > 11) then
  write(6,'(A)') 'k = 11 may not be the highest rank of the ITO for this case, but it '
  write(6,'(A)') 'is the maximal k implemented in the "EasySpin" function in MATLAB.'
end if
write(6,'(A)') 'Knm are proportionality coefficients between the ESO and operators defined in '
write(6,'(A)') 'J. Chem. Phys., 137, 064112 (2012).'
write(6,'(100A)') ('-',i=1,51),'|'
write(6,'(A)') '  k |  q  | i |   (Knm)^2  |         B(k,q)        |'

if ((dim-1) > 11) then
  Nmax = 11
else
  nmax = dim-1
end if

do N=1,nmax,2
  write(6,'(A)') '----|-----|---|------------|-----------------------|'
  do M=-N,N
    if (M < 0) then
      write(6,'((1x,I2,1x,A),(1x,I3,1x,A),1x,A,1x,A,12x,A,2(ES22.14,1x,A))') N,'|',M,'|','X','|','|', &
                                                                             dble(BNMS(1,N,abs(M)))*knm(n,abs(m)),'|'
      write(6,'((1x,I2,1x,A),(1x,I3,1x,A),1x,A,1x,A,1x,F10.3,1x,A,2(ES22.14,1x,A))') N,'|',M,'|','Y','|', &
                                                                                     knm(n,abs(m))*knm(n,abs(m)),'|', &
                                                                                     dble(BNMS(2,N,abs(M)))*knm(n,abs(m)),'|'
      write(6,'((1x,I2,1x,A),(1x,I3,1x,A),1x,A,1x,A,12x,A,2(ES22.14,1x,A))') N,'|',M,'|','Z','|','|', &
                                                                             dble(BNMS(3,N,abs(M)))*knm(n,abs(m)),'|'
      write(6,'(A)') '    |-----|---|------------|-----------------------|'
    else
      write(6,'((1x,I2,1x,A),(1x,I3,1x,A),1x,A,1x,A,12x,A,2(ES22.14,1x,A))') N,'|',M,'|','X','|','|', &
                                                                             dble(BNMC(1,N,abs(M)))*knm(n,abs(m)),'|'
      write(6,'((1x,I2,1x,A),(1x,I3,1x,A),1x,A,1x,A,1x,F10.3,1x,A,2(ES22.14,1x,A))') N,'|',M,'|','Y','|', &
                                                                                     knm(n,abs(m))*knm(n,abs(m)),'|', &
                                                                                     dble(BNMC(2,N,abs(M)))*knm(n,abs(m)),'|'
      write(6,'((1x,I2,1x,A),(1x,I3,1x,A),1x,A,1x,A,12x,A,2(ES22.14,1x,A))') N,'|',M,'|','Z','|','|', &
                                                                             dble(BNMC(3,N,abs(M)))*knm(n,abs(m)),'|'
      if (M /= N) write(6,'(A)') '    |-----|---|------------|-----------------------|'
    end if !M<0
  end do !M
end do !N
write(6,'(100A)') ('-',i=1,51),'|'

! Calculation of the ZFS tensors and the coefficients of the higher order spin-operators Enm and Fnm
if (dim > 2) then
  ESUM = 0.0_wp
  do I=1,dim
    ELOC(i) = 0.0_wp
    ESUM = ESUM+ESOM(I)
  end do
  E0 = ESUM/dble(dim)
  do I=1,dim
    ELOC(I) = ESOM(I)-E0
  end do

  call zcopy_(dim*dim,[(0.0_wp,0.0_wp)],0,HZFS,1)
  do i1=1,dim
    do i2=1,dim
      do i=1,dim
        HZFS(i1,i2) = HZFS(i1,i2)+ELOC(i)*conjg(ZOUT(i,I1))*ZOUT(i,I2)
      end do
    end do
  end do

  if (iprint > 2) then
    write(6,*)
    write(6,'(5X,A)') 'SPIN-ORBIT ENERGIES OF THE FIRST MOLECULAR TERM SHIfTED TO THE MASS CENTER'
    write(6,*)
    write(6,'(5X,A,F10.4)') 'E0 = ',E0
    write(6,'(15X,A,11X,A)') 'ESOM','ESO_LOC'

    do I=1,dim
      write(6,'(5X,F15.6,2X,F15.6)') ESOM(I),ELOC(I)
    end do
  end if

  write(6,*)
  write(6,'(100A)') ('-',i=1,87)
  write(6,'(A)') 'DECOMPOSITION OF THE ZERO-FIELD SPLITTING (ZFS) IN IRREDUCIBLE TENSOR OPERATORS (ITO):'
  write(6,'(100A)') ('-',i=1,87)
  write(6,*)
  write(6,'(A)') 'Ab Initio Calculated Zero-Field Splitting Matrix written in the basis of Pseudospin Eigenfunctions'
  if (mod(dim,2) == 0) then
    write(6,'(950A)') ('-',i=1,10),(('-',i=1,24),j=1,dim),'|'
    write(6,'(10x,A,50(8x,A,I3,A,7x,A))') '|',('|',2*i-dim-1,'/2 >','|',i=1,dim)
    write(6,'(950A)') ('-',i=1,10),'|',(('-',i=1,23),'|',j=1,dim)
    do i=1,dim
      write(6,'(1x,A,I3,A,1x,A,50(2F11.5,1x,A))') '<',2*i-dim-1,'/2','| |',(HZFS(j,i),'|',j=1,dim)
    end do
    write(6,'(950A)') ('-',i=1,10),(('-',i=1,24),j=1,dim),'|'
  else
    write(6,'(950A)') ('-',i=1,8),(('-',i=1,24),j=1,dim),'|'
    write(6,'(8x,A,50(8x,A,I3,A,9x,A))') '|',('|',-(dim-1)/2-1+i,' >','|',i=1,dim)
    write(6,'(950A)') ('-',i=1,8),'|',(('-',i=1,23),'|',j=1,dim)
    do I=1,dim
      write(6,'(1x,A,I3,1x,A,50(2F11.5,1x,A))') '<',-(dim-1)/2-1+i,'| |',(HZFS(j,i),'|',j=1,dim)
    end do
    write(6,'(950A)') ('-',i=1,8),(('-',i=1,24),j=1,dim),'|'
  end if

  C(1:dim,-dim:dim) = (0.0_wp,0.0_wp)
  do N=2,dim-1,2
    do M=0,N
      call zcopy_(dim*dim,[(0.0_wp,0.0_wp)],0,DIP_O,1)
      call zcopy_(dim*dim,[(0.0_wp,0.0_wp)],0,DIP_W,1)
      call zcopy_(dim*dim,[(0.0_wp,0.0_wp)],0,HZFS_MONM,1)
      call zcopy_(dim*dim,[(0.0_wp,0.0_wp)],0,HZFS_MWNM,1)
      call zcopy_(dim*dim,[(0.0_wp,0.0_wp)],0,DIP_MOW,1)

      call Stewens_matrixel(N,M,dim,DIP_O,DIP_W,IPRINT)

      if (IPRINT > 5) then
        write(6,'(/)')
        write(6,'(5x,a,i3,3x,A,I3)') 'DIP_STEWENS_O  N = ',N,'M =',M
        write(6,*)
        do i=1,dim
          write(6,'(20(2X,2ES20.10))') (DIP_O(i,j),j=1,dim)
        end do
        write(6,'(/)')
        write(6,'(5x,a,i3,3x,A,I3)') 'DIP_STEWENS_W  N = ',N,'M =',M
        write(6,*)
        do i=1,dim
          write(6,'(20(2X,2ES20.10))') (DIP_W(i,j),j=1,dim)
        end do
      end if

      SP_HZFSO = (0.0_wp,0.0_wp)
      SP_HZFSW = (0.0_wp,0.0_wp)
      SP_MOW = (0.0_wp,0.0_wp)
      SP_MOW = trace(dim,DIP_O,DIP_W)
      SP_HZFSO = trace(dim,HZFS,DIP_O)
      SP_HZFSW = trace(dim,HZFS,DIP_W)

      C(N,-M) = SP_HZFSO/SP_MOW
      C(N,M) = SP_HZFSW/SP_MOW

      if (IPRINT > 5) then
        write(6,'(/)')
        write(6,'( 5x,a)') 'HZFS_MONM(i,j)'
        write(6,*)
        do i=1,dim
          write(6,'(20(2F18.10,2x))') (HZFS_MONM(i,j),j=1,dim)
        end do
        write(6,*)
        write(6,'(5X,a,2F18.10)') 'SP_HZFSO = ',SP_HZFSO
        write(6,'(/)')
        write(6,'( 5x,a)') 'HZFS_MWNM(i,j)'
        write(6,*)
        do i=1,dim
          write(6,'(20(2F18.10,2x))') (HZFS_MWNM(i,j),j=1,dim)
        end do
        write(6,*)
        write(6,'(5X,a,2F18.10)') 'SP_HZFSW = ',SP_HZFSW
        write(6,'(/)')
        write(6,'( 5x,a)') 'HZFS_MOW(i,j) (i,j)'
        write(6,*)
        do i=1,dim
          write(6,'(20(2F18.10,2x))') (DIP_MOW(i,j),j=1,dim)
        end do
        write(6,*)
        write(6,'(5X,a,2F18.10)') 'SP_MOW = ',SP_MOW
      end if
    end do !M
  end do !N

  CNMC(1:dim,0:dim) = (0.0_wp,0.0_wp)
  CNMS(1:dim,0:dim) = (0.0_wp,0.0_wp)
  do N=2,dim-1,2
    do M=0,N
      if (M == 0) then
        CNMC(N,M) = (0.5_wp,0.0_wp)*(C(N,M)+C(N,-M))
      else
        m_fact = cmplx((-1)**M,0,wp)
        CNMC(N,M) = C(N,M)+m_fact*C(N,-M)
        CNMS(N,M) = (C(N,M)-m_fact*C(N,-M))*(0.0_wp,-1.0_wp)
      end if
    end do
  end do

  write(6,'(A)') 'The ZFS Hamiltonian:'
  write(6,'(A)') '   ZFS = SUM_{n,m}: [ E(n,m) * O(n,m) +  F(n,m) * W(n,m) ]'
  write(6,'(A)') 'where:'
  write(6,'(A)') '   O(n,m) =  0.5 * ( (-1)**m * Y(n,+m) + Y(n,-m) );'
  write(6,'(A)') '   W(n,m) = -0.5 * ( (-1)**m * Y(n,+m) - Y(n,-m) ) * I;    (I = imaginary unit)'
  write(6,'(A)') '   n - the rank of the ITO, = 2, 4, 6, ... 2*spin;'
  write(6,'(A)') '   m - the component of the ITO, = 0, 1, ... n;'
  write(6,'(A)') 'The quantization axis is the main magnetic axis of this multiplet (Zm).'
  write(6,'(100A)') ('-',i=1,59),'|'
  write(6,'(A)') '  n  |  m  |         E(n,m)        |         F(n,m)        |'
  do N=2,dim-1,2
    write(6,'(A)') '-----|-----|-----------------------|-----------------------|'
    do M=0,N
      write(6,'(2(1x,I2,2x,A),2(ES22.14,1x,A))') N,'|',M,'|',dble(CNMC(N,M)),'|',dble(CNMS(N,M)),'|'
    end do
  end do
  write(6,'(100A)') ('-',i=1,59),'|'

! decomposition of the ZFS matrix in ExtEnded Stevens Operators
  write(6,'(100A)') ('*',i=1,80)
  write(6,'(A)') 'The ZFS Hamiltonian:'
  write(6,'(A)') '   ZFS = SUM_{k,q} * [ B(k,q) * O(k,q) ];'
  write(6,'(A)') 'where:'
  write(6,'(A)') '   O(k,q) =  Extended Stevens Operators (ESO) as defined in:'
  write(6,'(10x,A)') '1. Rudowicz, C.; J.Phys.C: Solid State Phys.,18(1985) 1415-1430.'
  write(6,'(10x,A)') '2. Implemented in the "EasySpin" function in MATLAB, www.easyspin.org.'
  write(6,'(A    )') '   k - the rank of the ITO, = 2, 4, 6, 8, 10, 12.'
  write(6,'(A)') '   q - the component of the ITO, = -k, -k+1, ... 0, 1, ... k;'

  if ((dim-1) > 12) then
    write(6,'(A)') 'k = 12 may not be the highest rank of the ITO for this case, but it '
    write(6,'(A)') 'is the maximal k implemented in the "EasySpin" function in MATLAB.'
  end if

  write(6,'(A)') 'Knm are proportionality coefficients between the ESO and operators defined in '
  write(6,'(A)') 'J. Chem. Phys., 137, 064112 (2012).'
  write(6,'(100A)') ('-',i=1,48),'|'
  write(6,'(A)') '  k |  q  |    (Knm)^2  |         B(k,q)        |'
  if ((dim-1) > 12) then
    Nmax = 12
  else
    Nmax = dim-1
  end if
  do N=2,Nmax,2
    write(6,'(A)') '----|-----|-------------|-----------------------|'
    do M=-N,N
      if (M < 0) then
        write(6,'((1x,I2,1x,A),(1x,I3,1x,A),F11.2,2x,A,2(ES22.14,1x,A))') N,'|',M,'|',knm(n,abs(m))*knm(n,abs(m)),'|', &
                                                                          dble(CNMS(N,abs(M)))*knm(n,abs(m)),'|'
      else
        write(6,'((1x,I2,1x,A),(1x,I3,1x,A),F11.2,2x,A,2(ES22.14,1x,A))') N,'|',M,'|',knm(n,abs(m))*knm(n,abs(m)),'|', &
                                                                          dble(CNMC(N,abs(M)))*Knm(n,abs(m)),'|'
      end if
    end do
  end do
  write(6,'(100A)') ('-',i=1,48),'|'
  !-----------------------------
  ! for the interface related to CF gradient calculation:
  if (GRAD) then
    LuDgrad = IsFreeUnit(81)
    call molcas_open(LuDgrad,'DMAT')
    do N=2,Nmax,2
      write(6,'(A)') '----|-----|-------------|-----------------------|'
      do M=-N,N
        if (M < 0) then
          write(LuDgrad,'(I4,I4,1x,2(ES25.15))') N,M,dble(CNMS(n,abs(n)))*knm(n,abs(m))
        else
          write(LuDgrad,'(I4,I4,1x,2(ES25.15))') N,M,dble(CNMC(n,abs(m)))*Knm(n,abs(m))
        end if
      end do
    end do
    close(LuDgrad)
  end if
  !-----------------------------
!
  do i=0,2
    ES(i) = 0.0_wp
    FS(i) = 0.0_wp
    ES(i) = CNMC(2,i)
    FS(i) = CNMS(2,i)
  end do
  call DMATRIX(Es,Fs,maxes,2)
end if ! decomposition of ZFS in higher-order ITO operators

write(6,*)
write(6,'(A,I2,A,I3,A)') 'ANGULAR MOMENTS ALONG THE MAIN MAGNETIC AXES'
call moments(dim,s_so2,dipso2,iprint)

!----------------------------------------------------------------------
call mma_deallocate(ELOC)
call mma_deallocate(axes_in_abc)
call mma_deallocate(DIP_O)
call mma_deallocate(DIP_W)
call mma_deallocate(MUX)
call mma_deallocate(MUY)
call mma_deallocate(MUZ)
call mma_deallocate(MUXZ)
call mma_deallocate(MUZX)
call mma_deallocate(HZFS)
call mma_deallocate(DIP_MOW)
call mma_deallocate(HZFS_MONM)
call mma_deallocate(HZFS_MWNM)
call mma_deallocate(ZOUT)
call mma_deallocate(AMS)
call mma_deallocate(AMSSPIN)
call mma_deallocate(DIPSO2)
call mma_deallocate(S_SO2)
call mma_deallocate(HCF2)
call mma_deallocate(SP_DIPO)
call mma_deallocate(SP_DIPW)

return

end subroutine G_HIGH_1
