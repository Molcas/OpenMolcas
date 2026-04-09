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

subroutine G_HIGH_1(iMLTPL,d,ESOM,GRAD,S_SOM,dipsom,Do_structure_abc,cryst,coord,gtens,maxes,iprint)
! This routine calculates the g-tensor and D-tensor in the basis of the any effective spin,
! (COMING FROM 1 MOLECULAR TERM)
!
! d ---  the multiplicity of the effective spin

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half, cZero, cOne, Onei
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iMLTPL, d, iprint
real(kind=wp), intent(in) :: ESOM(d), cryst(6), coord(3)
logical(kind=iwp), intent(in) :: GRAD, Do_structure_abc
complex(kind=wp), intent(in) :: s_som(3,d,d), dipsom(3,d,d)
real(kind=wp), intent(out) :: gtens(3), maxes(3,3)
integer(kind=iwp) :: I, I1, I2, J, L, LuDgrad, M, N, nmax, rc
real(kind=wp) :: axes_in_abc(3,3), CHECK_SGN2, E0, ESUM, knm(12,0:12), m_fact
complex(kind=wp) :: CHECK_SGN, ES(0:2), FS(0:2), SP_HZFSO, SP_HZFSW, SP_MOW
real(kind=wp), allocatable :: ELOC(:)
complex(kind=wp), allocatable :: AMS(:,:,:), AMS_TMP(:,:), AMSSPIN(:,:,:), B(:,:,:), BNMC(:,:,:), BNMS(:,:,:), C(:,:), CNMC(:,:), &
                                 CNMS(:,:), DIP_MOW(:,:), DIP_O(:,:), DIP_W(:,:), DIPSO2(:,:,:), HCF2(:,:,:,:), HZFS(:,:), &
                                 HZFS_MONM(:,:), HZFS_MWNM(:,:), MUX(:,:), MUXZ(:,:), MUY(:,:), MUZ(:,:), MUZX(:,:), S_SO2(:,:,:), &
                                 SP_DIPO(:), SP_DIPW(:), ZOUT(:,:)
integer(kind=iwp), external :: IsFreeUnit
complex(kind=wp), external :: trace
!-----------------------------------------------------------------------

call mma_allocate(ELOC,d,'ELOC')

call mma_allocate(DIP_O,d,d,'DIP_O')
call mma_allocate(DIP_W,d,d,'DIP_W')
call mma_allocate(MUX,d,d,'MUX')
call mma_allocate(MUY,d,d,'MUY')
call mma_allocate(MUZ,d,d,'MUZ')
call mma_allocate(MUXZ,d,d,'MUXZ')
call mma_allocate(MUZX,d,d,'MUZX')
call mma_allocate(HZFS,d,d,'HZFS')
call mma_allocate(DIP_MOW,d,d,'DIP_MOW')
call mma_allocate(HZFS_MONM,d,d,'HZFS_MONM')
call mma_allocate(HZFS_MWNM,d,d,'HZFS_MWNM')
call mma_allocate(ZOUT,d,d,'ZOUT')
call mma_allocate(AMS,3,d,d,'AMS')
call mma_allocate(AMS_TMP,d,d,'AMS_TMP')
call mma_allocate(AMSSPIN,3,d,d,'AMSSPIN')
call mma_allocate(DIPSO2,3,d,d,'DIPSO2')
call mma_allocate(S_SO2,3,d,d,'S_SO2')
call mma_allocate(HCF2,d,3,d,d,'HCF2')
call mma_allocate(SP_DIPO,3,'SP_DIPO')
call mma_allocate(SP_DIPW,3,'SP_DIPW')

!-----------------------------------------------------------------------
call atens(dipsom,d,gtens,maxes,2)
! save data for construction of the blocking barriers
!axes(iMLTPL,:,:) = maxes(:,:)
! compute the magnetic axes in the crystalographic coordinate system, If requested:
if (do_structure_abc) then
  axes_in_abc(:,:) = Zero
  if (iprint > 4) then
    write(u6,'(A, 6F12.6)') 'cryst = ',(cryst(i),i=1,6)
    write(u6,'(A, 3F12.6)') 'coord = ',(coord(i),i=1,3)
  end if

  rc = 0
  call abc_axes(cryst,coord,maxes,axes_in_abc,1,rc)

  write(u6,'(19x,11a,3x,a)') '|',repeat('-',4),'|',repeat('-',5),' a ',repeat('-',7),' b ',repeat('-',7),' c ',repeat('-',3),'|', &
                             'a , b , c  -- crystallographic axes'
  write(u6,'(A,F12.9,A,3F10.6,1x,A,16x,a)') ' gX = ',gtens(1),' | Xm |',(axes_in_abc(j,1),j=1,3),'|','(defined in the input)'
  write(u6,'(A,F12.9,A,3F10.6,1x,A)') ' gY = ',gtens(2),' | Ym |',(axes_in_abc(j,2),j=1,3),'|'
  write(u6,'(A,F12.9,A,3F10.6,1x,A)') ' gZ = ',gtens(3),' | Zm |',(axes_in_abc(j,3),j=1,3),'|'
  write(u6,'(2a)') repeat('-',56),'|'
end if ! do_structure_abc
! Compute the matrix elements of the magnetic moment in the coordinate system
! of magnetic axes.  ==> I.e. ROTATE the matrix DipSO to the coordinate system of magnetic axes
!call unitmat(maxes,3)
call rotmom2(dipsom,d,maxes,dipso2)
call rotmom2(s_som,d,maxes,s_so2)

if (iprint > 2) then
  call prMom('G_HIGH_1:  DIPSO2(l,i,j):',dipso2,d)
  call prMom('G_HIGH_1:   S_SO2(l,i,j):',s_so2,d)
end if

ZOUT(:,:) = cZero
AMS(:,:,:) = cZero
AMSSPIN(:,:,:) = cZero
HCF2(:,:,:,:) = cZero

call mu_order(d,s_so2,dipso2,gtens,1,HCF2,AMS,AMSSPIN,ZOUT,iprint)

check_sgn = cZero
check_sgn2 = Zero
mux(:,:) = HCF2(1,1,:,:)
muy(:,:) = HCF2(1,2,:,:)
muz(:,:) = HCF2(1,3,:,:)
call ZGEMM_('N','N',d,d,d,cOne,mux,d,muz,d,cZero,muxz,d)
call ZGEMM_('N','N',d,d,d,cOne,muz,d,mux,d,cZero,muzx,d)

if (abs(muy(1,2)) > 1.0e-25_wp) then
  check_sgn = -Onei*(muxz(1,2)-muzx(1,2))/muy(1,2)
  check_sgn2 = real(check_sgn)
else
  write(u6,'(A)') 'Is it an Ising doublet?'
  write(u6,'(A)') 'For an Ising doublet gX=gY=0, therefore, the product gX * gY * gZ is also zero'
end if

if (iprint > 2) then
  write(u6,'(5x,A,2F20.14)') 'check_sgn  = ',check_sgn
  write(u6,'(5x,A,F20.14)') 'check_sgn2 = ',check_sgn2
end if

write(u6,'(A,F11.6)') 'CHECK-SIGN parameter = ',check_sgn2
if (check_sgn2 < Zero) then
  write(u6,'(A,i2,a)') 'The sign of the product gX * gY * gZ for multiplet',iMLTPL,': < 0.'
else if (check_sgn2 > Zero) then
  write(u6,'(A,i2,a,F9.6)') 'The sign of the product gX * gY * gZ for multiplet',iMLTPL,': > 0.'
end if
! Obtain the b3m and c3m coefficients:
call mma_allocate(B,[1,3],[1,d],[-d,d],label='B')
call mma_allocate(BNMC,[1,3],[1,d],[0,d],label='BNMC')
call mma_allocate(BNMS,[1,3],[1,d],[0,d],label='BNMS')
B(:,:,:) = cZero
do N=1,d-1,2
  do M=0,N
    DIP_O(:,:) = cZero
    DIP_W(:,:) = cZero

    call Stewens_matrixel(N,M,d,DIP_O,DIP_W,IPRINT)

    if (iprint > 3) then
      write(u6,*)
      write(u6,'( 5x,a,i2,a,i3)') 'DIP_O, N = ',N,', M =',m
      write(u6,*)
      do i=1,d
        write(u6,'(20(2F10.6,2x))') (DIP_O(i,j),j=1,d)
      end do

      write(u6,*)
      write(u6,'( 5x,a,i2,a,i3)') 'DIP_W, N = ',N,', M =',m
      write(u6,*)
      do i=1,d
        write(u6,'(20(2F10.6,2x))') (DIP_W(i,j),j=1,d)
      end do
    end if

    SP_DIPO(:) = cZero
    SP_DIPW(:) = cZero
    SP_MOW = cZero
    SP_MOW = trace(d,DIP_O,DIP_W)
    do l=1,3
      AMS_TMP(:,:) = AMS(l,:,:)
      SP_DIPO(l) = trace(d,AMS_TMP,DIP_O)
      SP_DIPW(l) = trace(d,AMS_TMP,DIP_W)

      B(l,n,-m) = SP_DIPO(l)/SP_MOW
      B(l,n,m) = SP_DIPW(l)/SP_MOW
    end do ! l
  end do !m
end do !n

BNMC(:,:,:) = cZero
BNMS(:,:,:) = cZero
do n=1,d-1,2
  BNMC(:,n,0) = B(:,n,0)
  do m=1,N
    m_fact = (-One)**M
    BNMC(:,n,m) = B(:,n,m)+m_fact*B(:,n,-m)
    BNMS(:,n,m) = -Onei*(B(:,n,m)-m_fact*B(:,n,-m))
  end do
end do !n

write(u6,*)
write(u6,'(A)') repeat('-',80)
write(u6,'(A)') 'DECOMPOSITION OF THE MAGNETIC MOMENT Mu_i IN IRREDUCIBLE TENSOR OPERATORS (ITO):'
write(u6,'(A)') repeat('-',80)
write(u6,*)
write(u6,'(A)') 'The quantization axis is the main magnetic axis of this multiplet (Zm).'
write(u6,'(A)') repeat('*',80)
write(u6,'(A)') '   Mu_i = SUM_{n,m}: [ B(i,n,m) * O(n,m) +  C(i,n,m) * W(n,m) ]'
write(u6,'(A)') 'where:'
write(u6,'(A)') '   O(n,m) =  0.5 * ( (-1)**m * Y(n,+m) + Y(n,-m) );'
write(u6,'(A)') '   W(n,m) = -0.5 * ( (-1)**m * Y(n,+m) - Y(n,-m) ) * I;   (I = imaginary unit)'
write(u6,'(A)') '   n - the rank of the ITO, = 1, 3, 5, ... 2*spin;'
write(u6,'(A)') '   m - the component of the ITO, = 0, 1, ... n;'
write(u6,'(A)') '   i - the Cartesian projection of the magnetic moment, i = x,y,z;'
write(u6,'(A)') 'These operators have been defined in: '
write(u6,'(A)') '  L. F. Chibotaru, L.Ungur, J. Chem. Phys., 137, 064112 (2012).'
write(u6,'(2A)') repeat('-',63),'|'
write(u6,'(A)') '  n  |  m  | i |        B(i,n,m)       |        C(i,n,m)       |'
do N=1,d-1,2
  write(u6,'(A)') '-----|-----|---|-----------------------|-----------------------|'
  do M=0,N
    if (M /= 0) write(u6,'(A)') '     |-----|---|-----------------------|-----------------------|'
    write(u6,'(2(1x,I2,2x,A),1x,A,1x,A,2(ES22.14,1x,A))') N,'|',M,'|','X','|',real(BNMC(1,N,M)),'|',real(BNMS(1,N,M)),'|'
    write(u6,'(2(1x,I2,2x,A),1x,A,1x,A,2(ES22.14,1x,A))') N,'|',M,'|','Y','|',real(BNMC(2,N,M)),'|',real(BNMS(2,N,M)),'|'
    write(u6,'(2(1x,I2,2x,A),1x,A,1x,A,2(ES22.14,1x,A))') N,'|',M,'|','Z','|',real(BNMC(3,N,M)),'|',real(BNMS(3,N,M)),'|'
  end do
end do
write(u6,'(2A)') repeat('-',63),'|'
! decomposition of the magnetic moment in Extended Stevens Operators

call Set_knm(knm)

write(u6,'(/)')
write(u6,'(A)') repeat('*',80)
write(u6,'(A)') '   Mu_i = SUM_{k,q} * [ B(i,k,q) * O(k,q) ];'
write(u6,'(A)') 'where:'
write(u6,'(A)') '   O(k,q) =  Extended Stevens Operators (ESO) as defined in:'
write(u6,'(10x,A)') '1. Rudowicz, C.; J.Phys.C: Solid State Phys.,18(1985) 1415-1430.'
write(u6,'(10x,A)') '2. Implemented in the "EasySpin" function in MATLAB, www.easyspin.org.'
write(u6,'(A    )') '   k - the rank of the ITO, = 1, 3, 5, 7, 9, 11;'
write(u6,'(A    )') '   q - the component of the ITO, = -k, -k+1, ... 0, 1, ... k;'
if (d-1 > 11) then
  write(u6,'(A)') 'k = 11 may not be the highest rank of the ITO for this case, but it '
  write(u6,'(A)') 'is the maximal k implemented in the "EasySpin" function in MATLAB.'
end if
write(u6,'(A)') 'Knm are proportionality coefficients between the ESO and operators defined in '
write(u6,'(A)') 'J. Chem. Phys., 137, 064112 (2012).'
write(u6,'(2A)') repeat('-',51),'|'
write(u6,'(A)') '  k |  q  | i |   (Knm)^2  |         B(k,q)        |'

if (d-1 > 11) then
  Nmax = 11
else
  nmax = d-1
end if

do N=1,nmax,2
  write(u6,'(A)') '----|-----|---|------------|-----------------------|'
  do M=-N,N
    if (M < 0) then
      write(u6,'((1x,I2,1x,A),(1x,I3,1x,A),1x,A,1x,A,12x,A,2(ES22.14,1x,A))') N,'|',M,'|','X','|','|', &
                                                                              real(BNMS(1,N,abs(M)))*knm(n,abs(m)),'|'
      write(u6,'((1x,I2,1x,A),(1x,I3,1x,A),1x,A,1x,A,1x,F10.3,1x,A,2(ES22.14,1x,A))') N,'|',M,'|','Y','|', &
                                                                                      knm(n,abs(m))*knm(n,abs(m)),'|', &
                                                                                      real(BNMS(2,N,abs(M)))*knm(n,abs(m)),'|'
      write(u6,'((1x,I2,1x,A),(1x,I3,1x,A),1x,A,1x,A,12x,A,2(ES22.14,1x,A))') N,'|',M,'|','Z','|','|', &
                                                                              real(BNMS(3,N,abs(M)))*knm(n,abs(m)),'|'
      write(u6,'(A)') '    |-----|---|------------|-----------------------|'
    else
      write(u6,'((1x,I2,1x,A),(1x,I3,1x,A),1x,A,1x,A,12x,A,2(ES22.14,1x,A))') N,'|',M,'|','X','|','|', &
                                                                              real(BNMC(1,N,abs(M)))*knm(n,abs(m)),'|'
      write(u6,'((1x,I2,1x,A),(1x,I3,1x,A),1x,A,1x,A,1x,F10.3,1x,A,2(ES22.14,1x,A))') N,'|',M,'|','Y','|', &
                                                                                      knm(n,abs(m))*knm(n,abs(m)),'|', &
                                                                                      real(BNMC(2,N,abs(M)))*knm(n,abs(m)),'|'
      write(u6,'((1x,I2,1x,A),(1x,I3,1x,A),1x,A,1x,A,12x,A,2(ES22.14,1x,A))') N,'|',M,'|','Z','|','|', &
                                                                              real(BNMC(3,N,abs(M)))*knm(n,abs(m)),'|'
      if (M /= N) write(u6,'(A)') '    |-----|---|------------|-----------------------|'
    end if !M<0
  end do !M
end do !N
write(u6,'(2A)') repeat('-',51),'|'

call mma_deallocate(B)
call mma_deallocate(BNMC)
call mma_deallocate(BNMS)

! Calculation of the ZFS tensors and the coefficients of the higher order spin-operators Enm and Fnm
if (d > 2) then
  ESUM = sum(ESOM(:))
  E0 = ESUM/real(d,kind=wp)
  ELOC(:) = ESOM(:)-E0

  HZFS(:,:) = cZero
  do i1=1,d
    do i2=1,d
      do i=1,d
        HZFS(i1,i2) = HZFS(i1,i2)+ELOC(i)*conjg(ZOUT(i,I1))*ZOUT(i,I2)
      end do
    end do
  end do

  if (iprint > 2) then
    write(u6,*)
    write(u6,'(5X,A)') 'SPIN-ORBIT ENERGIES OF THE FIRST MOLECULAR TERM SHIfTED TO THE MASS CENTER'
    write(u6,*)
    write(u6,'(5X,A,F10.4)') 'E0 = ',E0
    write(u6,'(15X,A,11X,A)') 'ESOM','ESO_LOC'

    do I=1,d
      write(u6,'(5X,F15.6,2X,F15.6)') ESOM(I),ELOC(I)
    end do
  end if

  write(u6,*)
  write(u6,'(A)') repeat('-',87)
  write(u6,'(A)') 'DECOMPOSITION OF THE ZERO-FIELD SPLITTING (ZFS) IN IRREDUCIBLE TENSOR OPERATORS (ITO):'
  write(u6,'(A)') repeat('-',87)
  write(u6,*)
  write(u6,'(A)') 'Ab Initio Calculated Zero-Field Splitting Matrix written in the basis of Pseudospin Eigenfunctions'
  if (mod(d,2) == 0) then
    write(u6,'(52A)') repeat('-',10),(repeat('-',24),j=1,d),'|'
    write(u6,'(10x,A,50(8x,A,I3,A,7x,A))') '|',('|',2*i-d-1,'/2 >','|',i=1,d)
    write(u6,'(102A)') repeat('-',10),'|',(repeat('-',23),'|',j=1,d)
    do i=1,d
      write(u6,'(1x,A,I3,A,1x,A,50(2F11.5,1x,A))') '<',2*i-d-1,'/2','| |',(HZFS(j,i),'|',j=1,d)
    end do
    write(u6,'(52A)') repeat('-',10),(repeat('-',24),j=1,d),'|'
  else
    write(u6,'(52A)') repeat('-',8),(repeat('-',24),j=1,d),'|'
    write(u6,'(8x,A,50(8x,A,I3,A,9x,A))') '|',('|',-(d-1)/2-1+i,' >','|',i=1,d)
    write(u6,'(102A)') repeat('-',8),'|',(repeat('-',23),'|',j=1,d)
    do I=1,d
      write(u6,'(1x,A,I3,1x,A,50(2F11.5,1x,A))') '<',-(d-1)/2-1+i,'| |',(HZFS(j,i),'|',j=1,d)
    end do
    write(u6,'(52A)') repeat('-',8),(repeat('-',24),j=1,d),'|'
  end if

  call mma_allocate(C,[1,d],[-d,d],label='C')
  call mma_allocate(CNMC,[1,d],[0,d],label='CNMC')
  call mma_allocate(CNMS,[1,d],[0,d],label='CNMS')
  C(:,:) = cZero
  do N=2,d-1,2
    do M=0,N
      DIP_O(:,:) = cZero
      DIP_W(:,:) = cZero
      HZFS_MONM(:,:) = cZero
      HZFS_MWNM(:,:) = cZero
      DIP_MOW(:,:) = cZero

      call Stewens_matrixel(N,M,d,DIP_O,DIP_W,IPRINT)

      if (IPRINT > 5) then
        write(u6,'(/)')
        write(u6,'(5x,a,i3,3x,A,I3)') 'DIP_STEWENS_O  N = ',N,'M =',M
        write(u6,*)
        do i=1,d
          write(u6,'(20(2X,2ES20.10))') (DIP_O(i,j),j=1,d)
        end do
        write(u6,'(/)')
        write(u6,'(5x,a,i3,3x,A,I3)') 'DIP_STEWENS_W  N = ',N,'M =',M
        write(u6,*)
        do i=1,d
          write(u6,'(20(2X,2ES20.10))') (DIP_W(i,j),j=1,d)
        end do
      end if

      SP_HZFSO = cZero
      SP_HZFSW = cZero
      SP_MOW = cZero
      SP_MOW = trace(d,DIP_O,DIP_W)
      SP_HZFSO = trace(d,HZFS,DIP_O)
      SP_HZFSW = trace(d,HZFS,DIP_W)

      C(N,-M) = SP_HZFSO/SP_MOW
      C(N,M) = SP_HZFSW/SP_MOW

      if (IPRINT > 5) then
        write(u6,'(/)')
        write(u6,'( 5x,a)') 'HZFS_MONM(i,j)'
        write(u6,*)
        do i=1,d
          write(u6,'(20(2F18.10,2x))') (HZFS_MONM(i,j),j=1,d)
        end do
        write(u6,*)
        write(u6,'(5X,a,2F18.10)') 'SP_HZFSO = ',SP_HZFSO
        write(u6,'(/)')
        write(u6,'( 5x,a)') 'HZFS_MWNM(i,j)'
        write(u6,*)
        do i=1,d
          write(u6,'(20(2F18.10,2x))') (HZFS_MWNM(i,j),j=1,d)
        end do
        write(u6,*)
        write(u6,'(5X,a,2F18.10)') 'SP_HZFSW = ',SP_HZFSW
        write(u6,'(/)')
        write(u6,'( 5x,a)') 'HZFS_MOW(i,j) (i,j)'
        write(u6,*)
        do i=1,d
          write(u6,'(20(2F18.10,2x))') (DIP_MOW(i,j),j=1,d)
        end do
        write(u6,*)
        write(u6,'(5X,a,2F18.10)') 'SP_MOW = ',SP_MOW
      end if
    end do !M
  end do !N

  CNMC(:,:) = cZero
  CNMS(:,:) = cZero
  do N=2,d-1,2
    do M=0,N
      if (M == 0) then
        CNMC(N,M) = Half*(C(N,M)+C(N,-M))
      else
        m_fact = (-One)**M
        CNMC(N,M) = C(N,M)+m_fact*C(N,-M)
        CNMS(N,M) = -Onei*(C(N,M)-m_fact*C(N,-M))
      end if
    end do
  end do

  write(u6,'(A)') 'The ZFS Hamiltonian:'
  write(u6,'(A)') '   ZFS = SUM_{n,m}: [ E(n,m) * O(n,m) +  F(n,m) * W(n,m) ]'
  write(u6,'(A)') 'where:'
  write(u6,'(A)') '   O(n,m) =  0.5 * ( (-1)**m * Y(n,+m) + Y(n,-m) );'
  write(u6,'(A)') '   W(n,m) = -0.5 * ( (-1)**m * Y(n,+m) - Y(n,-m) ) * I;    (I = imaginary unit)'
  write(u6,'(A)') '   n - the rank of the ITO, = 2, 4, 6, ... 2*spin;'
  write(u6,'(A)') '   m - the component of the ITO, = 0, 1, ... n;'
  write(u6,'(A)') 'The quantization axis is the main magnetic axis of this multiplet (Zm).'
  write(u6,'(2A)') repeat('-',59),'|'
  write(u6,'(A)') '  n  |  m  |         E(n,m)        |         F(n,m)        |'
  do N=2,d-1,2
    write(u6,'(A)') '-----|-----|-----------------------|-----------------------|'
    do M=0,N
      write(u6,'(2(1x,I2,2x,A),2(ES22.14,1x,A))') N,'|',M,'|',real(CNMC(N,M)),'|',real(CNMS(N,M)),'|'
    end do
  end do
  write(u6,'(2A)') repeat('-',59),'|'

  ! decomposition of the ZFS matrix in Extended Stevens Operators
  write(u6,'(A)') repeat('*',80)
  write(u6,'(A)') 'The ZFS Hamiltonian:'
  write(u6,'(A)') '   ZFS = SUM_{k,q} * [ B(k,q) * O(k,q) ];'
  write(u6,'(A)') 'where:'
  write(u6,'(A)') '   O(k,q) =  Extended Stevens Operators (ESO) as defined in:'
  write(u6,'(10x,A)') '1. Rudowicz, C.; J.Phys.C: Solid State Phys.,18(1985) 1415-1430.'
  write(u6,'(10x,A)') '2. Implemented in the "EasySpin" function in MATLAB, www.easyspin.org.'
  write(u6,'(A    )') '   k - the rank of the ITO, = 2, 4, 6, 8, 10, 12.'
  write(u6,'(A)') '   q - the component of the ITO, = -k, -k+1, ... 0, 1, ... k;'

  if ((d-1) > 12) then
    write(u6,'(A)') 'k = 12 may not be the highest rank of the ITO for this case, but it '
    write(u6,'(A)') 'is the maximal k implemented in the "EasySpin" function in MATLAB.'
  end if

  write(u6,'(A)') 'Knm are proportionality coefficients between the ESO and operators defined in '
  write(u6,'(A)') 'J. Chem. Phys., 137, 064112 (2012).'
  write(u6,'(2A)') repeat('-',48),'|'
  write(u6,'(A)') '  k |  q  |    (Knm)^2  |         B(k,q)        |'
  if (d-1 > 12) then
    Nmax = 12
  else
    Nmax = d-1
  end if
  do N=2,Nmax,2
    write(u6,'(A)') '----|-----|-------------|-----------------------|'
    do M=-N,N
      if (M < 0) then
        write(u6,'((1x,I2,1x,A),(1x,I3,1x,A),F11.2,2x,A,2(ES22.14,1x,A))') N,'|',M,'|',knm(n,abs(m))*knm(n,abs(m)),'|', &
                                                                           real(CNMS(N,abs(M)))*knm(n,abs(m)),'|'
      else
        write(u6,'((1x,I2,1x,A),(1x,I3,1x,A),F11.2,2x,A,2(ES22.14,1x,A))') N,'|',M,'|',knm(n,abs(m))*knm(n,abs(m)),'|', &
                                                                           real(CNMC(N,abs(M)))*Knm(n,abs(m)),'|'
      end if
    end do
  end do
  write(u6,'(2A)') repeat('-',48),'|'
  !-----------------------------
  ! for the interface related to CF gradient calculation:
  if (GRAD) then
    LuDgrad = IsFreeUnit(81)
    call molcas_open(LuDgrad,'DMAT')
    do N=2,Nmax,2
      write(u6,'(A)') '----|-----|-------------|-----------------------|'
      do M=-N,N
        if (M < 0) then
          write(LuDgrad,'(I4,I4,1x,2(ES25.15))') N,M,real(CNMS(n,abs(n)))*knm(n,abs(m))
        else
          write(LuDgrad,'(I4,I4,1x,2(ES25.15))') N,M,real(CNMC(n,abs(m)))*Knm(n,abs(m))
        end if
      end do
    end do
    close(LuDgrad)
  end if
  !-----------------------------

  ES(:) = CNMC(2,0:2)
  FS(:) = CNMS(2,0:2)

  call mma_deallocate(C)
  call mma_deallocate(CNMC)
  call mma_deallocate(CNMS)

  call DMATRIX(Es,Fs,maxes,2)
end if ! decomposition of ZFS in higher-order ITO operators

write(u6,*)
write(u6,'(A,I2,A,I3,A)') 'ANGULAR MOMENTS ALONG THE MAIN MAGNETIC AXES'
call moments(d,s_so2,dipso2,iprint)

!----------------------------------------------------------------------
call mma_deallocate(ELOC)
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
call mma_deallocate(AMS_TMP)
call mma_deallocate(AMSSPIN)
call mma_deallocate(DIPSO2)
call mma_deallocate(S_SO2)
call mma_deallocate(HCF2)
call mma_deallocate(SP_DIPO)
call mma_deallocate(SP_DIPW)

return

end subroutine G_HIGH_1
