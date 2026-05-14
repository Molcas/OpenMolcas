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
! Copyright (C) 1990,2020, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

!#define _DEBUGPRINT_
module Real_Spherical

use Index_Functions, only: C_Ind3, nTri_Elem1, nTri3_Elem1
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
private

integer(kind=iwp) :: lmax_internal = -1
logical(kind=iwp) :: Condon_Shortley_phase_factor = .false.
integer(kind=iwp), allocatable :: ipSph(:), iSphCr(:)
real(kind=wp), allocatable :: RSph(:)
character(len=8), allocatable :: LblCBs(:), LblSBs(:)

public :: ipSph, iSphCr, LblCBs, LblSBs, lmax_internal, RSph, Sphere, Sphere_Free

!                                                                      *
!***********************************************************************
!                                                                      *
contains
!                                                                      *
!***********************************************************************
!                                                                      *
subroutine Sphere_Free()

  call mma_deallocate(RSph,safe='*')
  call mma_deallocate(ipSph,safe='*')
  call mma_deallocate(iSphCr,safe='*')
  call mma_deallocate(LblCBs,safe='*')
  call mma_deallocate(LblSBs,safe='*')
  lmax_internal = -1

end subroutine Sphere_Free
!                                                                      *
!***********************************************************************
!                                                                      *
subroutine Sphere(lMax)
!***********************************************************************
!                                                                      *
! Object: create the transformation matrices from cartesian gaussians  *
!         to spherical gaussians. By having these matricies being      *
!         defined dynamical we ensure that any extension of the        *
!         program to higher angular momentum is simply done by chang-  *
!         ing iTabMx to the appropiate value.                          *
!         In addition, this will also allow us to have any order of    *
!         vectors in the matrix, i.e. we can have our own format or    *
!         any other odd order (MOLECULE).                              *
!                                                                      *
! Called from: Input                                                   *
!                                                                      *
! Calling    : Real_Sphere                                             *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March 1990                                               *
!***********************************************************************
!               Credits.                                               *
!               2020, R. Lindh; P. R. Taylor; L. Birnoschi; A. Dzubak; *
!                     M. Navarrete; C. Gonzalez-Espinoza; G. Raggi;    *
!                     N. F. Chilton at OpenMolcas2020                  *
!***********************************************************************

  use define_af, only: iTabMx

  integer(kind=iwp), intent(in) :: lMax
  integer(kind=iwp) :: iAng, iElem, ii, iii, jjj, l, m, MxFnc, n, nElem, nSphCr, nSphr
# ifdef _DEBUGPRINT_
  integer(kind=iwp) :: i, iLbl, j
# endif

  ! check if required ang mom is greater than hard-coded limit
  if (lMax > iTabMx) then
    call WarningMessage(2,' Sphere: Increase iTabMx!')
    call Abend()
  end if

  if (lmax < 0) then
    write(u6,*) 'Sphere: lmax<0'
    call Abend()
  end if
  if (lmax > lmax_internal) then
    call Sphere_Free()
    lmax_internal = lMax
  else
    return
  end if
  call Get_lScalar('CSPF',Condon_Shortley_phase_factor)

  nSphCr = nTri3_Elem1(lmax)
  call mma_allocate(iSphCr,nSphCr,Label='iSphCr')
  iSphCr(:) = 0

  !write(u6,*) 'C&S',Condon_Shortley_phase_factor

  ! Make the labels
  ! Gives info on basis function angular momenta
  ! n, l, ml or assigns it as a diffuse/polarising function with '*'

  MxFnc = nSphCr
  call mma_allocate(LblCBs,MxFnc,Label='LblCBs')
  call mma_allocate(LblSBs,MxFnc,Label='LblSBs')

  call Make_Labels(LblCbs,LblSbs,MxFnc,lMax)

  ! Allocate memory for transformation matrices
  ! Here, ipSph are the pointers to memory locations in RSph, for the
  ! transformation matrices of given ang mom
  nSphr = 0
  do iAng=0,lMax
    nSphr = nSphr+nTri_Elem1(iAng)**2
  end do
  call mma_allocate(RSph,nSphr,label='RSph')
  call mma_allocate(ipSph,[0,lMax],label='ipSph')
  ipSph(0) = 1
  do iAng=0,lMax-1
    ipSph(iAng+1) = ipSph(iAng)+nTri_Elem1(iAng)**2
  end do

  ! Here the transformation matrices from cartesian to spherical are made
  call Real_Sphere(ipSph,lMax,RSph,nSphr)

  ! Set up the symmetry properties of the spherical gaussians
  iii = 0
  jjj = 0
  do n=0,lMax
    nElem = nTri_Elem1(n)
    ii = 0
    do m=n,0,-2
      do l=-m,m
        iii = iii+1
        do iElem=1,nElem
          if (RSph(iElem-1+ii+ipSph(n)) /= Zero) exit
        end do
        iSphCr(iii) = iElem+jjj
        ii = ii+nElem
      end do
    end do
    jjj = jjj+nElem
  end do

# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) ' Spherical Harmonic expansions '
  write(u6,*)
  iLbl = 1
  do n=0,lMax
    nElem = nTri_Elem1(n)
    ii = 0
    write(u6,*)
    write(u6,'(8X,31(2X,I1,I1,I1))') ((i,j,n-i-j,j=n-i,0,-1),i=n,0,-1)
    write(u6,*)
    do m=n,0,-2
      do l=-m,m
        write(u6,'(1X,A6,1X,31F5.2)') LblSbs(iLbl),(RSph(i+ii+ipSph(n)),i=0,nElem-1)
        ii = ii+nElem
        iLbl = iLbl+1
      end do
      write(u6,*)
    end do
    write(u6,*)
  end do
# endif

end subroutine Sphere

subroutine Real_Sphere(ipSph,lMax,RSph,nSphr)

  integer(kind=iwp), intent(in) :: lMax, ipSph(0:lMax), nSphr
  real(kind=wp), intent(out) :: RSph(nSphr)
  integer(kind=iwp) :: i, i00, i10, i2, i20, iCont, iOff, j, l, mElem, nElem

  i00 = ipSph(0)
  i10 = ipSph(0)
  do i=0,lMax
    i2 = ipSph(i)
    nElem = nTri_Elem1(i)
    i20 = i2+i*nElem
    ! First generate the coefficients for Y(i,0) -- always real
    call Recurse(RSph(i00),RSph(i10),RSph(i20),i)
    ! Use ladder operators to generate Y(i,m)
    call Ladder(RSph(i2),i)

    ! Now do the contaminant, by simply multiply with r**2

    j = i-2
    if (j >= 0) then
      iCont = i2+(2*i+1)*nElem
      iOff = ipSph(j)
      mElem = nTri_Elem1(j)
      do l=j,0,-2
        call Contaminant(RSph(iCont),i,RSph(iOff),j,l)
        iCont = iCont+(2*l+1)*nElem
        iOff = iOff+(2*l+1)*mElem
      end do
    end if

    i00 = i10
    i10 = i20
  end do

  ! Normalize

  do i=0,lMax
    call NrmSph(RSph(ipSph(i)),i)
  end do

end subroutine Real_Sphere

subroutine Recurse(P0,P1,P2,n2)
!***********************************************************************
!                                                                      *
!     The Legendre polynomial is identical to Y(l,0).                  *
!     Note that it is real and that there is no Condon-Shortley phase  *
!     factor to consider.                                              *
!                                                                      *
!***********************************************************************

  integer(kind=iwp), intent(in) :: n2
  real(kind=wp), intent(in) :: P0(nTri_Elem1(n2-2)), P1(nTri_Elem1(n2-1))
  real(kind=wp), intent(out) :: P2(nTri_Elem1(n2))
  integer(kind=iwp) :: ix, iy, iz, n0, n1
  real(kind=wp) :: Fact_1, Fact_2

  P2(:) = Zero

  ! Use recurrence relation for Legendre polynomials
  !
  !   (n+1) P_{n+1} = (2n+1) z P_n - n r^2 P_{n-1}

  if (n2 == 0) then

    P2(1) = One

  else

    ! P_{n+1} = (2n+1)/(n+1) z P_n

    Fact_1 = real(2*n2-1,kind=wp)/real(n2,kind=wp)
    n1 = n2-1
    do ix=n1,0,-1
      do iy=n1-ix,0,-1
        iz = n1-ix-iy
        P2(C_Ind3(ix,iy,iz+1)) = P2(C_Ind3(ix,iy,iz+1))+Fact_1*P1(C_Ind3(ix,iy,iz))
      end do
    end do

    ! P_{n+1} = - n/(n+1) (x^2+y^2+z^2) P_{n-1}

    Fact_2 = real(n2-1,kind=wp)/real(n2,kind=wp)
    n0 = n1-1
    do ix=n0,0,-1
      do iy=n0-ix,0,-1
        iz = n0-ix-iy
        P2(C_Ind3(ix+2,iy,iz)) = P2(C_Ind3(ix+2,iy,iz))-Fact_2*P0(C_Ind3(ix,iy,iz))
        P2(C_Ind3(ix,iy+2,iz)) = P2(C_Ind3(ix,iy+2,iz))-Fact_2*P0(C_Ind3(ix,iy,iz))
        P2(C_Ind3(ix,iy,iz+2)) = P2(C_Ind3(ix,iy,iz+2))-Fact_2*P0(C_Ind3(ix,iy,iz))
      end do
    end do

  end if

end subroutine Recurse

subroutine Ladder(P0,n)

  integer(kind=iwp), intent(in) :: n
  real(kind=wp), intent(inout) :: P0(nTri_Elem1(n),-n:n)
  integer(kind=iwp) :: ix, iy, iz, m, m_m, m_p
  real(kind=wp) :: Fact

  ! Generate Y(l,m) from Y(l,m-1), starting the process from Y(l,0)

  do m=0,n-1
    m_p = m+1
    m_m = -(m+1)
    P0(:,m_p) = Zero
    P0(:,m_m) = Zero
    Fact = One/(Two*sqrt(real(n*(n+1)-m*(m-1),kind=wp)))

    ! The spherical harmonic is a two component (real,imaginary)
    ! function.
    !
    !   Y(n, m) =(-1)**  m  x (S(+,m), S(-,m)) and
    !   Y(n,-m) =(-1)**(-m) x (S(+,m),-S(-,m))
    !
    ! with S(-,0)=0
    !
    ! The ladder operator is subdivided in a similar way
    !
    !   L(+)=(Lr,Li),  L(-)=(Lr,-Li)
    !
    ! Hence
    !
    !   L(+) Y(n,m)= C x Y(n,m+1)
    !
    ! or
    !
    !   C(S(+,m+1),S(-,m+1))=(Lr S(+,m)-Li S(-,m),Li S(+,m)+Lr S(-,m))

    do ix=n,0,-1
      do iy=n-ix,0,-1
        iz = n-ix-iy

        ! Generating the real part

        if (iz >= 1) P0(C_Ind3(ix+1,iy,iz-1),m_p) = P0(C_Ind3(ix+1,iy,iz-1),m_p)+Fact*real(iz,kind=wp)*P0(C_Ind3(ix,iy,iz),m)
        if (ix >= 1) P0(C_Ind3(ix-1,iy,iz+1),m_p) = P0(C_Ind3(ix-1,iy,iz+1),m_p)-Fact*real(ix,kind=wp)*P0(C_Ind3(ix,iy,iz),m)
        if (m /= 0) then
          if (iz >= 1) P0(C_Ind3(ix,iy+1,iz-1),m_p) = P0(C_Ind3(ix,iy+1,iz-1),m_p)-Fact*real(iz,kind=wp)*P0(C_Ind3(ix,iy,iz),-m)
          if (iy >= 1) P0(C_Ind3(ix,iy-1,iz+1),m_p) = P0(C_Ind3(ix,iy-1,iz+1),m_p)+Fact*real(iy,kind=wp)*P0(C_Ind3(ix,iy,iz),-m)
        end if

        ! Generating the imaginary part

        if (iz >= 1) P0(C_Ind3(ix,iy+1,iz-1),m_m) = P0(C_Ind3(ix,iy+1,iz-1),m_m)+Fact*real(iz,kind=wp)*P0(C_Ind3(ix,iy,iz),m)
        if (iy >= 1) P0(C_Ind3(ix,iy-1,iz+1),m_m) = P0(C_Ind3(ix,iy-1,iz+1),m_m)-Fact*real(iy,kind=wp)*P0(C_Ind3(ix,iy,iz),m)
        if (m /= 0) then
          if (iz >= 1) P0(C_Ind3(ix+1,iy,iz-1),m_m) = P0(C_Ind3(ix+1,iy,iz-1),m_m)+Fact*real(iz,kind=wp)*P0(C_Ind3(ix,iy,iz),-m)
          if (ix >= 1) P0(C_Ind3(ix-1,iy,iz+1),m_m) = P0(C_Ind3(ix-1,iy,iz+1),m_m)-Fact*real(ix,kind=wp)*P0(C_Ind3(ix,iy,iz),-m)
        end if

      end do
    end do

    ! Up to this point we have been operating on the Legendre and
    ! associated Legendre polynomials. Let us now put in the
    ! Condon-Shortley phase factor

    if (Condon_Shortley_phase_factor .and. (mod(m+1,2) /= 0)) then
      !write(u6,*) 'C&S phase factor included.'
      P0(:,m_p) = -P0(:,m_p)
      P0(:,m_m) = -P0(:,m_m)
    end if

  end do ! m

end subroutine Ladder

subroutine Contaminant(P0,i,Px,j,l)
! This subroutine generates the lower ang mom contaminants for the
! given ang mom

  integer(kind=iwp), intent(in) :: i, j, l
  real(kind=wp), intent(inout) :: P0(nTri_Elem1(i),-l:l)
  real(kind=wp), intent(in) :: Px(nTri_Elem1(j),-l:l)
  integer(kind=iwp) :: ix, iy, iz, m

  ! Px = (x^2+y^2+z^2) x P0

  do m=-l,l
    P0(:,m) = Zero
    do ix=j,0,-1
      do iy=j-ix,0,-1
        iz = j-ix-iy
        P0(C_Ind3(ix+2,iy,iz),m) = P0(C_Ind3(ix+2,iy,iz),m)+Px(C_Ind3(ix,iy,iz),m)
        P0(C_Ind3(ix,iy+2,iz),m) = P0(C_Ind3(ix,iy+2,iz),m)+Px(C_Ind3(ix,iy,iz),m)
        P0(C_Ind3(ix,iy,iz+2),m) = P0(C_Ind3(ix,iy,iz+2),m)+Px(C_Ind3(ix,iy,iz),m)
      end do
    end do
  end do

end subroutine Contaminant

subroutine NrmSph(P,n)

  integer(kind=iwp), intent(in) :: n
  real(kind=wp), intent(inout) :: P(nTri_Elem1(n),nTri_Elem1(n))
  integer(kind=iwp) :: ijx, ijy, ijz, ix, iy, iz, jx, jy, jz, k, m, nt
  real(kind=wp) :: DF, rMax, temp, tmp
  real(kind=wp), external :: DblFac

  nt = nTri_Elem1(n)
  do m=1,nt
    rMax = Zero
    do k=1,nt
      if (abs(P(k,m)) > rMax) rMax = abs(P(k,m))
    end do
    do k=1,nt
      if (abs(P(k,m)) < 1.0e-12_wp*rMax) P(k,m) = Zero
    end do
    tmp = Zero
    do ijx=2*n,0,-2
      do ijy=2*n-ijx,0,-2
        ijz = 2*n-ijx-ijy
        DF = DblFac(ijx-1)*DblFac(ijy-1)*DblFac(ijz-1)
        temp = Zero
        do ix=min(n,ijx),max(0,ijx-n),-1
          jx = ijx-ix
          do iy=min(n-ix,ijy),max(0,ijy-n+jx),-1
            jy = ijy-iy
            iz = n-ix-iy
            jz = n-jx-jy
            temp = temp+P(C_Ind3(ix,iy,iz),m)*P(C_Ind3(jx,jy,jz),m)
          end do
        end do
        tmp = tmp+DF*temp
      end do
    end do
    P(:,m) = P(:,m)/sqrt(tmp)
  end do

end subroutine NrmSph
!                                                                      *
!***********************************************************************
!                                                                      *
end module Real_Spherical
