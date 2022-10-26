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

subroutine JACOB(ARRAY,VECS,NDIM,LENVEC)
!subroutine JACOBI(ARRAY,VECS,NDIM,LENVEC)

implicit real*8(A-H,O-Z)
dimension ARRAY(*)
dimension VECS(LENVEC,*)
parameter(EPS=1.0D-16)
parameter(EPS2=1.0D-30)
logical IFTEST

IfTest = .false.
#ifdef _DEBUGPRINT_
IfTest = .true.
#endif

if (NDIM <= 1) return
! A shift is applied. Use a representative diagonal value:
NDTRI = (NDIM*(NDIM+1))/2
SHIFT = 0.5d0*(ARRAY(1)+ARRAY(NDTRI))
SHIFT = dble(nint(SHIFT))
II = 0
do I=1,NDIM
  II = II+I
  ARRAY(II) = ARRAY(II)-SHIFT
end do

! Sanity test:
call CHK4NAN(NDIM*(NDIM+1)/2,Array,Ierr)
if (IERR /= 0) then
  call ABEND()
end if

if (IFTEST) then
  write(6,*) ' JACOBI test prints:'
  write(6,*) ' NSWEEP = Nr of sweeps'
  write(6,*) ' NR     = Rotations this sweep'
  write(6,*) ' NROT   = Nr of 2-rotations (Accum)'
  write(6,*) ' VNSUM  = von Neumann''s sum'
  write(6,*) ' SBDMAX = Largest subdiag element in this sweep.'
  write(6,*) '   NSWEEP   NR     NROT       VNSUM               SBDMAX'

  SBDMAX = 0.0d0
  VNSUM = 0.0d0
  do I=2,NDIM
    II = (I*(I-1))/2
    do J=1,I-1
      ARIJ = ARRAY(II+J)
      VNSUM = VNSUM+ARIJ**2
      SBDMAX = max(SBDMAX,abs(ARIJ))
    end do
  end do
  write(6,'(A,2G20.12)') ' Initial values:        ',VNSUM,SBDMAX
end if

! NSWEEP counts number of sweeps over subdiagonal elements.
NSWEEP = 0
! NROT counts total number of 2x2 rotations
NROT = 0

! Head of loop over sweeps
10 continue
NSWEEP = NSWEEP+1
! NR: Nr of 2x2 rotations in this sweep
NR = 0

SUBDAC = 0d0
NSUBD = 0
do I=2,NDIM
  II = (I*(I-1))/2
  do J=1,I-1
    !TEST do IMJ=1,NDIM-1
    !TEST   do I=IMJ+1,NDIM
    !TEST     II = I*(I-1)/2
    !TEST     J = I-IMJ
    JJ = (J*(J-1))/2
    ARIJ = ARRAY(II+J)
    AAIJ = abs(ARIJ)
    ARII = ARRAY(II+I)
    ARJJ = ARRAY(JJ+J)
    DIFF = ARII-ARJJ
    SGN = 1.0d0
    if (DIFF < 0.0d0) then
      DIFF = -DIFF
      SGN = -SGN
    end if
    SUBDAC = SUBDAC+AAIJ
    NSUBD = NSUBD+1
    ! Decide if we should rotate: SUBDAC is accumulated sum of abs of subdiag
    ! values. Therefore, we are certain that they are distributed around the
    ! value SUBDAC/NSUBD. Skip rotations that are too small compared to average:
    if (dble(NSUBD)*AAIJ <= 0.5*SUBDAC) goto 100
    ! Or, if the resulting rotation would be insignificant compared to DIFF:
    if (AAIJ <= EPS*DIFF) goto 100
    ! Or, if the resulting rotation would be insignificant in absolute size:
    if (AAIJ <= EPS2) goto 100
    ! Determine size of 2x2 rotation:
    NR = NR+1
    DUM = DIFF+sqrt(DIFF**2+4.0d0*AAIJ**2)
    TN = 2.0d0*SGN*ARIJ/DUM
    CS = 1.0d0/sqrt(1.0d0+TN**2)
    SN = CS*TN
    ! TN,CS,SN=TAN,COS AND SIN OF ROTATION ANGLE.
    ! The following partially unrolled loops are written by Markus Fuelscher
    kmin = 1
    kmax = j-1
    kleft = mod(kmax-kmin+1,4)
    if (kleft == 1) then
      Aii1 = ARRAY(ii+1)
      Ajj1 = ARRAY(jj+1)
      ARRAY(ii+1) = SN*Ajj1+CS*Aii1
      ARRAY(jj+1) = CS*Ajj1-SN*Aii1
    else if (kleft == 2) then
      Ajj1 = ARRAY(jj+1)
      Ajj2 = ARRAY(jj+2)
      Aii1 = ARRAY(ii+1)
      Aii2 = ARRAY(ii+2)
      ARRAY(ii+1) = SN*Ajj1+CS*Aii1
      ARRAY(ii+2) = SN*Ajj2+CS*Aii2
      ARRAY(jj+1) = CS*Ajj1-SN*Aii1
      ARRAY(jj+2) = CS*Ajj2-SN*Aii2
    else if (kleft == 3) then
      Ajj1 = ARRAY(jj+1)
      Ajj2 = ARRAY(jj+2)
      Ajj3 = ARRAY(jj+3)
      Aii1 = ARRAY(ii+1)
      Aii2 = ARRAY(ii+2)
      Aii3 = ARRAY(ii+3)
      ARRAY(ii+1) = SN*Ajj1+CS*Aii1
      ARRAY(ii+2) = SN*Ajj2+CS*Aii2
      ARRAY(ii+3) = SN*Ajj3+CS*Aii3
      ARRAY(jj+1) = CS*Ajj1-SN*Aii1
      ARRAY(jj+2) = CS*Ajj2-SN*Aii2
      ARRAY(jj+3) = CS*Ajj3-SN*Aii3
    end if
    kmin = kmin+kleft
    do k=kmin,kmax,4
      Ajj0 = ARRAY(jj+k+0)
      Ajj1 = ARRAY(jj+k+1)
      Ajj2 = ARRAY(jj+k+2)
      Ajj3 = ARRAY(jj+k+3)
      Aii0 = ARRAY(ii+k+0)
      Aii1 = ARRAY(ii+k+1)
      Aii2 = ARRAY(ii+k+2)
      Aii3 = ARRAY(ii+k+3)
      ARRAY(ii+k+0) = SN*Ajj0+CS*Aii0
      ARRAY(ii+k+1) = SN*Ajj1+CS*Aii1
      ARRAY(ii+k+2) = SN*Ajj2+CS*Aii2
      ARRAY(ii+k+3) = SN*Ajj3+CS*Aii3
      ARRAY(jj+k+0) = CS*Ajj0-SN*Aii0
      ARRAY(jj+k+1) = CS*Ajj1-SN*Aii1
      ARRAY(jj+k+2) = CS*Ajj2-SN*Aii2
      ARRAY(jj+k+3) = CS*Ajj3-SN*Aii3
    end do
    kmin = j+1
    kmax = i-1
    kleft = mod(kmax-kmin+1,4)
    kk = jj+j
    if (kleft == 1) then
      k0 = kk
      Ak0j = ARRAY(k0+j)
      Aii0 = ARRAY(ii+kmin+0)
      ARRAY(ii+kmin+0) = SN*Ak0j+CS*Aii0
      ARRAY(kk+j) = CS*Ak0j-SN*Aii0
      kk = k0+kmin
    else if (kleft == 2) then
      k0 = kk
      k1 = k0+kmin
      Ak0j = ARRAY(k0+j)
      Ak1j = ARRAY(k1+j)
      Aii0 = ARRAY(ii+kmin+0)
      Aii1 = ARRAY(ii+kmin+1)
      ARRAY(k0+j) = CS*Ak0j-SN*Aii0
      ARRAY(k1+j) = CS*Ak1j-SN*Aii1
      ARRAY(ii+kmin+0) = SN*Ak0j+CS*Aii0
      ARRAY(ii+kmin+1) = SN*Ak1j+CS*Aii1
      kk = k1+kmin+1
    else if (kleft == 3) then
      k0 = kk
      k1 = k0+kmin
      k2 = k1+kmin+1
      Ak0j = ARRAY(k0+j)
      Ak1j = ARRAY(k1+j)
      Ak2j = ARRAY(k2+j)
      Aii0 = ARRAY(ii+kmin+0)
      Aii1 = ARRAY(ii+kmin+1)
      Aii2 = ARRAY(ii+kmin+2)
      ARRAY(k0+j) = CS*Ak0j-SN*Aii0
      ARRAY(k1+j) = CS*Ak1j-SN*Aii1
      ARRAY(k2+j) = CS*Ak2j-SN*Aii2
      ARRAY(ii+kmin+0) = SN*Ak0j+CS*Aii0
      ARRAY(ii+kmin+1) = SN*Ak1j+CS*Aii1
      ARRAY(ii+kmin+2) = SN*Ak2j+CS*Aii2
      kk = k2+kmin+2
    end if
    kmin = kmin+kleft
    do k=kmin,kmax,4
      k0 = kk
      k1 = k0+k
      k2 = k1+k+1
      k3 = k2+k+2
      Ak0j = ARRAY(k0+j)
      Ak1j = ARRAY(k1+j)
      Ak2j = ARRAY(k2+j)
      Ak3j = ARRAY(k3+j)
      Aii0 = ARRAY(ii+k+0)
      Aii1 = ARRAY(ii+k+1)
      Aii2 = ARRAY(ii+k+2)
      Aii3 = ARRAY(ii+k+3)
      ARRAY(k0+j) = CS*Ak0j-SN*Aii0
      ARRAY(k1+j) = CS*Ak1j-SN*Aii1
      ARRAY(k2+j) = CS*Ak2j-SN*Aii2
      ARRAY(k3+j) = CS*Ak3j-SN*Aii3
      ARRAY(ii+k+0) = SN*Ak0j+CS*Aii0
      ARRAY(ii+k+1) = SN*Ak1j+CS*Aii1
      ARRAY(ii+k+2) = SN*Ak2j+CS*Aii2
      ARRAY(ii+k+3) = SN*Ak3j+CS*Aii3
      kk = k3+k+3
    end do
    kmin = i+1
    kmax = NDIM
    kleft = mod(kmax-kmin+1,4)
    kk = ii+i
    if (kleft == 1) then
      k0 = kk
      Ak0j = ARRAY(k0+j)
      Ak0i = ARRAY(k0+i)
      ARRAY(k0+j) = CS*Ak0j-SN*Ak0i
      ARRAY(k0+i) = SN*Ak0j+CS*Ak0i
      kk = k0+kmin
    else if (kleft == 2) then
      k0 = kk
      k1 = k0+kmin
      Ak0j = ARRAY(k0+j)
      Ak0i = ARRAY(k0+i)
      Ak1j = ARRAY(k1+j)
      Ak1i = ARRAY(k1+i)
      ARRAY(k0+j) = CS*Ak0j-SN*Ak0i
      ARRAY(k0+i) = SN*Ak0j+CS*Ak0i
      ARRAY(k1+j) = CS*Ak1j-SN*Ak1i
      ARRAY(k1+i) = SN*Ak1j+CS*Ak1i
      kk = k1+kmin+1
    else if (kleft == 3) then
      k0 = kk
      k1 = k0+kmin
      k2 = k1+kmin+1
      Ak0j = ARRAY(k0+j)
      Ak0i = ARRAY(k0+i)
      Ak1j = ARRAY(k1+j)
      Ak1i = ARRAY(k1+i)
      Ak2j = ARRAY(k2+j)
      Ak2i = ARRAY(k2+i)
      ARRAY(k0+j) = CS*Ak0j-SN*Ak0i
      ARRAY(k0+i) = SN*Ak0j+CS*Ak0i
      ARRAY(k1+j) = CS*Ak1j-SN*Ak1i
      ARRAY(k1+i) = SN*Ak1j+CS*Ak1i
      ARRAY(k2+j) = CS*Ak2j-SN*Ak2i
      ARRAY(k2+i) = SN*Ak2j+CS*Ak2i
      kk = k2+kmin+2
    end if
    kmin = kmin+kleft
    do k=kmin,kmax,4
      k0 = kk
      k1 = k0+k
      k2 = k1+k+1
      k3 = k2+k+2
      Ak0j = ARRAY(k0+j)
      Ak0i = ARRAY(k0+i)
      Ak1j = ARRAY(k1+j)
      Ak1i = ARRAY(k1+i)
      Ak2j = ARRAY(k2+j)
      Ak2i = ARRAY(k2+i)
      Ak3j = ARRAY(k3+j)
      Ak3i = ARRAY(k3+i)
      ARRAY(k0+j) = CS*Ak0j-SN*Ak0i
      ARRAY(k0+i) = SN*Ak0j+CS*Ak0i
      ARRAY(k1+j) = CS*Ak1j-SN*Ak1i
      ARRAY(k1+i) = SN*Ak1j+CS*Ak1i
      ARRAY(k2+j) = CS*Ak2j-SN*Ak2i
      ARRAY(k2+i) = SN*Ak2j+CS*Ak2i
      ARRAY(k3+j) = CS*Ak3j-SN*Ak3i
      ARRAY(k3+i) = SN*Ak3j+CS*Ak3i
      kk = k3+k+3
    end do
    ! Update the diagonal elements of A
    Temp = 2.0d0*CS*SN*ARIJ
    CS2 = CS**2
    SN2 = SN**2
    ARRAY(jj+j) = SN2*ARII+CS2*ARJJ-Temp
    ARRAY(ii+j) = 0.0d0
    ARRAY(ii+i) = CS2*ARII+SN2*ARJJ+Temp
    ! Update rows/columns of the eigenvectors VECS
    kmin = 1
    kmax = LENVEC
    kleft = mod(kmax-kmin+1,4)
    if (kleft == 1) then
      C1j = VECS(1,j)
      C1i = VECS(1,i)
      VECS(1,i) = SN*C1j+CS*C1i
      VECS(1,j) = CS*C1j-SN*C1i
    else if (kleft == 2) then
      C1j = VECS(1,j)
      C2j = VECS(2,j)
      C1i = VECS(1,i)
      C2i = VECS(2,i)
      VECS(1,i) = SN*C1j+CS*C1i
      VECS(2,i) = SN*C2j+CS*C2i
      VECS(1,j) = CS*C1j-SN*C1i
      VECS(2,j) = CS*C2j-SN*C2i
    else if (kleft == 3) then
      C1j = VECS(1,j)
      C2j = VECS(2,j)
      C3j = VECS(3,j)
      C1i = VECS(1,i)
      C2i = VECS(2,i)
      C3i = VECS(3,i)
      VECS(1,i) = SN*C1j+CS*C1i
      VECS(2,i) = SN*C2j+CS*C2i
      VECS(3,i) = SN*C3j+CS*C3i
      VECS(1,j) = CS*C1j-SN*C1i
      VECS(2,j) = CS*C2j-SN*C2i
      VECS(3,j) = CS*C3j-SN*C3i
    end if
    kmin = kmin+kleft
    do k=kmin,kmax,4
      C0j = VECS(k+0,j)
      C1j = VECS(k+1,j)
      C2j = VECS(k+2,j)
      C3j = VECS(k+3,j)
      C0i = VECS(k+0,i)
      C1i = VECS(k+1,i)
      C2i = VECS(k+2,i)
      C3i = VECS(k+3,i)
      VECS(k+0,i) = SN*C0j+CS*C0i
      VECS(k+1,i) = SN*C1j+CS*C1i
      VECS(k+2,i) = SN*C2j+CS*C2i
      VECS(k+3,i) = SN*C3j+CS*C3i
      VECS(k+0,j) = CS*C0j-SN*C0i
      VECS(k+1,j) = CS*C1j-SN*C1i
      VECS(k+2,j) = CS*C2j-SN*C2i
      VECS(k+3,j) = CS*C3j-SN*C3i
    end do
    ! End of code by Markus Fuelscher

100 continue
  end do
end do
NROT = NROT+NR

if (IFTEST) then
  SBDMAX = 0.0d0
  VNSUM = 0.0d0
  do I=2,NDIM
    II = (I*(I-1))/2
    do J=1,I-1
      ARIJ = ARRAY(II+J)
      VNSUM = VNSUM+ARIJ**2
      SBDMAX = max(SBDMAX,abs(ARIJ))
    end do
  end do
  write(6,'(3I8,2G20.12)') NSWEEP,NR,NROT,VNSUM,SBDMAX
end if

! CHECK IF CONVERGED:
if (NR > 0) goto 10

! The shifted diagonal values are restored:
II = 0
do I=1,NDIM
  II = II+I
  ARRAY(II) = ARRAY(II)+SHIFT
end do

! STOP
return

end subroutine JACOB
