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
! Copyright (C) 2026, Lila Zapp                                        *
!***********************************************************************

subroutine expkap_localisation(kappa,nOrb2Loc,Umat)
! analogous to exp_series in scf, only for the specified orbital subspace
!
! purpose: returns exp(kappa) as Umat
!
! kappa_cnt and xkappa_cnt are auxiliary matrices of dim (norb2loc,norb2loc) required for the Taylor expansion.
! they need to be allocated outside of this subroutine

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nOrb2Loc
real(kind=wp), intent(inout) :: kappa(nOrb2Loc,nOrb2Loc), Umat(nOrb2Loc,nOrb2Loc)
integer(kind=iwp) :: i, j

!Call RecPrt('Kappa',' ',Kappa,nOrb2Loc,nOrb2Loc)
do i=1,nOrb2Loc
  do j=1,i-1
    if (abs(Kappa(i,j)+Kappa(j,i)) > 1.0e-16_wp) then
      write(u6,*) ' Mismatch of off-diagonal value in Kappa matrix'
      call Abend()
    end if
  end do
  if (Kappa(i,i) /= Zero) then
    write(u6,*) ' Error in diagonal value in Kappa matrix'
    call Abend()
  end if
end do

#define _OLD_
#ifdef _OLD_
call expkap_localisation_Local()
#else
block
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: One
  integer(kind=iwp) ideg, m, ldh, lwsp, iexph, ns, iflag, i, j
  real(kind=wp) t
  real(kind=wp), allocatable :: wsp(:)
  integer(kind=iwp), allocatable :: ipiv(:)
  ideg = 13
  m = nOrb2Loc
  t = One
  ldh = m
  lwsp = 4*m**2+ideg+1
  call mma_allocate(wsp,lwsp,Label='wsp')
  call mma_allocate(ipiv,m,Label='ipiv')
  call dgpadm(ideg,m,t,kappa,ldh,wsp,lwsp,ipiv,iexph,ns,iflag)

  if (iflag == 0) then
    do j=1,m
      do i=1,m
        Umat(i,j) = wsp(iexph)
        iexph = iexph+1
      end do
    end do
  else
    write(u6,*) 'DGPADM: iflag=',iflag
    call Abend()
  end if
  call mma_deallocate(ipiv)
  call mma_deallocate(wsp)
end block
#endif

!#define _TEST_
#ifdef _TEST_
block
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: One
  real(kind=wp), allocatable :: UnitMatrix(:,:)
  call mma_allocate(UnitMatrix,nOrb2Loc,nOrb2Loc,Label='UM')
  UnitMatrix(:,:) = Zero

  call dgemm_('N','T',nOrb2Loc,nOrb2Loc,nOrb2Loc, &
              One,UMat,nOrb2Loc, &
                  UMat,nOrb2Loc, &
              Zero,UnitMatrix,nOrb2Loc)

  do i=1,nOrb2Loc
    do j=1,i-1
      if (abs(UnitMatrix(i,j)) > 1.0e-12_wp) then
        write(u6,*) 'Too large off diagonal value'
        call Abend()
      end if
      if (abs(UnitMatrix(j,i)) > 1.0e-12_wp) then
        write(u6,*) 'Too large off diagonal value'
        call Abend()
      end if
    end do
    if (abs(One-UnitMatrix(i,i)) > 1.0e-12_wp) then
      write(u6,*) 'Diagonal not accurate!'
      write(u6,*) 'UM(i,i)=',UnitMatrix(i,i)
      call Abend()
    end if
  end do
  call mma_deallocate(UnitMatrix)
end block
#endif

contains

subroutine expkap_localisation_Local()
  ! analogous to exp_series in scf, only for the specified orbital subspace
  !
  ! purpose: returns exp(kappa) as Umat
  !
  ! kappa_cnt and xkappa_cnt are auxiliary matrices of dim (norb2loc,norb2loc) required for the Taylor expansion.
  ! they need to be allocated outside of this subroutine

  use Localisation_globals, only: kappa_cnt, xkappa_cnt
  use Constants, only: One

  integer(kind=iwp) :: cnt, maxel(2)
  real(kind=wp) :: factor, ithrsh
  real(kind=wp), parameter :: thrsh_taylor = 1.0e-16_wp
  logical(kind=iwp), parameter :: debug_exp = .false.
  real(kind=wp), external :: DDot_

  kappa_cnt(:,:) = kappa(:,:) ! kappa^cnt = kappa since cnt=1
  xkappa_cnt(:,:) = kappa(:,:)

  call unitmat(Umat,nOrb2Loc)

  cnt = 1
  factor = One
  ithrsh = One
  maxel(:) = 0

  Umat(:,:) = Umat(:,:)-kappa(:,:)

  if (debug_exp) then
    write(u6,*) 'Taylor expansion: n=1'
    call RecPrt('Umat = I + kappa^1',' ',Umat(:,:),nOrb2Loc,nOrb2Loc)
    call RecPrt('kappa',' ',kappa(:,:),nOrb2Loc,nOrb2Loc)
    write(u6,*) 'Taylor expansion: more terms'
  end if

  do while (ithrsh > thrsh_taylor)

    ! the number of the term = the exponent for kappa in that term
    cnt = cnt+1

    ! the faculty value that the matrix will be divided by
    factor = factor*real(cnt,kind=wp)

    ! calculate the cnt'th exponent of the kappa matrix
    ! initial kappa matrix (just kappa^1)
    ! C <= alpha*A*B + beta*C
    ! kappa_cnt <= 1*kappa_cnt*kappa + 0*kappa_cnt
    call dgemm_('N','N',nOrb2Loc,nOrb2Loc,nOrb2Loc, &
                One/real(cnt,kind=wp),xkappa_cnt,nOrb2Loc, &
                                      kappa,nOrb2Loc, &
                Zero,kappa_cnt,norb2Loc)

    xkappa_cnt(:,:) = kappa_cnt(:,:)

    ! differentiation of odd and even cases, because this expands exp(-kappa)
    ! all terms starting at n=2
    if (mod(cnt,2) == 0) then
      Umat(:,:) = Umat(:,:)+kappa_cnt(:,:)
      if (debug_exp) write(u6,'(A,ES10.1,A,I2,A,ES12.4)') 'term: +',One/factor,' * kappa^',cnt,', current ithrsh = ',ithrsh
    else
      Umat(:,:) = Umat(:,:)-kappa_cnt(:,:)
      if (debug_exp) write(u6,'(A,ES10.1,A,I2,A,ES12.4)') 'term: -',One/factor,' * kappa^',cnt,', current ithrsh = ',ithrsh
    end if

    ithrsh = 1.0e6_wp*sqrt(DDot_(nOrb2Loc**2,Kappa_Cnt,1,Kappa_Cnt,1))/(Norb2loc**2)
    !ithrsh = maxval(abs(Kappa_Cnt(:,:))/(abs(Umat(:,:))+thrsh_taylor))

    ! sanity check for divergence
    if (ithrsh/720.0e6_wp > One) then
      write(u6,*) 'Bug: the Taylor expansion of exp(kappa) diverges - numerical error'
      write(u6,*) 'Stopping Taylor expansion at ',cnt,'-th term. Rescale kappa before feeding it this subroutine. ithrsh',ithrsh
      call Abend()
    end if

    if (debug_exp) then
      write(u6,'(A,ES10.1,A,I2,A,ES12.4)') 'term:  ',One/factor,' * kappa^',cnt,', new ithrsh     = ',ithrsh
      call RecPrt('kappa^cnt',' ',kappa_cnt,nOrb2Loc,nOrb2Loc)
      call RecPrt('Umat',' ',Umat,nOrb2Loc,nOrb2Loc)
    end if

  end do

  maxel(:) = maxloc(abs(Umat(:,:)))
  if (Umat(maxel(1),maxel(2)) > One+1.0e-8_wp) then
    write(u6,*) 'element of U bigger cannot be larger than 1, as U is supposed to be orthogonal/unitary:',Umat(maxel(1),maxel(2))
    call Abend()
  end if

# ifdef _DEBUGPRINT_
  kappa_cnt(:,:) = Zero

  call dgemm_('N','T',nOrb2Loc,nOrb2Loc,nOrb2Loc, &
              One,Umat,nOrb2Loc, &
                  Umat,nOrb2Loc, &
              Zero,kappa_cnt,norb2Loc)
  call RecPrt('kappa',' ',kappa(:,:),nOrb2Loc,nOrb2Loc)
  call RecPrt('Umat',' ',Umat(:,:),nOrb2Loc,nOrb2Loc)
  call RecPrt('U^(T)U = I',' ',kappa_cnt(:,:),nOrb2Loc,nOrb2Loc)

# endif

end subroutine expkap_localisation_local

end subroutine expkap_localisation
