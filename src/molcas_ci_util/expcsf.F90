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

subroutine EXPCSF(ICS,NLEV,IMS,LEX,coef,LuVecDet)
! Expand IMS component of ICS(NLEV) in determinants using the
! procedure from Shavitt in "The Unitary Group", "Lecture Notes in
! Chemistry" Vol. 22, pp. 55.

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NLEV, ICS(NLEV), IMS, LuVecDet
integer(kind=iwp), intent(out) :: LEX(NLEV)
real(kind=wp), intent(in) :: coef
integer(kind=iwp) :: I, IA, IALPHA, IB, IBETA, ICOEF(2), ILEX, ISOMO, J, K, NALPHA, NSOMO
real(kind=wp) :: qdet
character(len=256) :: LINE
character(len=6) :: STRING
logical(kind=iwp) :: LAST, lod39, LPHASE

! find number of singly occupied orbitals
NSOMO = 0
do I=1,NLEV
  if ((ICS(I) == 1) .or. (ICS(I) == 2)) NSOMO = NSOMO+1
end do
NALPHA = (NSOMO+IMS)/2
call INIT_LEX(NALPHA,LEX)
! Loop over possible determinants
LAST = .false.
do while (.not. LAST)
  ICOEF(1) = 1
  ICOEF(2) = 1
  LPHASE = .true.
  IALPHA = 0
  IBETA = 0
  ISOMO = 0
  ILEX = 1
  IA = 0
  IB = 0
  LINE = ' '
  K = 26
  write(LINE(K:K),'(A)') '|'
  K = K+NLEV+1
  write(LINE(K:K),'(A)') '|'
  K = 26
  do I=1,NLEV
    K = K+1
    if (ICS(I) == 0) then
      LINE(K:K) = '0'
    else if (ICS(I) == 1) then
      IB = IB+1
      ISOMO = ISOMO+1
      if (ISOMO == LEX(ILEX)) then
        ILEX = ILEX+1
        ICOEF(1) = ICOEF(1)*(IA+IB-IBETA)
        IALPHA = IALPHA+1
        LINE(K:K) = 'a'
      else
        ICOEF(1) = ICOEF(1)*(IA+IB-IALPHA)
        IBETA = IBETA+1
        LINE(K:K) = 'b'
      end if
      ICOEF(2) = ICOEF(2)*IB
    else if (ICS(I) == 2) then
      IA = IA+1
      IB = IB-1
      ISOMO = ISOMO+1
      if (ISOMO == LEX(ILEX)) then
        ILEX = ILEX+1
        ICOEF(1) = ICOEF(1)*(IBETA-IA+1)
        if (mod(IB,2) == 0) LPHASE = .not. LPHASE
        IALPHA = IALPHA+1
        LINE(K:K) = 'a'
      else
        ICOEF(1) = ICOEF(1)*(IALPHA-IA+1)
        if (mod(IB,2) /= 0) LPHASE = .not. LPHASE
        IBETA = IBETA+1
        LINE(K:K) = 'b'
      end if
      ICOEF(2) = ICOEF(2)*(IB+2)
    else
      IA = IA+1
      IALPHA = IALPHA+1
      IBETA = IBETA+1
      if (mod(IB,2) /= 0) LPHASE = .not. LPHASE
      LINE(K:K) = '2'
    end if
    if (ICOEF(1) == 0) exit
  end do
  if (ICOEF(1) /= 0) then
    ! If non-zero coefficient, print the thing
    call SIMPLIFY(ICOEF)
    if (LPHASE) then
      write(LINE(1:9),'(2X,A7)') '+ sqrt('
    else
      write(LINE(1:9),'(2X,A7)') '- sqrt('
    end if
    write(LINE(10:16),'(I6,A1)') ICOEF(1),'/'
    write(STRING(1:6),'(I6)') ICOEF(2)
    J = 17
    do I=1,6
      if (STRING(I:I) /= ' ') then
        LINE(J:J) = STRING(I:I)
        J = J+1
      end if
    end do
    write(LINE(23:23),'(A1)') ')'
    write(u6,*) LINE(1:K+1)
    ! Write to GronOR vecdet files
    if (LuVecDet > 0) then
      inquire(unit=LuVecDet,opened=lod39)
      if (lod39) then
        qdet = coef*sqrt(real(icoef(1),kind=wp)/real(icoef(2),kind=wp))
        if (LPHASE) then
          write(LuVecDet,'(es15.8,6x,a)') qdet,LINE(27:K)
        else
          write(LuVecDet,'(es15.8,6x,a)') -qdet,LINE(27:K)
        end if
      end if
    end if
    ! End of print-out
  end if
  ! Get the next determinant
  call LEX_ITER(NSOMO,NALPHA,LEX,LAST)
end do

end subroutine EXPCSF
