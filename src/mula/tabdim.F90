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
! Copyright (C) 1995, Niclas Forsberg                                  *
!               1999, Anders Bernhardsson                              *
!***********************************************************************

subroutine TabDim(nDim,nOsc,nTabDim)
!  Purpose:
!    Fill a nDim*nOsc matrix with binomial coefficients.
!    This matrix is then used to calculate the dimension
!    of a table containing excitations for nOsc oscillators.
!    nDim is the maximum sum of the oscillator quanta.
!
!  Input:
!    nDim    : Integer variable - maximun excitation.
!    nOsc    : Integer variable - number of dimensions.
!
!  Result:
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nDim, nOsc
integer(kind=iwp), intent(out) :: nTabDim
integer(kind=iwp) :: iRow, jCol
integer(kind=iwp), allocatable :: binomCoef(:,:)

! Initialize.

! Calculate binomial coefficients.
if (nDim > 0) then
  call mma_allocate(binomCoef,[0,nDim],[1,nOsc],label='binomCoef')
  do jCol=1,nOsc
    binomCoef(0,jCol) = 1
  end do
  do iRow=0,nDim
    binomCoef(iRow,1) = 1
  end do
  do jCol=2,nOsc
    do iRow=1,nDim
      binomCoef(iRow,jCol) = binomCoef(iRow-1,jCol)+binomCoef(iRow,jCol-1)
    end do
  end do

  ! Sum all elements of the nOsc'th column.
  nTabDim = 0
  do iRow=0,nDim
    nTabDim = nTabDim+binomCoef(iRow,nOsc)
  end do
  call mma_deallocate(binomCoef)
else
  nTabDim = 1
end if

end subroutine TabDim
!####
function iDetNr(iocc,graph2,nOsc,m)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iDetNr
integer(kind=iwp), intent(in) :: nOsc, iocc(nOsc), m, graph2(0:m,0:m,nOsc)
integer(kind=iwp) :: i, iqnew, iqold, n
! iocc: occupation vector
! graph2: vertex graph table
! nOsc: number of nodes
! m: number of quantas
! number of determinants in with lower number of quantas.

! Calculate the index of occupation string iocc

iqnew = 0
iqold = 0
n = 0
do i=1,nOsc
  iqnew = iqnew+iocc(i)
  n = n+graph2(iqnew,iqold,i)
  iqold = iqnew
end do

iDetNr = n

end function iDetNr

! Muln  is a set of subroutines that calculates
! <i|H|j> where H is a operator described
! M_1,2..,n(a_1+a^t_1)*(a_2+a^t_2)...(a_n+a^t_n)
! you can find n=1,2,3,4 in this file.
!
! nmat : Occupation of slater det.
! F    : Output <i|A|j>
! iCre,iAnn : Gives the resulting slater determinant if a^t (a) is acting on SD
! mat  : Matrix describing the operator expanded in Normal modes
! m_Ord: Number of slater determinants.
!
! Anders Bernhardsson Friday the 13th august 1999

!####
subroutine Mul1(nMat,F,iCre,iAnn,mat,m_Ord,nOsc,rdx)

use mula_global, only: mdim1, ndim1, ndim2
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: m_Ord, nOsc, nMat(0:m_Ord,nOsc), iCre(0:ndim1,ndim2), iAnn(0:ndim1,ndim2)
real(kind=wp), intent(out) :: F(0:mdim1,0:ndim1)
real(kind=wp), intent(in) :: Mat(nOsc), rdx(1)
integer(kind=iwp) :: i, iOrd, iOsc, jOrd
real(kind=wp) :: fact, sqr(0:50)

do i=0,50
  sqr(i) = sqrt(Half*i)
end do
do iOrd=0,m_Ord
  do iOsc=1,nOsc
    jOrd = iAnn(iOrd,iOsc)
    if (jOrd >= 0) then
      Fact = sqr(nmat(iord,iosc))*rdx(1)
      F(iOrd,jOrd) = F(iOrd,jOrd)+Fact*Mat(iOsc)
    end if
  end do
end do
do iOrd=0,m_Ord
  do iOsc=1,nOsc
    jOrd = iCre(iOrd,iOsc)
    if (jOrd >= 0) then
      Fact = sqr(nmat(jord,iosc))
      F(iOrd,jOrd) = F(iOrd,jOrd)+Fact*Mat(iOsc)
    end if
  end do
end do

end subroutine Mul1
!####
subroutine Mul2(nMat,F,iCre,iAnn,mat,m_Ord,nOsc,rdx)

use mula_global, only: mdim1, ndim1, ndim2
use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: m_Ord, nOsc, nMat(0:m_Ord,nOsc), iCre(0:ndim1,ndim2), iAnn(0:ndim1,ndim2)
real(kind=wp), intent(inout) :: F(0:mdim1,0:ndim1)
real(kind=wp), intent(in) :: Mat(nOsc,nOsc), rdx(2)
integer(kind=iwp) :: i, iOrd, iOsc, jOrd, jOsc, kOrd
real(kind=wp) :: fact, r, rsym, sqr(0:50)

rsym = Half
do i=0,50
  sqr(i) = sqrt(Half*i)
end do
r = rdx(1)*rdx(2)
do iOrd=0,m_Ord
  do iOsc=1,nOsc
    jOrd = iAnn(iOrd,iOsc)
    if (jOrd >= 0) then
      do jOsc=1,nOsc
        kOrd = iAnn(jOrd,jOsc)
        if (kOrd >= 0) then
          Fact = sqr(nmat(iord,iosc))*sqr(nmat(jord,josc))*r*rsym
          F(iOrd,kOrd) = F(iOrd,kOrd)+Fact*Mat(iOsc,josc)
        end if
      end do
    end if
  end do
end do
if (r /= Zero) then
  do iOrd=0,m_Ord
    do iOsc=1,nOsc
      jOrd = iann(iord,iosc)
      if (jOrd >= 0) then
        do jOsc=1,nOsc
          kOrd = iCre(jord,josc)
          if (kOrd >= 0) then
            r = rdx(1)*Mat(iosc,jOsc)+rdx(2)*Mat(josc,iOsc)
            Fact = sqr(nmat(iord,iosc))*sqr(nmat(kord,josc))*r*rsym
            F(iOrd,kOrd) = F(iOrd,kOrd)+Fact
          end if
        end do
      end if
    end do
  end do
end if
do iOrd=0,m_Ord
  do iOsc=1,nOsc
    jOrd = iCre(iord,iosc)
    if (jOrd >= 0) then
      do jOsc=1,nOsc
        kOrd = iCre(jord,josc)
        if (kOrd >= 0) then
          Fact = sqr(nmat(jord,iosc))*sqr(nmat(kord,josc))*rsym
          F(iOrd,kOrd) = F(iOrd,kOrd)+Fact*Mat(iosc,jOsc)
        end if
      end do
    end if
  end do
end do
r = rdx(2)
do iOrd=0,m_Ord
  do iOsc=1,nOsc
    F(iOrd,iOrd) = F(iOrd,iOrd)+Mat(iosc,iOsc)*r*rsym*Half
  end do
end do

end subroutine Mul2
!####
subroutine Mul3(nMat,F,iCre,iAnn,mat,m_Ord,nOsc,rdx)

use mula_global, only: mdim1, ndim1, ndim2
use Constants, only: One, Six, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: m_Ord, nOsc, nMat(0:m_Ord,nOsc), iCre(0:ndim1,ndim2), iAnn(0:ndim1,ndim2)
real(kind=wp), intent(inout) :: F(0:mdim1,0:ndim1)
real(kind=wp), intent(in) :: Mat(nOsc,nOsc,nOsc), rdx(3)
integer(kind=iwp) :: i, iOrd, iOsc, jOrd, jOsc, kOrd, kOsc, lOrd
real(kind=wp) :: fact, r, relem, rsym, sqr(0:50)

rsym = One/Six
do i=0,50
  sqr(i) = sqrt(Half*i)
end do
do iOrd=0,m_Ord
  do iOsc=1,nOsc
    jOrd = iAnn(iOrd,iOsc)
    if (jOrd >= 0) then
      do jOsc=1,nOsc
        kOrd = iAnn(jOrd,jOsc)
        if (kOrd >= 0) then
          do kOsc=1,nOsc
            lOrd = iAnn(kOrd,kOsc)
            if (lOrd >= 0) then
              r = rdx(1)*rdx(2)*rdx(3)
              Fact = sqr(nmat(iord,iosc))*sqr(nmat(jord,josc))*sqr(nmat(kord,kosc))*r*rsym
              F(iOrd,lOrd) = F(iOrd,lOrd)+Fact*Mat(iOsc,josc,kosc)
            end if
          end do
        end if
      end do
    end if
  end do
end do

do iOrd=0,m_Ord
  do iOsc=1,nOsc
    jOrd = iAnn(iOrd,iOsc)
    if (jOrd >= 0) then
      do jOsc=1,nOsc
        kOrd = iAnn(jOrd,jOsc)
        if (kOrd >= 0) then
          do kOsc=1,nOsc
            lOrd = iCre(kOrd,kOsc)
            if (lOrd >= 0) then
              r = rdx(1)*rdx(2)*Mat(iOsc,josc,kosc)+rdx(2)*rdx(3)*Mat(kOsc,iosc,josc)+rdx(3)*rdx(1)*Mat(iOsc,kosc,josc)
              Fact = sqr(nmat(iord,iosc))*sqr(nmat(jord,josc))*sqr(nmat(lord,kosc))*r*rsym
              F(iOrd,lOrd) = F(iOrd,lOrd)+Fact
            end if
          end do
        end if
      end do
    end if
  end do
end do

do iOrd=0,m_Ord
  do iOsc=1,nOsc
    jOrd = iAnn(iOrd,iOsc)
    if (jOrd >= 0) then
      do jOsc=1,nOsc
        kOrd = iCre(jOrd,jOsc)
        if (kOrd >= 0) then
          do kOsc=1,nOsc
            lOrd = iCre(kOrd,kOsc)
            if (lOrd >= 0) then
              r = rdx(1)*mat(iosc,josc,kosc)+rdx(2)*mat(josc,iosc,kosc)+rdx(3)*mat(josc,kosc,iosc)
              Fact = sqr(nmat(iord,iosc))*sqr(nmat(kord,josc))*sqr(nmat(lord,kosc))*rsym*r
              F(iOrd,lOrd) = F(iOrd,lOrd)+Fact
            end if
          end do
        end if
      end do
    end if
  end do
end do

do iOrd=0,m_Ord
  do iOsc=1,nOsc
    jOrd = iCre(iOrd,iOsc)
    if (jOrd >= 0) then
      do jOsc=1,nOsc
        kOrd = iCre(jOrd,jOsc)
        if (kOrd >= 0) then
          do kOsc=1,nOsc
            lOrd = iCre(kOrd,kOsc)
            if (lOrd >= 0) then
              Fact = sqr(nmat(jord,iosc))*sqr(nmat(kord,josc))*sqr(nmat(lord,kosc))*rsym
              F(iOrd,lOrd) = F(iOrd,lOrd)+Fact*Mat(iOsc,josc,kosc)
            end if
          end do
        end if
      end do
    end if
  end do
end do

do iOrd=0,m_Ord
  do iOsc=1,nOsc
    jOrd = iCre(iOrd,iOsc)
    if (jOrd >= 0) then
      do jOsc=1,nOsc
        Fact = sqr(nmat(jord,iosc))
        relem = Mat(josc,josc,iosc)*rdx(2)+rdx(3)*(mat(iosc,josc,josc)+mat(josc,iosc,josc))
        F(iOrd,jOrd) = F(iOrd,jOrd)+Fact*relem*rsym*Half
      end do
    end if
  end do
end do

do iOrd=0,m_Ord
  do iOsc=1,nOsc
    jOrd = iann(iOrd,iOsc)
    if (jOrd >= 0) then
      do jOsc=1,nOsc
        Fact = sqr(nmat(iord,iosc))
        relem = Mat(iosc,josc,josc)*rdx(1)*rdx(3)+rdx(2)*rdx(3)*(mat(josc,josc,iosc)+mat(josc,iosc,josc))
        F(iOrd,jOrd) = F(iOrd,jOrd)+Fact*relem*rsym*Half
      end do
    end if
  end do
end do

end subroutine Mul3
!####
subroutine Mul4(nMat,F,iCre,iAnn,mat,m_Ord,nOsc,rdx)

use mula_global, only: mdim1, ndim1, ndim2
use Constants, only: One, Half, Quart
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: m_Ord, nOsc, nMat(0:m_Ord,nOsc), iCre(0:ndim1,ndim2), iAnn(0:ndim1,ndim2)
real(kind=wp), intent(inout) :: F(0:mdim1,0:ndim1)
real(kind=wp), intent(in) :: Mat(nOsc,nOsc,nOsc,nOsc), rdx(4)
integer(kind=iwp) :: i, iOrd, iOsc, jOrd, jOsc, kOrd, kOsc, lOrd, lOsc, mOrd
real(kind=wp) :: fact, r, relem, rsym, sqr(0:50)

do i=0,50
  sqr(i) = sqrt(Half*i)
end do
rsym = One/24.0_wp
r = Rdx(1)*rdx(2)*rdx(3)*rdx(4)
do iOrd=0,m_Ord
  do iOsc=1,nOsc
    jOrd = iAnn(iOrd,iOsc)
    if (jOrd >= 0) then
      do jOsc=1,nOsc
        kOrd = iAnn(jOrd,jOsc)
        if (kOrd >= 0) then
          do kOsc=1,nOsc
            lOrd = iAnn(kOrd,kOsc)
            if (lOrd >= 0) then
              do lOsc=1,nOsc
                mOrd = iAnn(lOrd,lOsc)
                if (mOrd >= 0) then
                  Fact = sqr(nmat(iord,iosc))*sqr(nmat(jord,josc))*sqr(nmat(kord,kosc))*sqr(nmat(lord,losc))*r*rsym
                  F(iOrd,mOrd) = F(iOrd,mOrd)+Fact*Mat(iOsc,josc,kosc,losc)
                end if
              end do
            end if
          end do
        end if
      end do
    end if
  end do
end do
do iOrd=0,m_Ord
  do iOsc=1,nOsc
    jOrd = iAnn(iOrd,iOsc)
    if (jOrd >= 0) then
      do jOsc=1,nOsc
        kOrd = iAnn(jOrd,jOsc)
        if (kOrd >= 0) then
          do kOsc=1,nOsc
            lOrd = iAnn(kOrd,kOsc)
            if (lOrd >= 0) then
              do lOsc=1,nOsc
                mOrd = iCre(lOrd,lOsc)
                if (mOrd >= 0) then
                  r = rdx(1)*rdx(2)*rdx(3)*mat(iosc,josc,kosc,losc)+rdx(2)*rdx(3)*rdx(4)*mat(losc,iosc,josc,kosc)+ &
                      rdx(3)*rdx(4)*rdx(1)*mat(iosc,losc,josc,kosc)+rdx(4)*rdx(1)*rdx(2)*mat(iosc,josc,losc,kosc)
                  Fact = sqr(nmat(iord,iosc))*sqr(nmat(jord,josc))*sqr(nmat(kord,kosc))*sqr(nmat(mord,losc))*r*rsym
                  F(iOrd,mOrd) = F(iOrd,mOrd)+Fact
                end if
              end do
            end if
          end do
        end if
      end do
    end if
  end do
end do
do iOrd=0,m_Ord
  do iOsc=1,nOsc
    jOrd = iAnn(iOrd,iOsc)
    if (jOrd >= 0) then
      do jOsc=1,nOsc
        kOrd = iAnn(jOrd,jOsc)
        if (kOrd >= 0) then
          do kOsc=1,nOsc
            lOrd = iCre(kOrd,kOsc)
            if (lOrd >= 0) then
              do lOsc=1,nOsc
                mOrd = iCre(lOrd,lOsc)
                if (mOrd >= 0) then
                  r = rdx(1)*rdx(2)*mat(iosc,josc,kosc,losc)+rdx(2)*rdx(3)*mat(kosc,iosc,josc,losc)+ &
                      rdx(3)*rdx(4)*mat(kosc,losc,iosc,josc)+rdx(4)*rdx(1)*mat(iosc,kosc,losc,josc)+ &
                      rdx(1)*rdx(3)*mat(iosc,kosc,josc,losc)+rdx(4)*rdx(2)*mat(kosc,iosc,losc,josc)
                  Fact = sqr(nmat(iord,iosc))*sqr(nmat(jord,josc))*sqr(nmat(lord,kosc))*sqr(nmat(mord,losc))*r*rsym
                  F(iOrd,mOrd) = F(iOrd,mOrd)+Fact
                end if
              end do
            end if
          end do
        end if
      end do
    end if
  end do
end do
do iOrd=0,m_Ord
  do iOsc=1,nOsc
    jOrd = iAnn(iOrd,iOsc)
    if (jOrd >= 0) then
      do jOsc=1,nOsc
        kOrd = iCre(jOrd,jOsc)
        if (kOrd >= 0) then
          do kOsc=1,nOsc
            lOrd = icre(kOrd,kOsc)
            if (lOrd >= 0) then
              do lOsc=1,nOsc
                mOrd = iCre(lOrd,lOsc)
                if (mOrd >= 0) then
                  r = rdx(1)*mat(iosc,josc,kosc,losc)+rdx(2)*mat(josc,iosc,kosc,losc)+rdx(3)*mat(josc,kosc,iosc,losc)+ &
                      rdx(4)*mat(josc,kosc,losc,iosc)
                  Fact = sqr(nmat(iord,iosc))*sqr(nmat(kord,josc))*sqr(nmat(lord,kosc))*sqr(nmat(mord,losc))*r*rsym
                  F(iOrd,mOrd) = F(iOrd,mOrd)+Fact
                end if
              end do
            end if
          end do
        end if
      end do
    end if
  end do
end do
do iOrd=0,m_Ord
  do iOsc=1,nOsc
    jOrd = iCre(iOrd,iOsc)
    if (jOrd >= 0) then
      do jOsc=1,nOsc
        kOrd = iCre(jOrd,jOsc)
        if (kOrd >= 0) then
          do kOsc=1,nOsc
            lOrd = iCre(kOrd,kOsc)
            if (lOrd >= 0) then
              do lOsc=1,nOsc
                mOrd = iCre(lOrd,lOsc)
                if (mOrd >= 0) then
                  Fact = sqr(nmat(jord,iosc))*sqr(nmat(kord,josc))*sqr(nmat(lord,kosc))*sqr(nmat(mord,losc))*rsym
                  F(iOrd,mOrd) = F(iOrd,mOrd)+Fact*Mat(iOsc,josc,kosc,losc)
                end if
              end do
            end if
          end do
        end if
      end do
    end if
  end do
end do
do iOrd=0,m_Ord
  do iOsc=1,nOsc
    jOrd = iann(iOrd,iOsc)
    if (jOrd >= 0) then
      do jOsc=1,nOsc
        kOrd = iann(jOrd,jOsc)
        if (kOrd >= 0) then
          do kOsc=1,nOsc
            Fact = sqr(nmat(iord,iosc))*sqr(nmat(jord,josc))*rsym
            relem = rdx(1)*rdx(2)*rdx(4)*mat(iosc,josc,kosc,kosc)
            relem = relem+rdx(1)*rdx(3)*rdx(4)*(mat(iosc,kosc,kosc,josc)+mat(iosc,kosc,josc,kosc))
            relem = relem+rdx(2)*rdx(3)*rdx(4)*(mat(kosc,kosc,josc,iosc)+mat(kosc,iosc,kosc,josc)+mat(kosc,josc,iosc,kosc))
            F(iOrd,kOrd) = F(iOrd,kOrd)+Fact*relem*Half
          end do
        end if
      end do
    end if
  end do
end do

do iOrd=0,m_Ord
  do iOsc=1,nOsc
    jOrd = icre(iOrd,iOsc)
    if (jOrd >= 0) then
      do jOsc=1,nOsc
        kOrd = icre(jOrd,jOsc)
        if (kOrd >= 0) then
          do kOsc=1,nOsc
            Fact = sqr(nmat(jord,iosc))*sqr(nmat(kord,josc))*rsym
            relem = rdx(2)*mat(kosc,kosc,iosc,josc)
            relem = relem+rdx(3)*(mat(iosc,kosc,kosc,josc)+mat(kosc,iosc,kosc,josc))
            relem = relem+rdx(4)*(mat(iosc,josc,kosc,kosc)+mat(iosc,kosc,josc,kosc)+mat(kosc,iosc,josc,kosc))
            F(iOrd,kOrd) = F(iOrd,kOrd)+Fact*relem*Half
          end do
        end if
      end do
    end if
  end do
end do
do iOrd=0,m_Ord
  do iOsc=1,nOsc
    jOrd = iann(iOrd,iOsc)
    if (jOrd >= 0) then
      do jOsc=1,nOsc
        kOrd = iCre(jOrd,jOsc)
        if (kOrd >= 0) then
          do kOsc=1,nOsc
            Fact = sqr(nmat(iord,iosc))*sqr(nmat(kord,josc))*rsym
            relem = (Mat(jOsc,kosc,kosc,iosc)+Mat(kOsc,josc,kosc,iosc)+Mat(jOsc,kosc,iosc,kosc)+ &
                    Mat(kOsc,josc,iosc,kosc))*rdx(3)*rdx(4)
            relem = relem+(Mat(kOsc,kosc,josc,iosc)+Mat(jOsc,iosc,kosc,kosc)+Mat(kOsc,iosc,josc,kosc))*rdx(2)*rdx(4)
            relem = relem+(Mat(kOsc,kosc,iosc,josc)+Mat(kOsc,iosc,kosc,josc))*rdx(2)*rdx(3)
            relem = relem+(Mat(iOsc,kosc,josc,kosc)+Mat(iOsc,josc,kosc,kosc))*rdx(1)*rdx(4)
            relem = relem+Mat(iOsc,kosc,kosc,josc)*rdx(1)*rdx(3)
            F(iOrd,kOrd) = F(iOrd,kOrd)+Fact*relem*Half
          end do
        end if
      end do
    end if
  end do
end do
do iOrd=0,m_Ord
  do iOsc=1,nOsc
    do jOsc=1,nOsc
      relem = rdx(4)*rdx(3)*(mat(iosc,josc,josc,iosc)+mat(iosc,josc,iosc,josc))
      relem = relem+rdx(2)*rdx(4)*mat(iosc,iosc,josc,josc)
      F(iOrd,iOrd) = F(iOrd,iOrd)+relem*rsym*Quart
    end do
  end do
end do

end subroutine Mul4
