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

subroutine RTSWGH(TARR,NT,U2,WGH,NRYS,nOrdOp)

use vRys_RW, only: HerR2, HerW2, iHerR2, iHerW2
use abdata, only: atab, btab, p0, tvalue
use Gateway_global, only: asymptotic_Rys
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Five, Twelve, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NT, NRYS, nOrdOp
real(kind=wp), intent(in) :: TARR(NT)
real(kind=wp), intent(out) :: U2(NRYS,NT), WGH(NRYS,NT)
integer(kind=iwp) :: IDEG, iroot, IT, J, k, nx
real(kind=wp) :: a2, a3, a4, a5, a6, b1, b2, b3, b4, b5, BK, c1, c2, c3, c4, c5, c6, corr, DELTA, p, R, R1, R2, RSUM, T, tmp, x1, &
                 x2, x3, xn, Z, ZZ
real(kind=wp), allocatable :: ALPHA(:), BETA(:), BINV(:), ROOT(:,:), RYS(:), RYSD(:)
real(kind=wp), parameter :: coef1 = -One/120.0_wp, coef2 = Five/120.0_wp, coef3 = -One/Twelve, coef4 = One/Twelve, &
                            coef5 = -Five/120.0_wp, coef6 = One/120.0_wp

if (NRYS > ubound(atab,1)) then
  call WarningMessage(2,' Too many requested Rys roots.')
  call AbEnd()
end if

call mma_allocate(ALPHA,[0,max(NRYS,1)],label='ALPHA')
call mma_allocate(BETA,[0,max(NRYS,1)],label='BETA')
call mma_allocate(BINV,max(NRYS,2),label='BINV')
call mma_allocate(ROOT,max(NRYS,2),max(NRYS,2),label='ROOT')
call mma_allocate(RYS,[0,max(NRYS,2)],label='RYS')
call mma_allocate(RYSD,[0,max(NRYS,1)],label='RYSD')
RYSD(0) = Zero

do IT=1,NT
  T = TARR(IT)
  ! Use asymptotic formulae if outside table.
  !MAW if (t > TVALUE(ubound(TVALUE,1)-2)) then

  ! For the FMM we use the asymptotic limit to compute the
  ! multipole-component of the integrals

  if ((t > TVALUE(ubound(TVALUE,1)-2)) .or. asymptotic_Rys) then
    tmp = One/T
    U2(:,IT) = HerR2(iHerR2(nRys):iHerR2(nRys)+nRys-1)*tmp
    WGH(:,IT) = HerW2(iHerW2(nRys):iHerW2(nRys)+nRys-1)*sqrt(tmp)
    cycle
  end if
  ! translate to tabulation function for equidist. interp.
  ! xn=interpol. variable.
  ! Ex: T=0.0--0.05 gives xn=0.0--1.0 (approx.)
  xn = Five*T+200.0_wp*T/(14.0_wp+T)
  nx = int(xn)
  p = xn-real(nx,kind=wp)
  a2 = (p+Two)
  a3 = (p+One)*a2
  a4 = (p)*a3
  a5 = (p-One)*a4
  a6 = (p-Two)*a5
  b5 = (p-Three)
  b4 = (p-Two)*b5
  b3 = (p-One)*b4
  b2 = (p)*b3
  b1 = (p+One)*b2
  c1 = coef1*b1
  c2 = coef2*a2*b2
  c3 = coef3*a3*b3
  c4 = coef4*a4*b4
  c5 = coef5*a5*b5
  c6 = coef6*a6
  ALPHA(0) = c1*ATAB(0,nx-2)+c2*ATAB(0,nx-1)+c3*ATAB(0,nx)+c4*ATAB(0,nx+1)+c5*ATAB(0,nx+2)+c6*ATAB(0,nx+3)
  ALPHA(1:NRYS) = c1*ATAB(1:NRYS,nx-2)+c2*ATAB(1:NRYS,nx-1)+c3*ATAB(1:NRYS,nx)+c4*ATAB(1:NRYS,nx+1)+c5*ATAB(1:NRYS,nx+2)+ &
                  c6*ATAB(1:NRYS,nx+3)
  BETA(1:NRYS) = c1*BTAB(1:NRYS,nx-2)+c2*BTAB(1:NRYS,nx-1)+c3*BTAB(1:NRYS,nx)+c4*BTAB(1:NRYS,nx+1)+c5*BTAB(1:NRYS,nx+2)+ &
                 c6*BTAB(1:NRYS,nx+3)
  BINV(1:NRYS) = One/BETA(1:NRYS)
  rys(0) = c1*p0(nx-2)+c2*p0(nx-1)+c3*p0(nx)+c4*p0(nx+1)+c5*p0(nx+2)+c6*p0(nx+3)
  ROOT(1,1) = ALPHA(0)
  x1 = (ALPHA(0)+ALPHA(1))*Half
  x2 = (ALPHA(0)-ALPHA(1))*Half
  x3 = sqrt(x2**2+BETA(1)**2)
  ROOT(1,2) = x1-x3
  ROOT(2,2) = x1+x3
  ! LOOP OVER DEGREE OF RYS POLY
  do IDEG=3,NRYS
    ! ESTIMATE POSITION OF U2 OF THIS DEGREE:
    ROOT(1,IDEG) = (real(IDEG,kind=wp)-Half)*ROOT(1,IDEG-1)/real(IDEG,kind=wp)
    ROOT(IDEG,IDEG) = ONE-(real(IDEG,kind=wp)-Half)*(ONE-ROOT(IDEG-1,IDEG-1))/real(IDEG,kind=wp)
    do IROOT=2,IDEG-1
      R1 = ROOT(IROOT,IDEG-1)
      R2 = ROOT(IROOT-1,IDEG-1)
      x2 = (real(IROOT,kind=wp)-Half)/real(IDEG,kind=wp)
      x1 = ONE-x2
      R = x1*R1+x2*R2
      ROOT(IROOT,IDEG) = R
    end do
    !if (IDEG == NRYS) ITER = 0
    RYSD(1) = RYS(0)*BINV(1)
    do IROOT=1,IDEG
      Z = ROOT(IROOT,IDEG)
      ! Compute the correction coefficient:
      corr = Zero
      do J=1,iroot-1
        corr = corr+one/(root(iroot,ideg)-root(j,ideg))
      end do
      do J=iroot+1,ideg
        corr = corr+one/(root(iroot,ideg)-root(j,ideg))
      end do
      ! COMPUTE RYS AND FIRST DERIVATIVE, DO NEWTON-RAPHSON:
      do
        RYS(1) = (Z-ALPHA(0))*RYSD(1)
        ZZ = (Z-ALPHA(1))
        RYSD(2) = (ZZ*RYSD(1)+RYS(1))*BINV(2)
        RYS(2) = (ZZ*RYS(1)-BETA(1)*RYS(0))*BINV(2)
        do K=2,IDEG-1
          ZZ = Z-ALPHA(K)
          BK = BETA(K)
          RYSD(K+1) = (RYS(K)+ZZ*RYSD(K)-BK*RYSD(K-1))*BINV(K+1)
          RYS(K+1) = (ZZ*RYS(K)-BK*RYS(K-1))*BINV(K+1)
        end do
        DELTA = -RYS(IDEG)/(RYSD(IDEG)-CORR*RYS(IDEG))
        Z = Z+DELTA
        !if (IDEG == NRYS) ITER = ITER+1
        if (abs(DELTA) <= 1.0e-8_wp) exit
      end do
      ROOT(IROOT,IDEG) = Z
    end do
  end do
  !if (NRYS > 2) write(u6,'(1x,a,f8.2)') ' Avg. iter/root:',(iter*one)/NRYS
  do IROOT=1,NRYS
    Z = ROOT(IROOT,NRYS)
    ! COMPUTE RYS VALUES AND ADD SQUARES TO GET WGH:
    RSUM = RYS(0)**2
    if (NRYS /= 1) then
      RYS(1) = (Z-ALPHA(0))*RYS(0)*BINV(1)
      RSUM = RSUM+RYS(1)**2
      if (NRYS /= 2) then
        ZZ = (Z-ALPHA(1))
        RYS(2) = (ZZ*RYS(1)-BETA(1)*RYS(0))*BINV(2)
        RSUM = RSUM+RYS(2)**2
        if (NRYS /= 3) then
          do K=2,NRYS-2
            ZZ = Z-ALPHA(K)
            BK = BETA(K)
            RYS(K+1) = (ZZ*RYS(K)-BK*RYS(K-1))*BINV(K+1)
            RSUM = RSUM+RYS(K+1)**2
          end do
        end if
      end if
    end if
    WGH(IROOT,IT) = ONE/RSUM
    U2(IROOT,IT) = ROOT(IROOT,NRYS)
  end do

end do

call mma_deallocate(ALPHA)
call mma_deallocate(BETA)
call mma_deallocate(BINV)
call mma_deallocate(ROOT)
call mma_deallocate(RYS)
call mma_deallocate(RYSD)

if ((nOrdOp == 1) .or. (nOrdOp == 2)) WGH(:,:) = (U2(:,:)/(One-U2(:,:)))**nOrdOp*WGH(:,:)

#ifdef _DEBUGPRINT_
call RecPrt(' RTSWGH: Roots',' ',U2,nRys,nT)
call RecPrt(' RTSWGH: Weights',' ',Wgh,nRys,nT)
#endif

return

end subroutine RTSWGH
