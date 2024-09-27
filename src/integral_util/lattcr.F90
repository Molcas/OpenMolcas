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
! Copyright (C) 2000, Gunnar Karlstrom                                 *
!               2000, Roland Lindh                                     *
!***********************************************************************

subroutine lattcr(Grid,nGrid_,nGrid_Eff_,PolEff,DipEff,cord,maxato,atorad,nPolComp,XF,nXF,nOrd_XF,XEle,iXPolType)
!***********************************************************************
!                                                                      *
!     Object: to compute effective polarizabilities and dipole moments *
!             on the Langevin grid.                                    *
!                                                                      *
!     Authors: G. Karlstroem                                           *
!              Dept. of Theor. Chem., Univ. of Lund, Sweden.           *
!                                                                      *
!              and                                                     *
!                                                                      *
!              R. Lindh                                                *
!              Dept. of Chem. Phys., Univ. of Lund, Sweden.            *
!                                                                      *
!              March 2000                                              *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use rctfld_module, only: Cordsi, DieDel, Dipsi, DistSparse, LatAto, lRFCav, lSparse, MaxA, MaxB, MaxC, nExpO, nGrid_Eff, nSparse, &
                         Polsi, PreFac, RadLat, rds, RotAlpha, RotBeta, RotGamma, rSca, Scala, Scalc
use Constants, only: Zero, Three, Twelve, Half, OneHalf
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nGrid_, MaxAto, nPolComp, nXF, nOrd_XF, XEle(nXF), iXPolType
real(kind=wp), intent(inout) :: Grid(3,nGrid_), PolEff(nPolComp,nGrid_), DipEff(nGrid_)
integer(kind=iwp), intent(inout) :: nGrid_Eff_
real(kind=wp), intent(in) :: cord(3,maxato), atorad(maxato), XF(*)
integer(kind=iwp) :: i, ii, Inc, iOrdOp, ixf, j, jj, k, kk, l, m, nGridOld, ni, nj, nk
real(kind=wp) :: atrad, co(3), dGGX, dGGY, dGGZ, drp, Ener, Ener1, fac, rp, rp2, rpa2, rrr, tr(3,3), xa, xp, xs, ya, yp, ys, za, &
                 zp, zs
real(kind=wp), external :: CovRadT

! Calculate number of entries per XFIELD point
Inc = 3
do iOrdOp=0,nOrd_XF
  Inc = Inc+nTri_Elem1(iOrdOp)
end do
if (iXPolType > 0) Inc = Inc+6

!write(u6,*) 'lattcr: polsi,dipsi=',polsi,dipsi

! Rotation matrix for the grid
tr(1,1) = cos(rotGamma)*cos(rotBeta)*cos(rotAlpha)-sin(rotGamma)*sin(rotAlpha)
tr(1,2) = cos(rotGamma)*cos(rotBeta)*sin(rotAlpha)+sin(rotGamma)*cos(rotAlpha)
tr(1,3) = -cos(rotGamma)*sin(rotBeta)
tr(2,1) = -sin(rotGamma)*cos(rotBeta)*cos(rotAlpha)-cos(rotGamma)*sin(rotAlpha)
tr(2,2) = -sin(rotGamma)*cos(rotBeta)*sin(rotAlpha)+cos(rotGamma)*cos(rotAlpha)
tr(2,3) = sin(rotGamma)*sin(rotBeta)
tr(3,1) = sin(rotBeta)*cos(rotAlpha)
tr(3,2) = sin(rotBeta)*sin(rotAlpha)
tr(3,3) = cos(rotBeta)

do ii=-(maxa+1),maxa,nSparse
  do jj=-(maxb+1),maxb,nSparse
    do kk=-(maxc+1),maxc,nSparse
      ni = min(nSparse,maxa-ii+1)
      nj = min(nSparse,maxb-jj+1)
      nk = min(nSparse,maxc-kk+1)
      if (LSparse .and. (ni == nSparse) .and. (nj == nSparse) .and. (nk == nSparse)) then
        xs = (real(ii,kind=wp)+(real(nSparse-1,kind=wp)*Half))*scala
        ys = (real(jj,kind=wp)+(real(nSparse-1,kind=wp)*Half))*scala
        zs = (real(kk,kind=wp)+(real(nSparse-1,kind=wp)*Half))*scalc
        nGridOld = nGrid_Eff
        outer: do l=1,latato
          co(1) = xs+cordsi(1,l)
          co(2) = ys+cordsi(2,l)
          co(3) = zs+cordsi(3,l)
          xp = sum(co(:)*tr(1,:))
          yp = sum(co(:)*tr(2,:))
          zp = sum(co(:)*tr(3,:))
          rp2 = xp*xp+yp*yp+zp*zp
          if (rp2 > radlat**2) exit outer
          if (lRFCav) then
            rp = sqrt(rp2)
            drp = rds-rp
            if (drp <= Zero) exit outer
            ener1 = (diedel/drp)**2
          else
            ener1 = Zero
          end if
          ener = Zero

          ! Check if the QC system annhilates the grid point

          do m=1,maxato
            xa = cord(1,m)
            ya = cord(2,m)
            za = cord(3,m)
            rrr = (atorad(m)*rsca)**2
            dggx = xa-xp
            dggy = ya-yp
            dggz = za-zp
            rpa2 = dggx**2+dggy**2+dggz**2
            if (rpa2 < distSparse**2) exit outer
            ener = ener+(rrr/rpa2)**nexpo
          end do

          ! Check if the XFIELD multipoles annhilates the grid point

          do iXF=1,nXF
            xa = XF((iXF-1)*Inc+1)
            ya = XF((iXF-1)*Inc+2)
            za = XF((iXF-1)*Inc+3)
            if (XEle(iXF) <= 0) then
              atrad = -real(XEle(iXF),kind=wp)/1000.0_wp
            else
              atrad = CovRadT(XEle(iXF))
            end if
            rrr = (atrad*rsca)**2
            dggx = xa-xp
            dggy = ya-yp
            dggz = za-zp
            rpa2 = dggx**2+dggy**2+dggz**2
            if (rpa2 < distSparse**2) exit outer
            ener = ener+(rrr/rpa2)**nexpo
            !if (rpa2 < Six) write(u6,*) 'DIST',iGrid,iXF,sqrt(rpa2),atrad
          end do

          ener = prefac*ener*500.0_wp+ener1
          if (ener > Twelve) exit outer
          fac = exp(-ener)
          nGrid_Eff_ = nGrid_Eff_+1
          Grid(1,nGrid_Eff_) = xp
          Grid(2,nGrid_Eff_) = yp
          Grid(3,nGrid_Eff_) = zp
          PolEff(1,nGrid_Eff_) = polsi*real(nSparse,kind=wp)**Three*fac*fac
          DipEff(nGrid_Eff_) = dipsi*real(nSparse,kind=wp)**OneHalf*fac
          write(u6,*) 'DGRID',xp,yp,zp,fac
        end do outer
        ! Every sparse lattice point in this cell is ok, so skip dense grid

        if (l > latato) cycle

        ! Not every sparse lattice point is ok, so delete the sparse and do dense
        nGrid_Eff_ = nGridOld
      end if
      ! Start the normal, dense grid
      do i=0,ni-1
        do j=0,nj-1
          do k=0,nk-1
            xs = real(ii+i,kind=wp)*scala
            ys = real(jj+j,kind=wp)*scala
            zs = real(kk+k,kind=wp)*scalc
            do l=1,latato
              co(1) = xs+cordsi(1,l)
              co(2) = ys+cordsi(2,l)
              co(3) = zs+cordsi(3,l)
              xp = sum(co(:)*tr(1,:))
              yp = sum(co(:)*tr(2,:))
              zp = sum(co(:)*tr(3,:))
              rp2 = xp*xp+yp*yp+zp*zp
              if (rp2 > radlat**2) cycle
              if (lRFCav) then
                rp = sqrt(rp2)
                drp = rds-rp
                if (drp <= Zero) cycle
                ener1 = (diedel/drp)**2
              else
                ener1 = Zero
              end if
              ener = Zero

              ! Check if the QC system annihilates the grid point

              do m=1,maxato
                xa = cord(1,m)
                ya = cord(2,m)
                za = cord(3,m)
                rrr = (atorad(m)*rsca)**2
                dggx = xa-xp
                dggy = ya-yp
                dggz = za-zp
                rpa2 = dggx**2+dggy**2+dggz**2
                ener = ener+(rrr/rpa2)**nexpo
                !if (rpa2 < Six) write(u6,*) 'DIST0',iGrid,sqrt(rpa2),atorad(m)
              end do

              ! Check if the XFIELD multipoles annhilates the grid point

              do iXF=1,nXF
                xa = XF((iXF-1)*Inc+1)
                ya = XF((iXF-1)*Inc+2)
                za = XF((iXF-1)*Inc+3)
                if (XEle(iXF) <= 0) then
                  atrad = -real(XEle(iXF),kind=wp)/1000.0_wp
                else
                  atrad = CovRadT(XEle(iXF))
                end if
                rrr = (atrad*rsca)**2
                dggx = xa-xp
                dggy = ya-yp
                dggz = za-zp
                rpa2 = dggx**2+dggy**2+dggz**2
                ener = ener+(rrr/rpa2)**nexpo
                !if (rpa2 < Six) write(u6,*) 'DIST',iGrid,iXF,sqrt(rpa2),atrad
              end do

              !ener = prefac*ener*tk5+ener1
              ener = prefac*ener*500.0_wp+ener1
              !ener = prefac*ener*1000.0_wp+ener1
              if (ener > Twelve) then
                write(u6,*) 'REMOVED',xp,yp,zp
                cycle
              end if

              fac = exp(-ener)

              nGrid_Eff_ = nGrid_Eff_+1
              Grid(1,nGrid_Eff_) = xp
              Grid(2,nGrid_Eff_) = yp
              Grid(3,nGrid_Eff_) = zp
              PolEff(1,nGrid_Eff_) = polsi*fac*fac
              DipEff(nGrid_Eff_) = dipsi*fac
              write(u6,*) 'GRID',xp,yp,zp,fac
            end do  ! l
          end do  ! k
        end do  ! j
      end do  ! i
    end do  ! kk
  end do  ! jj
end do  ! ii

!write(u6,*) 'Grid...',Grid(1,1),Grid(2,1),Grid(3,1)
!write(u6,*) 'DipEff...',DipEff(1),DipEff(101)
!write(u6,*) 'PolEff...',PolEff(1),PolEff(101)
!write(u6,*) 'Lattcr: nGrid_Eff=',nGrid_Eff

return

end subroutine lattcr
