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

use Constants, only: Zero, Half
use rctfld_module, only: MaxA, nSparse, MaxB, lSparse, Scala, Scalc, nGrid_Eff, LatAto, RadLat, lRFCav, rds, DieDel, rSca, &
                         DistSparse, nExpO, PreFac, Polsi, Dipsi, Cordsi, MaxC, RotAlpha, RotBeta, RotGamma

implicit none
integer nGrid_, MaxAto, nPolComp, nXF, nOrd_XF, iXPolType
real*8 Grid(3,nGrid_), PolEff(nPolComp,nGrid_), DipEff(nGrid_)
real*8 cord(3,maxato), atorad(maxato), XF(*)
integer XEle(nXF)
integer ixyz, nElem
integer Inc, iOrdOp
real*8 tr(3,3), co(3)
integer ii, jj, kk, ni, nj, nk, nGridOld, l, m, ixf, i, j, k, nGrid_Eff_
real*8 xs, ys, zs, xp, yp, zp, rp2, rp, Ener1, Ener, xa, ya, za, drp, rrr, dGGX, dGGY, dGGZ, rpa2, atrad, fac
real*8, external :: CovRadT
! Statement function for Cartesian index
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

! Calculate number of entries per XFIELD point
Inc = 3
do iOrdOp=0,nOrd_XF
  Inc = Inc+nElem(iOrdOp)
end do
if (iXPolType > 0) Inc = Inc+6

!write(6,*) 'lattcr: polsi,dipsi=',polsi,dipsi

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
        xs = (dble(ii)+(dble(nSparse-1)*half))*scala
        ys = (dble(jj)+(dble(nSparse-1)*half))*scala
        zs = (dble(kk)+(dble(nSparse-1)*half))*scalc
        nGridOld = nGrid_Eff
        do l=1,latato
          co(1) = xs+cordsi(1,l)
          co(2) = ys+cordsi(2,l)
          co(3) = zs+cordsi(3,l)
          xp = Zero
          yp = Zero
          zp = Zero
          do m=1,3
            xp = xp+co(m)*tr(1,m)
            yp = yp+co(m)*tr(2,m)
            zp = zp+co(m)*tr(3,m)
          end do
          rp2 = xp*xp+yp*yp+zp*zp
          if (rp2 > radlat**2) go to 13
          if (lRFCav) then
            rp = sqrt(rp2)
            drp = rds-rp
            if (drp <= Zero) go to 13
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
            if (rpa2 < distSparse**2) goto 13
            ener = ener+(rrr/rpa2)**nexpo
          end do

          ! Check if the XFIELD multipoles annhilates the grid point

          do iXF=1,nXF
            xa = XF((iXF-1)*Inc+1)
            ya = XF((iXF-1)*Inc+2)
            za = XF((iXF-1)*Inc+3)
            if (XEle(iXF) <= 0) then
              atrad = -dble(XEle(iXF))/1000.0d0
            else
              atrad = CovRadT(XEle(iXF))
            end if
            rrr = (atrad*rsca)**2
            dggx = xa-xp
            dggy = ya-yp
            dggz = za-zp
            rpa2 = dggx**2+dggy**2+dggz**2
            if (rpa2 < distSparse**2) goto 13
            ener = ener+(rrr/rpa2)**nexpo
            !if (rpa2 < 6.0D0) write(6,*) 'DIST',iGrid,iXF,sqrt(rpa2),atrad
          end do

          ener = prefac*ener*500.0d0+ener1
          if (ener > 12.d0) goto 13
          fac = exp(-ener)
          nGrid_Eff_ = nGrid_Eff_+1
          Grid(1,nGrid_Eff_) = xp
          Grid(2,nGrid_Eff_) = yp
          Grid(3,nGrid_Eff_) = zp
          PolEff(1,nGrid_Eff_) = polsi*dble(nSparse)**3.0d0*fac*fac
          DipEff(nGrid_Eff_) = dipsi*dble(nSparse)**1.5d0*fac
          write(6,*) 'DGRID',xp,yp,zp,fac
        end do
        ! Every sparse lattice point in this cell is ok, so skip dense grid

        goto 14

13      continue
        ! Not every sparse lattice point is ok, so delete the sparse and do dense
        nGrid_Eff_ = nGridOld
      end if
      ! Start the normal, dense grid
      do i=0,ni-1
        do j=0,nj-1
          do k=0,nk-1
            xs = dble(ii+i)*scala
            ys = dble(jj+j)*scala
            zs = dble(kk+k)*scalc
            do l=1,latato
              co(1) = xs+cordsi(1,l)
              co(2) = ys+cordsi(2,l)
              co(3) = zs+cordsi(3,l)
              xp = Zero
              yp = Zero
              zp = Zero
              do m=1,3
                xp = xp+co(m)*tr(1,m)
                yp = yp+co(m)*tr(2,m)
                zp = zp+co(m)*tr(3,m)
              end do
              rp2 = xp*xp+yp*yp+zp*zp
              if (rp2 > radlat**2) go to 11
              if (lRFCav) then
                rp = sqrt(rp2)
                drp = rds-rp
                if (drp <= Zero) go to 11
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
                !if (rpa2 < 6.0D0) write(6,*) 'DIST0',iGrid,sqrt(rpa2),atorad(m)
              end do

              ! Check if the XFIELD multipoles annhilates the grid point

              do iXF=1,nXF
                xa = XF((iXF-1)*Inc+1)
                ya = XF((iXF-1)*Inc+2)
                za = XF((iXF-1)*Inc+3)
                if (XEle(iXF) <= 0) then
                  atrad = -dble(XEle(iXF))/1000.0d0
                else
                  atrad = CovRadT(XEle(iXF))
                end if
                rrr = (atrad*rsca)**2
                dggx = xa-xp
                dggy = ya-yp
                dggz = za-zp
                rpa2 = dggx**2+dggy**2+dggz**2
                ener = ener+(rrr/rpa2)**nexpo
                !if (rpa2 < 6.0D0) write(6,*) 'DIST',iGrid,iXF,sqrt(rpa2),atrad
              end do

              !ener = prefac*ener*tk5+ener1
              ener = prefac*ener*500.0d0+ener1
              !ener = prefac*ener*1000.0D0+ener1
              if (ener > 12.d0) then
                write(6,*) 'REMOVED',xp,yp,zp
                Go To 11
              end if

              fac = exp(-ener)

              nGrid_Eff_ = nGrid_Eff_+1
              Grid(1,nGrid_Eff_) = xp
              Grid(2,nGrid_Eff_) = yp
              Grid(3,nGrid_Eff_) = zp
              PolEff(1,nGrid_Eff_) = polsi*fac*fac
              DipEff(nGrid_Eff_) = dipsi*fac
              write(6,*) 'GRID',xp,yp,zp,fac
11            continue
            end do  ! l
          end do  ! k
        end do  ! j
      end do  ! i
14    continue
    end do  ! kk
  end do  ! jj
end do  ! ii

!write(6,*) 'Grid...',Grid(1,1),Grid(2,1),Grid(3,1)
!write(6,*) 'DipEff...',DipEff(1),DipEff(101)
!write(6,*) 'PolEff...',PolEff(1),PolEff(101)
!write(6,*) 'Lattcr: nGrid_Eff=',nGrid_Eff

return

end subroutine lattcr
