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

subroutine edip(EF,DipMom,dEF,PolEff,DipEff,Grid,nGrid_Eff,nPolComp,nAnisopol,nXF,iXPolType,nXMolnr,XMolnr)
!***********************************************************************
!                                                                      *
!     Object: to solve equation system iteratively.                    *
!                                                                      *
!     Input:                                                           *
!            dEF   : the electric field of the QM system               *
!            Cavxyz: the MM expansion of the QM system                 *
!            DipEff: Effective dipole moments                          *
!            PolEff: Effective polarizabilities                        *
!            Grid  : list of grid centers                              *
!            nGird_: effective list length                             *
!                                                                      *
!     Scratch:                                                         *
!            Ravxyz: incremental charge distribution on the boundary   *
!                    of the cavity                                     *
!                                                                      *
!     Output:                                                          *
!            EF    : Total EF                                          *
!            DipMom: Langevin dipole moments on the grid               *
!                                                                      *
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

use rctfld_module, only: cLim, DampIter, DipCutOff, EPS, EPSInF, FMax, lAmberPol, lDamping, lMax, lRFCav, rDS, Scal14, TK
use Langevin_arrays, only: Cavxyz, Ravxyz
use Constants, only: Zero, One, Three, Four, Six
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use rctfld_module, only: AFac, ScalA, ScalB, ScalC
use Constants, only: Half
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nGrid_Eff, nPolComp, nAnisoPol, nXF, iXPolType, nXMolNr, XMolnr(nXMolnr,nXF)
real(kind=wp), intent(inout) :: EF(4,nGrid_Eff), DipMom(3,nGrid_Eff), dEF(4,nGrid_Eff)
real(kind=wp), intent(in) :: PolEff(nPolComp,nGrid_Eff), DipEff(nGrid_Eff), Grid(3,nGrid_Eff)
integer(kind=iwp) :: i, iGrid, Iter, jGrid
real(kind=wp) :: aLang, D1, D2, Dip_Eff, Dist3, DistI, dx, dy, dz, emx, ex, Fax, Fay, Faz, ffTots, FTest, ftot, FTots, fx, fy, fz, &
                 ghx, ghx1, ghy, ghy1, ghz, ghz1, QQO, r2, r2I, rx, ry, rz, S, Scal, ska, Temp, Tr1, TR2, uInd, V, v_Dummy, x
logical(kind=iwp) :: lExcl, NonEq
#ifdef _DEBUGPRINT_
real(kind=wp) :: DDotR, Del, DipAbs, RadAbs, TestA
#endif

#ifdef _DEBUGPRINT_
call RecPrt('edip: dEF(permanent) ',' ',dEF,4,nGrid_Eff)
call RecPrt('edip: PolEff ',' ',PolEff,nPolComp,nGrid_Eff)
call RecPrt('edip: DipEff ',' ',DipEff,1,nGrid_Eff)
call RecPrt('edip: Grid ',' ',Grid,3,nGrid_Eff)
write(u6,*) 'nGrid_Eff,nPolComp,nAnisopol,tk,dampIter,dipCutoff,clim,lDamping',nGrid_Eff,nPolComp,nAnisopol,tk,dampIter,dipCutoff, &
            clim,lDamping
do i=1,nGrid_Eff
  write(u6,*) 'EDOTr ',i,Grid(1,i)*dEF(1,i)+Grid(2,i)*dEF(2,i)+Grid(3,i)*dEF(3,i)
end do
#endif

qqo = Zero

NonEq = .false.

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'Iter fmax             testa'
#endif
Iter = 0
do
# ifdef _DEBUGPRINT_
  testa = fmax*afac
# endif
  Iter = Iter+1
  ! Loop over Langevin grid and make EF and dipol moments at the
  ! grid self consistent.

  do iGrid=1,nGrid_Eff

    fx = dEF(1,iGrid)+EF(1,iGrid)
    fy = dEF(2,iGrid)+EF(2,iGrid)
    fz = dEF(3,iGrid)+EF(3,iGrid)
    ftot = fx*fx+fy*fy+fz*fz
    ! Update EF and square norm

    EF(1,iGrid) = fx
    EF(2,iGrid) = fy
    EF(3,iGrid) = fz
    EF(4,iGrid) = ftot

    ! Reset update vector

    dEF(:,iGrid) = Zero

  end do
  do iGrid=1,nGrid_Eff
    fx = EF(1,iGrid)
    fy = EF(2,iGrid)
    fz = EF(3,iGrid)
    ftot = EF(4,iGrid)

    ! Skip if square norm below threshold

    !if (dEF(4,iGrid) < testa) cycle

    ghx = Grid(1,iGrid)
    ghy = Grid(2,iGrid)
    ghz = Grid(3,iGrid)

    ! Pick up dipole moment at grid point

    dx = DipMom(1,iGrid)
    dy = DipMom(2,iGrid)
    dz = DipMom(3,iGrid)

    !Dip_Eff = DipEff(iGrid)*real(Min(Iter,100),kind=wp)/100.0_wp
    Dip_Eff = DipEff(iGrid)

    ! Compute new dipole moment as a function of the EF, effective dipole
    ! moment and effective polarizability.

    if (Dip_Eff < 1.0e-10_wp) then
      if (iGrid > nAnisoPol) then   ! isotropic
        DipMom(1,iGrid) = fx*PolEff(1,iGrid)
        DipMom(2,iGrid) = fy*PolEff(1,iGrid)
        DipMom(3,iGrid) = fz*PolEff(1,iGrid)
      else  ! anisotropic
        DipMom(1,iGrid) = fx*PolEff(1,iGrid)+fy*PolEff(2,iGrid)+fz*PolEff(3,iGrid)
        DipMom(2,iGrid) = fx*PolEff(2,iGrid)+fy*PolEff(4,iGrid)+fz*PolEff(5,iGrid)
        DipMom(3,iGrid) = fx*PolEff(3,iGrid)+fy*PolEff(5,iGrid)+fz*PolEff(6,iGrid)
      end if
    else   ! NB!! Only isotropic implemented for dip>0
      fftots = sqrt(ftot)
      ftots = One/fftots
      x = Dip_Eff*tk*fftots
      ex = exp(x)
      emx = One/ex
      alang = (ex+emx)/(ex-emx)-One/x
      !alang = x/Three  ! Linear approximation
      i = iGrid
      uind = Dip_Eff*alang+ftot*PolEff(1,iGrid)*ftots
      DipMom(1,iGrid) = uind*fx*ftots
      DipMom(2,iGrid) = uind*fy*ftots
      DipMom(3,iGrid) = uind*fz*ftots
    end if

    ! Grid

    ! Compute the change in the dipole moment between the old (dx,dy,dz) and
    ! the new (DipMom).

    ! Try damping the change in dipole moment for better convergence

    DipMom(1,iGrid) = (One-dampIter)*DipMom(1,iGrid)+dampIter*dx
    DipMom(2,iGrid) = (One-dampIter)*DipMom(2,iGrid)+dampIter*dy
    DipMom(3,iGrid) = (One-dampIter)*DipMom(3,iGrid)+dampIter*dz

    dx = DipMom(1,iGrid)-dx
    dy = DipMom(2,iGrid)-dy
    dz = DipMom(3,iGrid)-dz

    ! Given the charge (qqo=0.0) and the change of the dipole moment
    ! at this point modify the multipole expansion around the origin
    ! accordingly. On the first iteration we will have the MM of the
    ! QM in Cavxyz too, in subsequential iterations we will only deal
    ! with incremental contributions.

    call qlm(ghx,ghy,ghz,qqo,dx,dy,dz,lMax,Cavxyz)

    ! Loop over the whole grid and update the EF due to the change of
    ! the dipole moment at the grid point "iGrid".

    Tr1 = Zero
    if (lDamping) then
      if (iGrid > nAnisopol) then
        Tr1 = PolEff(1,iGrid)
      else
        Tr1 = (PolEff(1,iGrid)+PolEff(4,iGrid)+PolEff(6,iGrid))/Three
      end if
    end if
    do jGrid=1,nGrid_Eff
      if (iGrid == jGrid) cycle
      scal = One
      if (lAmberpol .and. (iXPolType > 0) .and. (iGrid <= nXF) .and. (jGrid <= nXF)) then
        lExcl = .false.
        do i=1,nXMolnr
          if (XMolnr(1,jGrid) == XMolnr(i,iGrid)) lExcl = .true.
          if (XMolnr(1,jGrid) == -XMolnr(i,iGrid)) scal = scal14
        end do
        ! exclude field from iGrid when calculating the field at jGrid
        if (lExcl) cycle
        !if (lExcl) then
        !  write(u6,*) 'EXCLUDE dip',iGrid,' at ',jGrid
        !  cycle
        !else if (scal < One) then
        !  write(u6,*) 'SCALE dip',iGrid,' at ',jGrid,' with ',scal
        !end if
      end if
      ghx1 = Grid(1,jGrid)
      ghy1 = Grid(2,jGrid)
      ghz1 = Grid(3,jGrid)
      rx = ghx-ghx1
      ry = ghy-ghy1
      rz = ghz-ghz1
      r2 = (rx*rx+ry*ry+rz*rz)
      if (r2 < dipCutoff**2) cycle

      r2i = One/r2
      ska = dx*rx+dy*ry+dz*rz
      disti = sqrt(r2i)
      dist3 = r2i*disti
      temp = Three*ska*r2i
      if (lDamping) then
        if (jGrid > nAnisopol) then
          Tr2 = PolEff(1,jGrid)
        else
          Tr2 = (PolEff(1,jGrid)+PolEff(4,jGrid)+PolEff(6,jGrid))/Three
        end if
        !FIXME: what is this number?
        s = 2.3268_wp*(Tr1*Tr2)**(One/Six)
        v = min(One,sqrt(r2)/s)
        d1 = Four*v**3-Three*v**4
        d2 = v**4
        !write(u6,*) 'DAMP',d1,d2,Tr1,Tr2,sqrt(r2)
        dEF(1,jGrid) = dEF(1,jGrid)-(dx*d1-temp*rx*d2)*dist3*scal
        dEF(2,jGrid) = dEF(2,jGrid)-(dy*d1-temp*ry*d2)*dist3*scal
        dEF(3,jGrid) = dEF(3,jGrid)-(dz*d1-temp*rz*d2)*dist3*scal
      else
        dEF(1,jGrid) = dEF(1,jGrid)-(dx-temp*rx)*dist3*scal
        dEF(2,jGrid) = dEF(2,jGrid)-(dy-temp*ry)*dist3*scal
        dEF(3,jGrid) = dEF(3,jGrid)-(dz-temp*rz)*dist3*scal
      end if
    end do         ! jGrid

  end do           ! iGrid

  if (lRFCav) then

    ! Compute the charge distribution on the boundary of the cavity due to the
    ! MM expansion at origin.

    Cavxyz(:) = Ravxyz(:)

    call AppFld(Ravxyz,rds,Eps,lMax,EpsInf,NonEq)

    ! Compute EF at the grid due to the charge distribution in MM expansion
    ! for the QM system plus the dipole moments on the grid.

    do iGrid=1,nGrid_Eff
      ghx1 = Grid(1,iGrid)
      ghy1 = Grid(2,iGrid)
      ghz1 = Grid(3,iGrid)
      fax = Zero
      fay = Zero
      faz = Zero

      ! Given the charge distribution on the boundary of the cavity
      ! compute EF at (ghx1,ghy1,ghz1).

      call hmod(ghx1,ghy1,ghz1,v_dummy,fax,fay,faz,Ravxyz,lmax)

      ! Accumulate in update vector

      dEF(1,iGrid) = dEF(1,iGrid)+fax
      dEF(2,iGrid) = dEF(2,iGrid)+fay
      dEF(3,iGrid) = dEF(3,iGrid)+faz
    end do   ! iGrid

  end if ! if (lRFCav)

  fmax = Zero
  do iGrid=1,nGrid_Eff
    ftest = dEF(1,iGrid)**2+dEF(2,iGrid)**2+dEF(3,iGrid)**2
    dEF(4,iGrid) = ftest
    fmax = max(ftest,fmax)
  end do     ! iGrid

  Cavxyz(:) = Zero

  ! Check convergence

  !call RecPrt('DipMom ',' ',DipMom,3,nGrid_Eff)

# ifdef _DEBUGPRINT_
  write(u6,*) Iter,fmax,testa
# endif
  if (fmax <= clim) exit
end do

! Now we have a MM from QM + Langevin grid which is consistent with the
! charge distribution on the boundary of the cavity. The Langevin
! distribution of dipole moments is also internally consistent!

#ifdef _DEBUGPRINT_
call RecPrt('edip: converged DipMom ',' ',DipMom,3,nGrid_Eff)

! Write out dipoles and a pointcharge representation of the dipoles
write(u6,*) 'QREP'
do i=1,nGrid_Eff
  dipabs = sqrt(DipMom(1,i)**2+DipMom(2,i)**2+DipMom(3,i)**2)
  del = 0.01_wp
  write(u6,*) Grid(1,i)+DipMom(1,i)/dipabs*del,Grid(2,i)+DipMom(2,i)/dipabs*del,Grid(3,i)+DipMom(3,i)/dipabs*del,dipabs/del*Half
  write(u6,*) Grid(1,i)-DipMom(1,i)/dipabs*del,Grid(2,i)-DipMom(2,i)/dipabs*del,Grid(3,i)-DipMom(3,i)/dipabs*del,-dipabs/del*Half
end do

do i=1,nGrid_Eff
  ddotr = Grid(1,i)*DipMom(1,i)+Grid(2,i)*DipMom(2,i)+Grid(3,i)*DipMom(3,i)
  dipabs = sqrt(DipMom(1,i)*DipMom(1,i)+DipMom(2,i)*DipMom(2,i)+DipMom(3,i)*DipMom(3,i))
  radabs = sqrt(Grid(1,i)*Grid(1,i)+Grid(2,i)*Grid(2,i)+Grid(3,i)*Grid(3,i))
  write(u6,*) 'RADPOL',radabs,dipabs/(scala*scalb*scalc),ddotr/(dipabs*radabs)
end do
#endif

return

end subroutine edip
