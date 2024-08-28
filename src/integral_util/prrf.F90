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
! Copyright (C) 1992,2000, Roland Lindh                                *
!***********************************************************************

subroutine PrRF(DSCF,NonEq,iCharge,jPrint)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!                                                                      *
!             Modified for Langevin polarizabilities, Marsk 2000 (RL)  *
!***********************************************************************

use External_Centers, only: nXF, iXPolType
use rctfld_module, only: lRF, PCM, lRFCav, Eps, EpsInf, rds, lMax, lLangevin, latato, RadLat, ScalA, ScalB, ScalC, ScaAA, PolSI, &
                         DipSI, gAtom, DieDel, TK, cLim, aFac, nExpo, PreFac, Solvent, Conductor, CordSI, RslPar
use Definitions, only: u6

implicit none
logical DSCF, NonEq
integer iCharge, jPrint
integer StrnLn, i, j, nSolvent
real*8 AArea, R_Min_Sphere

if (jPrint >= 2) then
  if (lRF .and. (.not. PCM) .and. lRFCav) then
    write(u6,*)
    write(u6,'(5X,A)') 'Reaction Field calculation: the Kirkwood model'
    write(u6,'(5X,A,ES10.3)') ' Dielectric Constant :',Eps
    write(u6,'(5X,A,ES10.3)') ' Eps_opt             :',EpsInf
    write(u6,'(5X,A,ES10.3)') ' Radius of Cavity(au):',rds
    write(u6,'(5X,A,I2)') ' l_Max               :',lMax
    if (NonEq) then
      write(u6,'(5X,A)') ' Calculation type    : non-equilibrium'
    else
      write(u6,'(5X,A)') ' Calculation type    : equilibrium'
    end if
    write(u6,*)
  end if
  if (iXPolType > 0) then
    write(u6,*)
    write(u6,'(5X,A)') ' Explicit polarisabilities activated'
    write(u6,'(5X,A)') ' -----------------------------------'
    write(u6,'(5X,A,I2)') ' Number of points    :',nXF
    if (iXPolType == 1) then
      write(u6,'(5X,A)') ' Polarisabilities are isotropic'
    else if (iXPolType == 2) then
      write(u6,'(5X,A)') ' Polarisabilities are anisotropic'
    end if
    write(u6,*)
  end if
  if (lLangevin) then
    write(u6,'(5X,A)') 'Langevin dipole moments activated'
    write(u6,'(5X,A,I2)') ' Gitter type         :',latato
    write(u6,'(5X,A)') ' Gitter centers'
    do i=1,latato
      write(u6,'(3(5x,F10.4))') (Cordsi(j,i),j=1,3)
    end do
    write(u6,'(5X,A,ES10.3)') ' Max. Latt. Extn(au) :',radlat
    write(u6,'(5X,A,F10.4)') ' Cell dimensions     :',scala
    write(u6,'(5X,A,F10.4)') '                      ',scalb
    write(u6,'(5X,A,F10.4)') '                      ',scalc
    write(u6,'(5X,A,F10.4)') ' Overal scaling      :',scaaa
    write(u6,'(5X,A,F10.4)') ' Site polarizability :',polsi
    write(u6,'(5X,A,F10.4)') ' Site dipole moment  :',dipsi
    write(u6,'(5X,A,F10.4)') ' Atoms in the latt.  :',gAtom
    write(u6,'(5X,A,F10.4)') ' Diel. delete param. :',diedel
    write(u6,'(5X,A,F10.4)') ' Inverse Boltzman f. :',tK
    write(u6,'(5X,A,ES10.1)') ' clim                :',clim
    write(u6,'(5X,A,F10.4)') ' afac                :',afac
    write(u6,'(5X,A,I10  )') ' nexp                :',nexpo
    write(u6,'(5X,A,F10.4)') ' prefac              :',prefac
    write(u6,*)
  end if
  if (PCM) then
    aArea = RSlPar(7)
    r_min_Sphere = RSlPar(3)
    write(u6,*)
    write(u6,'(5X,A)') ' Polarizable Continuum Model (PCM) activated'
    nSolvent = StrnLn(Solvent)
    write(u6,'(5X,A,A)') ' Solvent: ',Solvent(1:nSolvent)
    if (.not. Conductor) then
      write(u6,'(5X,A)') ' Version: Dielectric'
    else
      write(u6,'(5X,A)') ' Version: Conductor'
    end if
    write(u6,'(5X,A,F6.4,A)') ' Average area for surface element on the cavity boundary: ',aArea,' angstrom^2'
    write(u6,'(5X,A,F6.4,A)') ' Minimum radius for added spheres: ',r_min_Sphere,' angstrom'
    if (NonEq) then
      write(u6,'(5X,A)') ' Calculation type: non-equilibrium (slow component from JobOld)'
    else
      write(u6,'(5X,A)') ' Calculation type: equilibrium'
    end if
    write(u6,*)
  end if
end if

if (lRF) call Init_RctFld(NonEq,iCharge)
if (DSCF) call Allok2()

return

end subroutine PrRF
