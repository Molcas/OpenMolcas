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

implicit none
logical DSCF, NonEq
integer iCharge, jPrint
integer StrnLn, i, j, nSolvent
real*8 AArea, R_Min_Sphere

if (jPrint >= 2) then
  if (lRF .and. (.not. PCM) .and. lRFCav) then
    write(6,*)
    write(6,'(5X,A)') 'Reaction Field calculation: the Kirkwood model'
    write(6,'(5X,A,ES10.3)') ' Dielectric Constant :',Eps
    write(6,'(5X,A,ES10.3)') ' Eps_opt             :',EpsInf
    write(6,'(5X,A,ES10.3)') ' Radius of Cavity(au):',rds
    write(6,'(5X,A,I2)') ' l_Max               :',lMax
    if (NonEq) then
      write(6,'(5X,A)') ' Calculation type    : non-equilibrium'
    else
      write(6,'(5X,A)') ' Calculation type    : equilibrium'
    end if
    write(6,*)
  end if
  if (iXPolType > 0) then
    write(6,*)
    write(6,'(5X,A)') ' Explicit polarisabilities activated'
    write(6,'(5X,A)') ' -----------------------------------'
    write(6,'(5X,A,I2)') ' Number of points    :',nXF
    if (iXPolType == 1) then
      write(6,'(5X,A)') ' Polarisabilities are isotropic'
    else if (iXPolType == 2) then
      write(6,'(5X,A)') ' Polarisabilities are anisotropic'
    end if
    write(6,*)
  end if
  if (lLangevin) then
    write(6,'(5X,A)') 'Langevin dipole moments activated'
    write(6,'(5X,A,I2)') ' Gitter type         :',latato
    write(6,'(5X,A)') ' Gitter centers'
    do i=1,latato
      write(6,'(3(5x,F10.4))') (Cordsi(j,i),j=1,3)
    end do
    write(6,'(5X,A,ES10.3)') ' Max. Latt. Extn(au) :',radlat
    write(6,'(5X,A,F10.4)') ' Cell dimensions     :',scala
    write(6,'(5X,A,F10.4)') '                      ',scalb
    write(6,'(5X,A,F10.4)') '                      ',scalc
    write(6,'(5X,A,F10.4)') ' Overal scaling      :',scaaa
    write(6,'(5X,A,F10.4)') ' Site polarizability :',polsi
    write(6,'(5X,A,F10.4)') ' Site dipole moment  :',dipsi
    write(6,'(5X,A,F10.4)') ' Atoms in the latt.  :',gAtom
    write(6,'(5X,A,F10.4)') ' Diel. delete param. :',diedel
    write(6,'(5X,A,F10.4)') ' Inverse Boltzman f. :',tK
    write(6,'(5X,A,ES10.1)') ' clim                :',clim
    write(6,'(5X,A,F10.4)') ' afac                :',afac
    write(6,'(5X,A,I10  )') ' nexp                :',nexpo
    write(6,'(5X,A,F10.4)') ' prefac              :',prefac
    write(6,*)
  end if
  if (PCM) then
    aArea = RSlPar(7)
    r_min_Sphere = RSlPar(3)
    write(6,*)
    write(6,'(5X,A)') ' Polarizable Continuum Model (PCM) activated'
    nSolvent = StrnLn(Solvent)
    write(6,'(5X,A,A)') ' Solvent: ',Solvent(1:nSolvent)
    if (.not. Conductor) then
      write(6,'(5X,A)') ' Version: Dielectric'
    else
      write(6,'(5X,A)') ' Version: Conductor'
    end if
    write(6,'(5X,A,F6.4,A)') ' Average area for surface element on the cavity boundary: ',aArea,' angstrom^2'
    write(6,'(5X,A,F6.4,A)') ' Minimum radius for added spheres: ',r_min_Sphere,' angstrom'
    if (NonEq) then
      write(6,'(5X,A)') ' Calculation type: non-equilibrium (slow component from JobOld)'
    else
      write(6,'(5X,A)') ' Calculation type: equilibrium'
    end if
    write(6,*)
  end if
end if

if (lRF) call Init_RctFld(NonEq,iCharge)
if (DSCF) call Allok2()

return

end subroutine PrRF
