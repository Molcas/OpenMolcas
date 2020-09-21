************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1992,2000, Roland Lindh                                *
************************************************************************
      SubRoutine PrRF(DSCF,NonEq,iCharge,jPrint)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              Allok2                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*                                                                      *
*             Modified for Langevin polarizabilities, Marsk 2000 (RL)  *
************************************************************************
      use External_Centers, only: nXF, iXPolType
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
#include "rctfld.fh"
      Logical DSCF, NonEq
      Integer StrnLn
*
      Call qEnter('PrRF')
*
      IF (jPrint.GE.2) THEN
      If (lRF.and..Not.PCM.and.lRFCav) Then
         Write (6,*)
         Write (6,'(5X,A)')
     &       'Reaction Field calculation: the Kirkwood model'
         Write (6,'(5X,A,E10.3)') ' Dielectric Constant :',Eps
         Write (6,'(5X,A,E10.3)') ' Eps_opt             :',EpsInf
         Write (6,'(5X,A,E10.3)') ' Radius of Cavity(au):',rds
         Write (6,'(5X,A,I2)')    ' l_Max               :',lMax
         If (NonEq) Then
            Write (6,'(5X,A)')
     &               ' Calculation type    : non-equilibrium'
         Else
            Write (6,'(5X,A)')
     &               ' Calculation type    : equilibrium'
         End If
         Write (6,*)
      End If
      If (iXPolType.gt.0) Then
         Write (6,*)
         Write (6,'(5X,A)') ' Explicit polarisabilities activated'
         Write (6,'(5X,A)') ' -----------------------------------'
         Write (6,'(5X,A,I2)')     ' Number of points    :',nXF
         If (iXPolType.eq.1) Then
            Write (6,'(5X,A)')
     &           ' Polarisabilities are isotropic'
         ElseIf (iXPolType.eq.2) Then
            Write (6,'(5X,A)')
     &           ' Polarisabilities are anisotropic'
         EndIf
         Write (6,*)
      EndIf
      If (lLangevin) Then
         Write (6,'(5X,A)') 'Langevin dipole moments activated'
         Write (6,'(5X,A,I2)')     ' Gitter type         :',latato
         Write (6,'(5X,A)')        ' Gitter centers'
         Do i = 1, latato
            Write (6,'(3(5x,F10.4))') (Cordsi(j,i),j=1,3)
         End Do
         Write (6,'(5X,A,E10.3)')  ' Max. Latt. Extn(au) :',radlat
         Write (6,'(5X,A,F10.4)')  ' Cell dimensions     :',scala
         Write (6,'(5X,A,F10.4)')  '                      ',scalb
         Write (6,'(5X,A,F10.4)')  '                      ',scalc
         Write (6,'(5X,A,F10.4)')  ' Overal scaling      :',scaaa
         Write (6,'(5X,A,F10.4)')  ' Site polarizability :',polsi
         Write (6,'(5X,A,F10.4)')  ' Site dipole moment  :',dipsi
         Write (6,'(5X,A,F10.4)')  ' Atoms in the latt.  :',gAtom
         Write (6,'(5X,A,F10.4)')  ' Diel. delete param. :',diedel
         Write (6,'(5X,A,F10.4)')  ' Inverse Boltzman f. :',tK
         Write (6,'(5X,A,E10.1)')  ' clim                :',clim
         Write (6,'(5X,A,F10.4)')  ' afac                :',afac
         Write (6,'(5X,A,I10  )')  ' nexp                :',nexpo
         Write (6,'(5X,A,F10.4)')  ' prefac              :',prefac
         Write (6,*)
      End If
      If (PCM) Then
         aArea = RSlPar(7)
         r_min_Sphere = RSlPar(3)
         Write (6,*)
         Write (6,'(5X,A)')
     &         ' Polarizable Continuum Model (PCM) activated'
         nSolvent=StrnLn(Solvent)
         Write (6,'(5X,A,A)')
     &         ' Solvent: ',Solvent(1:nSolvent)
         If (.Not.Conductor) Then
            Write (6,'(5X,A)') ' Version: Dielectric'
         Else
            Write (6,'(5X,A)') ' Version: Conductor'
         End If
         Write (6,'(5X,A,F6.4,A)')
     &               ' Average area for surface element on the '
     &             //'cavity boundary: ', aArea, ' Angstrom^2'
         Write (6,'(5X,A,F6.4,A)')
     &               ' Minimum radius for added spheres: ',
     &               r_min_Sphere, ' Angstrom'
         If (NonEq) Then
            Write (6,'(5X,A)')
     &               ' Calculation type: non-equilibrium'
     &             //' (slow component from JobOld)'
         Else
            Write (6,'(5X,A)')
     &               ' Calculation type: equilibrium'
         End If
         Write (6,*)
      End If
      ENDIF

      If (lRF) Call Init_RctFld(NonEq,iCharge)
      If (DSCF) Call Allok2
*
      Call qExit('PrRF')
      Return
*
      End
