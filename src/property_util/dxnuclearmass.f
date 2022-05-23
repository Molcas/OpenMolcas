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
* Copyright (C) 2000, Per-Olof Widmark                                 *
*               2017, Ignacio Fdez. Galvan                             *
************************************************************************
      Real*8 Function dxNuclearMass(Z,A,Rc,Opt)
************************************************************************
*                                                                      *
* Routine: dxNuclearMass                                               *
* Purpose: To compute mass for isotope with Z protons and A-Z neutrons.*
*          Nontabulated isotopes are computed by the semiempirical     *
*          mass formula.                                               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University, Sweden.                                    *
* Written: March 2000                                                  *
* History: March 2017, use isotopes module, Ignacio Fdez. Galvan       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Parameters:                                                          *
* Z   - The nuclear charge for the isotope. Input.                     *
* A   - The number of nucleons for the isotope. Input.                 *
* Rc  - Return code. Output                                            *
* Opt - Options. Input.                                                *
*                                                                      *
************************************************************************
      Use Isotopes
      Implicit None
#include "proputil.fh"
#include "constants2.fh"
*----------------------------------------------------------------------*
* Parameters.                                                          *
*----------------------------------------------------------------------*
      Integer    StopOnWarning
      Parameter (StopOnWarning=_OPT_STOP_ON_WARNING_)
      Integer    StopOnError
      Parameter (StopOnError=_OPT_STOP_ON_ERROR_)
*----------------------------------------------------------------------*
* Dummy parameters.                                                    *
*----------------------------------------------------------------------*
      Integer Z
      Integer A
      Integer Rc
      Integer Opt
*----------------------------------------------------------------------*
* Local variables.                                                     *
*----------------------------------------------------------------------*
      Real*8  Coef(7)
      Data    Coef/1.00781360D0,1.00866184D0,0.01685183D0,0.01928950D0,
     &             0.00075636D0,0.10146129D0,0.02449108D0/
      Save    Coef
      Real*8  t
*----------------------------------------------------------------------*
* Search table.                                                        *
*----------------------------------------------------------------------*
      dxNuclearMass=NuclideMass(Z,A)
*----------------------------------------------------------------------*
* Optionally use the semi-empirical mass formula.                      *
*----------------------------------------------------------------------*
      If(dxNuclearMass.lt.0.0d0) Then
         Write(6,'(a)') '***'
         Write(6,'(a)') '*** dxNuclearMass: warning'
         Write(6,'(2a,2i6)') '*** semi empirical mass formula used for',
     &              ' nuclei (Z,A)=',Z,A
         Write(6,'(a)') '***'
         If(iAnd(StopOnWarning,Opt).ne.0) Call Quit_OnUserError
         t=0.0d0
         t=t+Coef(1)*Z
         t=t+Coef(2)*(A-Z)
         t=t-Coef(3)*A
         t=t+Coef(4)*A**(2.0d0/3.0d0)
         t=t+Coef(5)*Z*(Z-1)/Dble(A)**(1.0d0/3.0d0)
         t=t+Coef(6)*(Z-0.5d0*A)**2/Dble(A)
         If( mod(Z,2).eq.0 .and. mod(A,2).eq.0 ) Then
            t=t-Coef(7)/Dble(A)**(0.75)
         End If
         If( mod(Z+1,2).eq.0 .and. mod(A,2).eq.0 ) Then
            t=t+Coef(7)/Dble(A)**(0.75)
         End If
         dxNuclearMass=uToau*t
      End If
*----------------------------------------------------------------------*
* Done.                                                                *
*----------------------------------------------------------------------*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(Rc)
      End
