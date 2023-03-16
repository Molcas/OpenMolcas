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
! Info for LDF SCF.
Module LDFSCF
Logical LDF_Timing
Logical LDF_UseConventionalIntegrals, LDF_UseLDFIntegrals
Logical LDF_UsePSDIntegrals, LDF_UseExactIntegralDiagonal
Logical LDF_IntegralCheck, LDF_FitVerification
Logical LDF_CoefficientCheck, LDF_UBCNorm
Logical LDF_CoulombCheck, LDF_OverlapCheck
Logical LDF_ChargeCheck, LDF_ChargePrint
Logical LDF_ModeCheck
Integer LDF_IntegralMode, LDF_IntegralPSDCheck
Integer LDF_UseExactIntegrals
Real*8  LDFracMem
Real*8  LDF_IntegralPrescreening, LDF_ContributionPrescreening
End Module LDFSCF
