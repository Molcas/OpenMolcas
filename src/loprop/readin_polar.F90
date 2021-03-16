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

subroutine ReadIn_Polar(NoField,Delta,MpProp_Level,Bond_Threshold,iPlot,iPrint,Standard,Opt_Method,UserDen,PrintDen,SubtractDen, &
                        SubScale,Restart,TDensity,nStateI,nStateF,XHole,Diffuse,dLimmo,Thrs1,Thrs2,nThrs,ThrsMul,Alpha,LIonize)

implicit real*8(a-h,o-z)
logical NoField, Standard, UserDen, PrintDen, SubtractDen
logical TDensity, XHole, Diffuse(3)
logical LIonize, Restart
character*12 Opt_Method
dimension dLimmo(2)

! copy input from standard input to a local scratch file

LuSpool = 21
call SpoolInp(LuSpool)

! read input

call RdInp_Polar(LuSpool,NoField,Delta,MpProp_Level,Bond_Threshold,iPlot,iPrint,Standard,Opt_Method,UserDen,PrintDen,SubtractDen, &
                 SubScale,Restart,TDensity,nStateI,nStateF,XHole,Diffuse,dLimmo,Thrs1,Thrs2,nThrs,ThrsMul,Alpha,LIonize)

! remove local copy of standard input

call Close_LuSpool(LuSpool)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine ReadIn_Polar
