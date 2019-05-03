************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine ReadIn_Polar(NoField,Delta,MpProp_Level,Bond_Threshold,
     &                        iPlot,iPrint,Standard,Opt_Method,UserDen,
     &                        PrintDen,SubtractDen,SubScale,Restart,
     &                        TDensity,nStateI,nStateF,XHole,Diffuse,
     &                        dLimmo,Thrs1,Thrs2,nThrs,ThrsMul,Alpha,
     &                        LIonize)
      Implicit Real*8 (a-h,o-z)
      Logical NoField, Standard, UserDen, PrintDen, SubtractDen
      Logical TDensity, XHole, Diffuse(3)
      Logical LIonize,Restart
      Character*12  Opt_Method
      Dimension dLimmo(2)
*
      Call qEnter('ReadIn')
*
*     copy input from standard input to a local scratch file
*
      LuSpool=21
      Call SpoolInp(LuSpool)
*
*     read input
*
      Call RdInp_Polar(LuSpool,NoField,Delta,MpProp_Level,
     &          Bond_Threshold,iPlot,iPrint,Standard,Opt_Method,UserDen,
     &          PrintDen,SubtractDen,SubScale,Restart,TDensity,nStateI,
     &          nStateF,XHole,Diffuse,dLimmo,Thrs1,Thrs2,nThrs,ThrsMul,
     &          Alpha, LIonize)
*
*     remove local copy of standard input
*
      Call Close_LuSpool(LuSpool)
      Call qExit('ReadIn')
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
