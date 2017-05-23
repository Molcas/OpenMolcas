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
      Real*8 Function CovRad(i)
      Implicit Real*8 (a-h,o-z)
      Real*8 CovRad_(0:86)
*     Character*144 Warning
      Data CovRad_/
     &  0.0d0,
     &  0.643d0,0.643d0,2.457d0,1.909d0,1.587d0,1.436d0,1.209d0, !1-7
     &  1.096d0,1.020d0,0.945d0,2.986d0,2.646d0,2.400d0,2.192d0, !8-14
     &  2.060d0,1.890d0,1.795d0,1.701d0,3.836d0,3.288d0,2.721d0, !15-21
     &  2.494d0,2.305d0,2.230d0,2.211d0,2.211d0,2.192d0,2.173d0, !22-28
     &  2.211d0,2.362d0,2.381d0,2.305d0,2.268d0,2.192d0,2.154d0, !29-35
     &  2.116d0,4.082d0,3.609d0,3.061d0,2.740d0,2.532d0,2.457d0, !36-42
     &  2.400d0,2.362d0,2.362d0,2.419d0,2.532d0,2.797d0,2.721d0, !43-49
     &  2.665d0,2.646d0,2.570d0,2.513d0,2.476d0,4.441d0,3.742d0, !50-56
     &  3.194d0,3.118d0,3.118d0,3.099d0,3.080d0,3.061d0,3.496d0, !57-63
     &  3.042d0,3.005d0,3.005d0,2.986d0,2.967d0,2.948d0,2.948d0, !64-70
     &  2.948d0,2.721d0,2.532d0,2.457d0,2.419d0,2.381d0,2.400d0, !71-77
     &  2.457d0,2.532d0,2.816d0,2.797d0,2.778d0,2.759d0,2.759d0, !78-84
     &  2.740d0,2.710d0/                                         !85-
*
      If (i.gt.86) Then
*        Write (Warning,'(2A)') 'CovRad: Warning i.gt.86!,;'//
*    &               'Guestimate of 2.70 au is used!'
*        Call WarningMessage(1,Warning)
         CovRad=2.70D0
      Else
         CovRad=CovRad_(i)
      End If
*
      Return
      End
