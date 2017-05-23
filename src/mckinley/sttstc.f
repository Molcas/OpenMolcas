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
      SubRoutine Sttstc
      Implicit Real*8(a-h,o-z)
#include "cputime.fh"
      Character*50 NamFld(nTotal+1)
      Character*60 Fmt
      Data NamFld
     &    / '1) Calculation of one electron integrals        :',
     &      '2) Calculation of two electron integrals        :',
     &      '     a) Decontraction of two electron density   :',
     &      '     b) Integral evalution & 2nd derivatives    :',
     &      '     c) Screening                               :',
     &      '     d) Transfromation of integrals             :',
     &      '     e) Direct Fock matrix generation           :',
     &      '     f) Direct MO transformation                :',
     &      '3)  Control and input                           :',
     &      '   T O T A L                                    :' /
*---- Write out timing informations
      Fmt='(2x,A)'
      Write(6,*)
      Call CollapseOutput(1,'Statistics and timing')
      Write(6,'(3X,A)')     '---------------------'
      Write(6,*)
      Write(6,Fmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
     &          //' - - - - - - - - -'
      Write(6,Fmt)'   Part of the program                           '
     &          //'   CPU    fraction'
      Write(6,Fmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
     &          //' - - - - - - - - -'
      If (CPUStat(nTotal).gt.0.01D0) Then
         TotCpu=CPUStat(nTotal)
      Else
         TotCpu=0.01D0
      End If
*     TotCpu=max(0.01,CPUStat(nTotal))
      CPUStat(nTwoel)=CPUStat(nIntegrals)+CPUStat(nScreen)+
     &                CPUStat(nTrans)+CPUStat(nTwoDens)+
     &                CPUStat(nFckAck)+CPUStat(nMOTrans)
      Diverse=CPUStat(nTotal)-CPUStat(nTwoEl)-CPUStat(nOneel)

      Do iFld = 1, nTotal - 1
         Write(6,'(2x,A45,2f10.2)')NamFld(iFld),CPUStat(iFld),
     &                             CPUStat(iFld)/TotCpu
      End Do
      Write(6,'(2x,A45,2f10.2)')NamFld(nTotal),diverse,
     &                          diverse/TotCpu

      Write(6,*)
      Write(6,'(2x,A45,2F10.2)')NamFld(nTotal+1),TotCpu
      Write(6,Fmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
     &          //' - - - - - - - - -'
      Call CollapseOutput(0,'Statistics and timing')
      Write(6,*)
      Return
      End
