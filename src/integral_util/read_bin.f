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
      Subroutine Read_Bin(iShell_A,iShell_B,iShell_C,iShell_D,G_Toc,
     &                    nQuad,Gamma,nGamma,LuGamma,Bin,lBin)
      Implicit Real*8 (a-h,o-z)
#include "SysDef.fh"
#include "real.fh"
      Real*8 G_Toc(nQuad), Bin(2,lBin), Gamma(nGamma)
*
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*
      Call FZero(Gamma,nGamma)
*
      iShell_AB=iTri(iShell_A,iShell_B)
      iShell_CD=iTri(iShell_C,iShell_D)
      iShell_ABCD=iTri(iShell_AB,iShell_CD)
*
*     Write (*,*) 'Reading Gammas for shell quadruplet ', iShell_ABCD
*
      iDisk=Int(G_Toc(iShell_ABCD))
      iRead=2
99    Call dDaFile(LuGamma,iRead,Bin,2*lBin,iDisk)
      lGamma=Int(Bin(1,lBin))
      iDisk =Int(Bin(2,lBin))
*
      Do iGamma = 1, lGamma
         jGamma=Int(Bin(2,iGamma))
         If (jGamma.gt.nGamma) Then
            Call WarningMessage(2,'Read_Bin: jGamma.gt.nGamma')
            Call Abend()
         End If
         Gamma(jGamma)=Bin(1,iGamma)
*        Write (*,*) Gamma(jGamma), jGamma
      End Do
*
      If (iDisk.ge.0) Go To 99
*
      Return
      End

      Subroutine Read_Bin_Columbus
     &              (iShell_A,iShell_B,iShell_C,iShell_D,G_Toc,
     &                    nQuad,Gamma,nGamma,LuGamma,Bin,lBin)
      Implicit Real*8 (a-h,o-z)
#include "SysDef.fh"
#include "real.fh"
      Real*8 G_Toc(nQuad), Bin(2,lBin), Gamma(nGamma)
      Integer idisk_save
*
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*
      Call FZero(Gamma,nGamma)
*
      iShell_AB=iTri(iShell_A,iShell_B)
      iShell_CD=iTri(iShell_C,iShell_D)
      iShell_ABCD=iTri(iShell_AB,iShell_CD)
*
*     Write (*,*) 'Reading Gammas for shell quadruplet ', iShell_ABCD
*
      iDisk=Int(G_Toc(iShell_ABCD))
*     write(6,'(a,i6,a,i6)') ' Reading gtoc for ',ishell_abcd,'=',iDisk
      iRead=2
* ddafile overwrites idisk
99    idisk_save=idisk
      Call dDaFile(LuGamma,iRead,Bin,2,iDisk)
      lGamma=Int(Bin(1,1))
*98    format(a,'got header1=',i10,i10)
*     Write (*,98) '1st:',lgamma,Int(Bin(2,1))
      If (lGamma.gt.lBin) Then
            Call WarningMessage(2,'Read_Bin_Columbus: lGamma.gt.lbin')
            Call Abend()
      End If
      iDisk=idisk_save
*     write(*,'(a,i6)') ' second read from idisk=',iDisk
      Call dDaFile(LuGamma,iRead,Bin,2*lGamma,iDisk)
      iDisk =Int(Bin(2,1))
*     Write (*,98) '2nd:',Int(Bin(1,1)),Int(Bin(2,1))
*
      Do iGamma = 2, lGamma
         jGamma=Int(Bin(2,iGamma))
         If (jGamma.gt.nGamma) Then
            Call WarningMessage(2,'Read_Bin_Columbus: jGamma.gt.nGamma')
            Call Abend()
         End If
         Gamma(jGamma)=Bin(1,iGamma)
*        Write (*,*) Gamma(jGamma), jGamma
      End Do
*
      If (iDisk.ge.0) Go To 99
*
      Return
      End
