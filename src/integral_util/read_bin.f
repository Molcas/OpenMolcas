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
      Subroutine Read_Bin(iShell_A,iShell_B,iShell_C,iShell_D,G_Toc,
     &                    nQuad,Gamma,nGamma,LuGamma,Bin,lBin)
      use Constants, only: Zero
      Implicit None
#include "SysDef.fh"
      Integer iShell_A,iShell_B,iShell_C,iShell_D,
     &        nQuad,nGamma,LuGamma,lBin
      Real*8 G_Toc(nQuad), Bin(2,lBin), Gamma(nGamma)
!
      Integer iDisk, lGamma, iGamma, jGamma
      Integer iShell_AB, iShell_CD, iShell_ABCD
      Integer, Parameter:: iRead=2
      Integer i, j, iTri
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
!
      Gamma(:)=Zero
!
      iShell_AB=iTri(iShell_A,iShell_B)
      iShell_CD=iTri(iShell_C,iShell_D)
      iShell_ABCD=iTri(iShell_AB,iShell_CD)
!
!     Write (*,*) 'Reading Gammas for shell quadruplet ', iShell_ABCD
!
      iDisk=Int(G_Toc(iShell_ABCD))
      Do While (iDisk>=0)
         Call dDaFile(LuGamma,iRead,Bin,2*lBin,iDisk)
         lGamma=Int(Bin(1,lBin))
         iDisk =Int(Bin(2,lBin))
!
         Do iGamma = 1, lGamma
            jGamma=Int(Bin(2,iGamma))
            If (jGamma.gt.nGamma) Then
               Call WarningMessage(2,'Read_Bin: jGamma.gt.nGamma')
               Call Abend()
            End If
            Gamma(jGamma)=Bin(1,iGamma)
!           Write (*,*) Gamma(jGamma), jGamma
         End Do
!
      End Do
!
      End Subroutine Read_Bin

      Subroutine Read_Bin_Columbus
     &              (iShell_A,iShell_B,iShell_C,iShell_D,G_Toc,
     &                    nQuad,Gamma,nGamma,LuGamma,Bin,lBin)
      use Constants, only: Zero
      Implicit None
#include "SysDef.fh"
      Integer iShell_A,iShell_B,iShell_C,iShell_D,
     &        nQuad,nGamma,LuGamma,lBin
      Real*8 G_Toc(nQuad), Bin(2,lBin), Gamma(nGamma)

      Integer iShell_AB, iShell_CD, iShell_ABCD
      Integer, Parameter:: iRead=2
      Integer idisk_save
      Integer iDisk, lGamma, iGamma, jGamma
      Integer i, j, iTri
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
!
      Gamma(:)=Zero
!
      iShell_AB=iTri(iShell_A,iShell_B)
      iShell_CD=iTri(iShell_C,iShell_D)
      iShell_ABCD=iTri(iShell_AB,iShell_CD)
!
!     Write (*,*) 'Reading Gammas for shell quadruplet ', iShell_ABCD
!
      iDisk=Int(G_Toc(iShell_ABCD))
      Do While (iDisk>=0)
         idisk_save=idisk
         Call dDaFile(LuGamma,iRead,Bin,2,iDisk)
         lGamma=Int(Bin(1,1))
         If (lGamma.gt.lBin) Then
               Call WarningMessage(2,
     &                        'Read_Bin_Columbus: lGamma.gt.lbin')
               Call Abend()
         End If
         iDisk=idisk_save
         Call dDaFile(LuGamma,iRead,Bin,2*lGamma,iDisk)
         iDisk =Int(Bin(2,1))
!
         Do iGamma = 2, lGamma
            jGamma=Int(Bin(2,iGamma))
            If (jGamma.gt.nGamma) Then
               Call WarningMessage(2,
     &                     'Read_Bin_Columbus: jGamma.gt.nGamma')
               Call Abend()
            End If
            Gamma(jGamma)=Bin(1,iGamma)
!        Write (*,*) Gamma(jGamma), jGamma
         End Do
!
      End Do
!
      End Subroutine Read_Bin_Columbus
