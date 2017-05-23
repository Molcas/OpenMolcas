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
      Subroutine Freeze_Default(iANr,nShell,lMax)
      Integer nShell(0:lMax)
      Parameter (nAtoms=106)
      Integer iDefaults(0:3,0:nAtoms)
*---- First row, H-He: no frozen orbitals
      Data ((iDefaults(i,j),i=0,3),j=  0,  2)/0,0,0,0,
     &                                        0,0,0,0,
     &                                        0,0,0,0/
*---- Second row, Li-Ne: 1s
      Data ((iDefaults(i,j),i=0,3),j=  3, 10)/1,0,0,0,
     &                                        1,0,0,0,
     &                                        1,0,0,0,
     &                                        1,0,0,0,
     &                                        1,0,0,0,
     &                                        1,0,0,0,
     &                                        1,0,0,0,
     &                                        1,0,0,0/
*---- Third row, Na-Si: 1s2s
      Data ((iDefaults(i,j),i=0,3),j= 11, 14)/2,0,0,0,
     &                                        2,0,0,0,
     &                                        2,0,0,0,
     &                                        2,0,0,0/
*---- Third row, P-Ar: 1s2s2p
      Data ((iDefaults(i,j),i=0,3),j= 15, 18) /2,1,0,0,
     &                                        2,1,0,0,
     &                                        2,1,0,0,
     &                                        2,1,0,0/
*---- Fourth row, K-Cu: 1s2s2p3s
      Data ((iDefaults(i,j),i=0,3),j= 19, 30)/3,1,0,0,
     &                                        3,1,0,0,
     &                                        3,1,0,0,
     &                                        3,1,0,0,
     &                                        3,1,0,0,
     &                                        3,1,0,0,
     &                                        3,1,0,0,
     &                                        3,1,0,0,
     &                                        3,1,0,0,
     &                                        3,1,0,0,
     &                                        3,1,0,0,
     &                                        3,1,0,0/
*---- Fourth row, Zn-Kr: 1s2s2p3s3p3d
      Data ((iDefaults(i,j),i=0,3),j= 31, 36)/3,2,1,0,
     &                                        3,2,1,0,
     &                                        3,2,1,0,
     &                                        3,2,1,0,
     &                                        3,2,1,0,
     &                                        3,2,1,0/
*---- Fifth row, Rb-Ag: 1s2s2p3s3p4s3d
      Data ((iDefaults(i,j),i=0,3),j= 37, 47)/4,2,1,0,
     &                                        4,2,1,0,
     &                                        4,2,1,0,
     &                                        4,2,1,0,
     &                                        4,2,1,0,
     &                                        4,2,1,0,
     &                                        4,2,1,0,
     &                                        4,2,1,0,
     &                                        4,2,1,0,
     &                                        4,2,1,0,
     &                                        4,2,1,0/
*---- Fifth row, Cd-Xe: 1s2s2p3s3p4s3d4p
      Data ((iDefaults(i,j),i=0,3),j= 48, 54)/4,3,1,0,
     &                                        4,3,1,0,
     &                                        4,3,1,0,
     &                                        4,3,1,0,
     &                                        4,3,1,0,
     &                                        4,3,1,0,
     &                                        4,3,1,0/
*---- Sixth row, Cs-Yb: 1s2s2p3s3p4s3d4p5s4d5p
      Data ((iDefaults(i,j),i=0,3),j= 55, 70)/4,3,2,0,
     &                                        4,3,2,0,
     &                                        4,3,2,0,
     &                                        4,3,2,0,
     &                                        4,3,2,0,
     &                                        4,3,2,0,
     &                                        4,3,2,0,
     &                                        4,3,2,0,
     &                                        4,3,2,0,
     &                                        4,3,2,0,
     &                                        4,3,2,0,
     &                                        4,3,2,0,
     &                                        4,3,2,0,
     &                                        4,3,2,0,
     &                                        4,3,2,0,
     &                                        4,3,2,0/
*---- Sixth row, Lu-Au: 1s2s2p3s3p4s3d4p4d4f
      Data ((iDefaults(i,j),i=0,3),j= 71, 79)/4,3,2,1,
     &                                        4,3,2,1,
     &                                        4,3,2,1,
     &                                        4,3,2,1,
     &                                        4,3,2,1,
     &                                        4,3,2,1,
     &                                        4,3,2,1,
     &                                        4,3,2,1,
     &                                        4,3,2,1/
*---- Sixth row, Hg-Rn: 1s2s2p3s3p4s3d4p4d4f5s5p
      Data ((iDefaults(i,j),i=0,3),j= 80, 86)/5,4,2,1,
     &                                        5,4,2,1,
     &                                        5,4,2,1,
     &                                        5,4,2,1,
     &                                        5,4,2,1,
     &                                        5,4,2,1,
     &                                        5,4,2,1/
*---- Seventh row, Fr-Cm: 1s2s2p3s3p4s3d4p5s4d5p4f5d
      Data ((iDefaults(i,j),i=0,3),j= 87, 96)/5,4,3,1,
     &                                        5,4,3,1,
     &                                        5,4,3,1,
     &                                        5,4,3,1,
     &                                        5,4,3,1,
     &                                        5,4,3,1,
     &                                        5,4,3,1,
     &                                        5,4,3,1,
     &                                        5,4,3,1,
     &                                        5,4,3,1/
      Data ((iDefaults(i,j),i=0,3),j= 97,103)/6,5,3,1,
     &                                        6,5,3,1,
     &                                        6,5,3,1,
     &                                        6,5,3,1,
     &                                        6,5,3,1,
     &                                        6,5,3,1,
     &                                        6,5,3,1/
      Data ((iDefaults(i,j),i=0,3),j=104,106)/6,5,3,2,
     &                                        6,5,3,2,
     &                                        6,5,3,2/
*
      If (iANr.lt.0.or.iANr.gt.nAtoms) Then
         Write (6,*) 'Freeze_Defaults: iAnr is out of range!'
         Write (6,*) 'iANr=',iANr
         Call Abend
      End If
*
      Call iCopy(lMax+1,0,0,nShell,1)
*
      Do i = 0, Min(lMax,3)
         nShell(i)=iDefaults(i,iAnr)
      End Do
*
      Return
      End
