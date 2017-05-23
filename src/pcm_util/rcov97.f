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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Function RCov97(IA,IB)
      Implicit Real*8(A-H,O-Z)
C
C     This function returns an estimated covalent bond distance (in Ang)
C     between atoms of atomic numbers IA and IB. Setting IB to 0 returns
C     the covalent radius of IA. Parameters for atoms heavier than At are
C     taken from the UFF force field (A.K.Rappe',C.J.Casewit,K.S.Colwell,
C     W.A.Goddard III and W.M.Skiff, J.Am.Chem.Soc. 114,10024 (1992)).
C
      Parameter (MxAtN=104)
      Dimension Rii(0:MxAtN)
      Save Rii
      Data Rii/0.0d0,
     1  0.354D0,0.849D0,
     2  1.336D0,1.010D0,0.838D0,0.757D0,0.700D0,0.658D0,0.668D0,0.920D0,
     3  1.539D0,1.421D0,1.244D0,1.117D0,1.101D0,1.064D0,1.044D0,1.032D0,
     4  1.953D0,1.761D0,
     D  1.513D0,1.412D0,1.402D0,1.345D0,1.382D0,
     D  1.270D0,1.241D0,1.164D0,1.302D0,1.193D0,
     4                  1.260D0,1.197D0,1.211D0,1.190D0,1.192D0,1.147D0,
     5  2.260D0,2.052D0,
     D  1.698D0,1.564D0,1.473D0,1.467D0,1.322D0,
     D  1.478D0,1.332D0,1.338D0,1.386D0,1.403D0,
     5                  1.459D0,1.398D0,1.407D0,1.386D0,1.382D0,1.267D0,
     6  2.570D0,2.277D0,
     F  1.943D0,1.841D0,1.823D0,1.816D0,1.801D0,1.780D0,1.771D0,
     F  1.735D0,1.732D0,1.710D0,1.696D0,1.673D0,1.660D0,1.637D0,
     D  1.671D0,1.611D0,1.511D0,1.392D0,1.372D0,
     D  1.372D0,1.371D0,1.364D0,1.262D0,1.340D0,
     6                  1.518D0,1.459D0,1.512D0,1.500D0,1.545D0,1.420D0,
     7  2.880D0,2.512D0,1.983D0,1.721D0,1.711D0,1.684D0,1.666D0,
     7  1.657D0,1.660D0,1.801D0,1.761D0,1.750D0,1.724D0,1.712D0,
     7  1.689D0,1.679D0,1.698D0,
     $  1.850D0/
      RCov97 = Rii(Min(Max(IA,0),MxAtN)) + Rii(Min(Max(IB,0),MxAtN))
      Return
      End
