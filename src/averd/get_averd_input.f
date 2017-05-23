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
*---------------------------------------------------------------------
* Read input for averd.
*
* Title                -Title
* Wset                -Set of weights. Dimension is Nset.
* iPrint        -How much print.(1=minimal,2=print average orbitals,
*                5=print all orbitals,99=wacko!)
* Nset                -Number of orbitals to read
* DensityBased        -Is the procedure density or orbital based
* ThrOcc        -Print which orbitals have occupation number below
*                this threshold.
*-----------------------------------------------------------------------
      Subroutine Get_Averd_input(Title,Wset,iPrint,Nset,DensityBased
     &                          ,ThrOcc)
      Implicit Real*8 (a-h,o-z)

#include "mxave.fh"
#include "warnings.fh"

      Character*72 Title
      Dimension Wset(MxSets)
      Logical DensityBased

      Character*180 Key
      Character*4 Kword
      Character*180 Get_Ln
      Integer iCLast
      External Get_Ln,iCLast

*
*-- Call subroutines that handle the input.
*
      LuRd=21
      Call SpoolInp(LuRd)
      Rewind(LuRd)
      Call RdNLst(LuRd,'AVERD')

*
*-- Label 1000 is the top.
*
1000  Continue

*
*-- Get_Ln read the keyword and skips line starting with *
*   or is empty.
*
      Key=Get_Ln(LuRd)
      Kword=Trim(Key)
      Call UpCase(Kword)

*
*-- The keywords...
*
      If (Kword(1:4).eq.'WSET') Go To 101
      If (Kword(1:4).eq.'PRIN') Go To 102
      If (Kword(1:4).eq.'TITL') Go To 103
      If (Kword(1:4).eq.'ORBI') Go To 104
      If (Kword(1:4).eq.'OCCU') Go To 105
      If (Kword(1:4).eq.'END ') Go To 9999

*
*-- ...and what happens if something else is encountered.
*
      iChrct=Len(KWord)
      Last=iCLast(KWord,iChrct)
      Write(6,*)' '
      Write(6,'(1X,A,A)')Kword(1:Last),' is not a valid keyword!'
      Write(6,*)' ERROR!'
      Call Quit(_RC_INPUT_ERROR_)

*
*-- Read weights.
*
101   Continue
      Key=Get_Ln(LuRd)
      Call Get_I(1,Nset,1)
      Key=Get_Ln(LuRd)
      Call Get_F(1,Wset,Nset)
      Go To 1000
*
*-- How much print?
*
102   Continue
      Key=Get_Ln(LuRd)
      Call Get_I(1,iPrint,1)
      Go To 1000
*
*-- Title
*
103   Continue
      Key=Get_Ln(LuRd)
      Title=Key(1:Len(Title))
      Go To 1000
*
*-- Should it be density based, or orbital based.
*
104   Continue
      DensityBased=.false.
      Go To 1000

*
*-- I want to be told which orbitals have occ.num. below threshold.
*
105   Continue
      Key=Get_Ln(LuRd)
      Call Get_F(1,ThrOcc,1)
      Go To 1000

*
*-- A most Graceful Exit.
*
9999  Continue

      Return
      End
