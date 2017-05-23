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
      Subroutine MandatoryInp(YesNo)
      Implicit Real*8 (a-h,o-z)

#include "warnings.fh"

      Logical YesNo(20)

      If(.not.YesNo(8).and..not.YesNo(7)) then
        Write(6,*)
        Write(6,*)' You have not specified what type of calculation'
     &//' this is.'
        Write(6,*)' Use either the RUN keyword or the SINGle-point'
     &//' keyword.'
        Call Quit(_RC_INPUT_ERROR_)
      Endif

      If(YesNo(3).and.YesNo(4)) then
        Write(6,*)
        Write(6,*)' You have specified both a SCFSection and a'
     &//' RASSisection.'
        Write(6,*)' They are mutually exclusive. Remove one.'
        Call Quit(_RC_INPUT_ERROR_)
      Endif

      If(YesNo(7).and..not.YesNo(6)) then
        Write(6,*)
        Write(6,*)' You have requested a single-point calculation, but'
     &//' no input coordinates were given.'
        Write(6,*)' Provide these in the SOLVent section.'
        Call Quit(_RC_INPUT_ERROR_)
      Endif

      If(YesNo(5).and..not.YesNo(6)) then
        Write(6,*)
        Write(6,*)' You have specified that initial coordinates are to'
     &//' be given in input, but no coordinates are found.'
        Write(6,*)' Provide these in the SOLVent section.'
        Call Quit(_RC_INPUT_ERROR_)
      Endif

      If(.not.YesNo(2).and..not.YesNo(7)) then
        Write(6,*)
        Write(6,*)' You fail to specify where from initial'
     &//' configuration should be collected.'
        Write(6,*)' Do this with the CONFiguration keyword.'
        Call Quit(_RC_INPUT_ERROR_)
      Endif

      If(YesNo(9).and..not.YesNo(10)) then
        Write(6,*)
        Write(6,*)' Your file specification implies that an'
     &//' extraction file is to be generated.'
        Write(6,*)' However, you have no EXTRact section.'
        Call Quit(_RC_INPUT_ERROR_)
      Endif

      If(.not.YesNo(9).and.YesNo(10)) then
        Write(6,*)
        Write(6,*)' You have a EXTRact section, but the file to'
     &//' read from is not a sampfile.'
        Write(6,*)' Change this after the FILE keyword.'
        Call Quit(_RC_INPUT_ERROR_)
      Endif

      Return
      End
