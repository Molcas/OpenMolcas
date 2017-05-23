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
      Subroutine DecideOnDirect(CanDoDirect,FoundTwoEls,DoDirect,
     &                          DoCholesky)
      Implicit Integer (a-z)
      Logical CanDoDirect,FoundTwoEls,DoDirect,DoCholesky
      Logical Direct,Expert,NeverDirect,AlwaysDirect
#include "OneDat.fh"
c----------------------------------------------------------------------
c                                                                     -
c  All modules using two-elecont integrals should call this routine   -
c  to decide whether to do the calculation integral-direct.           -
c                                                                     -
c  On input:     CanDoDirect : Direct capability of module (T/F).     -
c                FoundTwoEls : Two-electron integral file available.  -
c                                                                     -
c  On exit:      DoDirect:     Do calculation integral-direct.        -
c                DoCholesky:   Do Cholesky calculation                -
c                                                                     -
c----------------------------------------------------------------------

c  Read option variable set in Seward
*
*     Call Get_iOption(iOptSeward)

      Call DecideOnCholesky(DoCholesky)
      If (DoCholesky) then
              DoDirect=.false.
         return
      Endif

      Call Get_iScalar('System BitSwitch',iOptSeward)
*
      Direct=Iand(iOptSeward,1).Eq.1
      Expert=Iand(iOptSeward,2).Eq.2
      AlwaysDirect=Direct.And..Not.Expert
      NeverDirect=(.Not.Direct).And..Not.Expert

      If (AlwaysDirect) Then
        If (.Not.CanDoDirect) Then
          Write(6,'(A)')' Error, cannot do integral-direct calculation!'
          Write(6,'(A)')' Turn off DIRECT option in SEWARD input.'
          Call Abend()
        Else
          DoDirect=.True.
        Endif
      ElseIf ((.Not.FoundTwoEls) .And.
     &        (NeverDirect.Or..Not.CanDoDirect)) Then
*------ No integrals, no direct calculation (allowed) - we have to crash!
        Write(6,'(2A)')' Two-electron integral file was not found!'
        If(CanDoDirect) Write(6,'(A)')' Try keyword DIRECT in SEWARD.'
        Call QTrace()
        Call Abend()
      Else If (.Not. FoundTwoEls) Then
        DoDirect = .True.
      Else
        DoDirect = .False.
      Endif
*
      Return
      End
