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
      Subroutine Funi_Print
      Implicit Real*8 (A-H,O-Z)
#include "nq_info.fh"
#include "print.fh"
      Logical Check
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function
*
      Check(i,j)=iAnd(i,2**(j-1)).ne.0
*                                                                      *
************************************************************************
*                                                                      *
      Call qEnter('Funi_Print')
      iRout=666
      iPrint=nPrint(iRout)
*                                                                      *
************************************************************************
*                                                                      *
      Call Get_dScalar('EThr',EThr)
      T_X=Min(T_X,EThr*1.0D-2)
      T_Y=Min(T_Y,EThr*1.0D-2)
      Threshold=Min(Threshold,EThr*1.0D-4)
*                                                                      *
************************************************************************
*                                                                      *
      If (iPrint.ge.10) Then
      Write (6,*)
      Write (6,'(6X,A)') 'Fast Utility for Numerical Integration'
      Write (6,'(6X,A)') 'Authors: R. Lindh, Y. Carissan, '
     &                 //'M. Wierzbowska and L. Gagliardi.'
      End If
*
      If (iPrint.ge.6) Then
      Write (6,*)
      Write (6,'(6X,A)') 'Numerical integration parameters'
      Write (6,'(6X,A)') '======================================'
      Write (6,'(6X,A,21X,A)') 'Radial quadrature type:    ',Quadrature
*
      If (Quadrature(1:3).eq.'LMG') Then
         Write (6,'(6X,A,E11.4)') 'Radial quadrature accuracy:',
     &                             Threshold
      Else
         Write (6,'(6X,A,18X,I5)')    'Size of radial grid:       ',
     &                             nR
      End If
*
      If (Check(iOpt_Angular,3)) Then
         Write (6,'(6X,A,25X,I4)') 'Lebedev angular grid:',L_Quad
      Else If (Check(iOpt_Angular,1)) Then
         Write (6,'(6X,A,I4)') 'Lobatto angular grid, l_max:',L_Quad
      Else
         Write (6,'(6X,A,I4)')
     &     'Gauss and Gauss-Legendre angular grid, l_max:',L_Quad
      End If
*
      If (Angular_Prunning.eq.On) Then
         Write (6,'(6X,A,1X,ES9.2)')
     &         'Angular grid prunned with the crowding factor:',
     &         Crowding
         Write (6,'(6X,A,1X,ES9.2)')
     &         '                            and fading factor:',
     &         Fade
      End If
      If (Check(iOpt_Angular,2)) Then
         Write (6,'(6X,A)') 'The whole atomic grid is scanned for each'
     &                    //' sub block.'
      End If
*
      Write (6,'(6X,A,2X,ES9.2)')
     &      'Screening threshold for integral computation:',
     &      T_Y
      Write (6,'(6X,A,2X,ES9.2)')
     &      'Screening threshold for density computation: ',
     &      T_X
      If (Quadrature(1:3).ne.'LMG') Then
         Write (6,'(6X,A,20X,ES9.2)') 'Radial quadrature accuracy:',
     &                             Threshold
      End If
*
      Write (6,'(6X,A,17X,I7)') 'Maximum batch size:        ',nGridMax
      If (NQ_Direct.eq.On) Then
         Write (6,'(6X,A)') 'AO values are recomputed each iteration'
      Else
         Write (6,'(6X,A)') 'AO values are stored on disk'
      End If
      Write (6,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Put flag on RUNFILE to indicate that we are doing DFT.
*
*     Call Get_iOption(iOpt)
      Call Get_iScalar('System BitSwitch',iOpt)
      iOpt=iOr(iOpt,2**6)
*     Call Put_iOption(iOpt)
      Call Put_iScalar('System BitSwitch',iOpt)
*                                                                      *
************************************************************************
*                                                                      *
      Call qExit('Funi_Print')
      Return
      End
