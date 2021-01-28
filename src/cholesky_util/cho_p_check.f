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
      SubRoutine Cho_P_Check(irc)
C
C     Purpose: check settings for parallel Cholesky.
C
C-TODO/FIXME: The features not allowed in parallel execution of the
C             Cholesky decomposition should be implemented later. Thus,
C             this subroutine is, in effect, a TODO-list.
C
      Use Para_Info, Only: nProcs, Is_Real_Par
      use ChoSubScr, only: Cho_SScreen
      Implicit None
      Integer irc
#include "cholesky.fh"
#include "chosimri.fh"
#include "cho_para_info.fh"

      Logical WriteBlank

      irc = 0
      WriteBlank=.True.

      If (Cho_Real_Par) Then ! TRUE PARALLEL
         If (Cho_DecAlg.ne.4 .and. Cho_DecAlg.ne.5 .and.
     &       Cho_DecAlg.ne.6) Then
            If (WriteBlank) Then
               Write(Lupri,*)
               WriteBlank=.False.
            End If
            Write(Lupri,'(A,A)')
     &      'Only possible parallel Cholesky decomposition algorithm ',
     &      'is "PARAllel".'
            Write(Lupri,'(A,I3,A)')
     &                     'Resetting Cho_DecAlg from ',Cho_DecAlg,
     &                     ' to 5 (parallel two-step algorithm),'
            Cho_DecAlg = 5
         End If
         If (MxShPr .ne. 1) Then
            If (WriteBlank) Then
               Write(Lupri,*)
               WriteBlank=.False.
            End If
            Write(Lupri,'(A,A,A)')
     &                     'Max. number of shell pair distributions ',
     &                     'calculated in each pass is 1 for ',
     &                     'parallel Cholesky.'
            Write(Lupri,'(A,I6,A)')
     &                     'Resetting MxShPr from ',MxShPr,' to 1'
            MxShPr = 1
         End If
         If (Cho_IntChk) Then
            If (WriteBlank) Then
               Write(Lupri,*)
               WriteBlank=.False.
            End If
            Write(Lupri,'(A)') 'You have requested integral checking.'
            Write(Lupri,'(A,A)')
     &                     'Integral checking is not possible for ',
     &                     'parallel Cholesky.'
            irc = irc + 1
         End If
         If (RstDia .or. RstCho) Then
            If (WriteBlank) Then
               Write(Lupri,*)
               WriteBlank=.False.
            End If
            If (RstDia) Then
               Write(Lupri,'(A)') 'You have requested diagonal restart.'
               irc = irc + 1
            End If
            If (RstCho) Then
               Write(Lupri,'(A)')
     &         'You have requested decomposition restart.'
               irc = irc + 1
            End If
            Write(Lupri,'(A,A)')
     &                     'Restart is not possible for parallel ',
     &                     'Cholesky.'
         End If
         If (Cho_ReOrd) Then
            If (WriteBlank) Then
               Write(Lupri,*)
               WriteBlank=.False.
            End If
            Write(Lupri,'(A,A)')
     &                     'Vector reordering is not possible for ',
     &                     'parallel Cholesky.'
            irc = irc + 1
         End If
         If (Cho_AdrVec .ne. 1) Then
            If (WriteBlank) Then
               Write(Lupri,*)
               WriteBlank=.False.
            End If
            Write(Lupri,'(A,A)')
     &                     'Address mode for vector I/O must be word-',
     &                     'addressable for parallel Cholesky.'
            Write(Lupri,'(A,I4,A)')
     &      'Resetting Cho_AdrVec from ',Cho_AdrVec,' to 1'
            Cho_AdrVec = 1
         End If
         If (IfcSew .ne. 2) Then
            If (WriteBlank) Then
               Write(Lupri,*)
               WriteBlank=.False.
            End If
            Write(Lupri,'(A,A)')
     &      'Seward interface must be directly in reduced ',
     &      'sets for parallel Cholesky.'
            Write(Lupri,'(A,I4,A)')
     &      'Resetting IfcSew from ',IfcSew,' to 2'
            IfcSew = 2
         End If
         If (Cho_TstScreen) Then
            If (WriteBlank) Then
               Write(Lupri,*)
               WriteBlank=.False.
            End If
            Write(Lupri,'(A,A)')
     &                     'Test of subtraction screening is not ',
     &                     'possible for parallel Cholesky.'
            Write(Lupri,'(A)') 'Turning Cho_TstScreen off.'
            Cho_TstScreen = .False.
         End If
         If (Cho_SScreen) Then
            If (WriteBlank) Then
               Write(Lupri,*)
               WriteBlank=.False.
            End If
            Write(Lupri,'(A,A)') 'Subtraction screening is not ',
     &                     'possible for parallel Cholesky.'
            irc = irc + 1
         End If
         If (Cho_SimRI) Then
            If (WriteBlank) Then
               Write(Lupri,*)
               WriteBlank=.False.
            End If
            Write(Lupri,'(A,A)') 'Simulation of RI is not ',
     &                     'possible for parallel Cholesky.'
            irc = irc + 1
         End If
      Else
         If (CHO_FAKE_PAR .and.
     &       nProcs.gt.1 .and. Is_Real_Par()) Then ! FAKE PARALLEL
            If (Cho_ReOrd) Then
               If (WriteBlank) Then
                  Write(Lupri,*)
                  WriteBlank=.False.
               End If
               Write(Lupri,'(A,A)')
     &                        'Vector reordering is not possible for ',
     &                        'parallel Cholesky.'
               irc = irc + 1
            End If
         End If
      End If

      End
