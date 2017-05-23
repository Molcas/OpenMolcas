************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1996, Markus P. Fuelscher                              *
************************************************************************
      Subroutine CiOvlp(jRoot,S1,S2,CI_vec)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Compute the overlap of the CI_vector, CI_vec, with the           *
*     a set of test vectors.                                           *
*                                                                      *
*     calling arguments:                                               *
*     jRoot   : integer                                                *
*               root identifier                                        *
*     S1      : array of real*8, input/output                          *
*               overlap matrix with test vectors                       *
*     S2      : array of real*8, input/output                          *
*               norm of the test configurations in the CI vector       *
*     CI_vec  : array of real*8, input                                 *
*               CI_vector                                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Real*8 (A-H,O-Z)


#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"

      Dimension S1(lRoots,lRoots)
      Dimension S2(lRoots,lRoots)
      Dimension CI_vec(nConf)

      If ( ITER.eq.1 ) Return

      Do kRoot = 1,nRoots
        Sum1 = 0.0d0
        Sum2 = 0.0d0
        Do kRef = 1,mxRef
          iConf = jCj(kRoot,kRef)
          If ( iConf.ne.0 ) then
            Sum1 = Sum1 + cCI(kRoot,kRef)*CI_vec(iConf)
            Sum2 = Sum2 + CI_vec(iConf)*CI_vec(iConf)
          End If
        End Do
        S1(jRoot,kRoot) = Abs(Sum1)
        S2(jRoot,kRoot) = Sum2
      End Do

      Return
      End

      Subroutine CiSelect(S1,S2)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Select CI_vector which matches best the test vectors             *
*                                                                      *
*     calling arguments:                                               *
*     S1      : array of real*8, input/output                          *
*               overlap matrix with test vectors                       *
*     S2      : array of real*8, input/output                          *
*               norm of the test configurations in the CI vector       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Real*8 (A-H,O-Z)


#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"

      Dimension S1(lRoots,lRoots)
      Dimension S2(lRoots,lRoots)

      Dimension iTemp(mxRoot)

      If ( ITER.eq.1 ) Return

*     make a local copy of the present selection vector of weights
      Call iCopy(mxRoot,0,1,iTemp,1)

*     make a new choice using the overlap ov the etst vector with
*     the CI vector in the subspace of the test vector
      Do kRoot = 1,nRoots
*...    search for the element of largest overlap
        maxS1  = 1
        S1max  = S1(1,kRoot)
        Do jRoot = 1,lRoots
          If( S1(jRoot,kRoot).gt.S1max ) then
            maxS1 = jRoot
            S1max = S1(jRoot,kRoot)
          End If
        End Do
        iTemp(kRoot) = maxS1
*...    cleanup,to ensure that we don't pick the same root twice
        Do jRoot = 1,nRoots
          S1(maxS1,jRoot) = S1(maxS1,jRoot)-999999.0d0
        End Do
      End Do

*     Restore S1
      Do kRoot = 1,nRoots
        Do jRoot = 1,nRoots
          S1(iTemp(kRoot),jRoot) = S1(iTemp(kRoot),jRoot)+999999.0d0
        End Do
      End Do

*     print
      If ( nRoots.eq.1 ) then
        Write (6,'(6X,A,T45,10I6)')
     &   'new root selected:',(iTemp(kRoot),kRoot=1,nRoots)
        Write (6,'(6X,A,T45,10F6.3)')
     &   'overlap           ',(S1(iTemp(kRoot),kRoot),kRoot=1,nRoots)
      Else
        Write (6,'(6X,A,T45,10I6)')
     &   'new roots selected:',(iTemp(kRoot),kRoot=1,nRoots)
        Write (6,'(6X,A,T45,10F6.3)')
     &   'overlap           ',(S1(iTemp(kRoot),kRoot),kRoot=1,nRoots)
      End If

*     Compare the overlap elemets <1|2> and <1|1> where
*     |1> is the CI vector in the subspace of the test vector and
*     |2> is the test vector.
*     The program breaks execution if the following roules apply:
*     If <1|2> is less than 0.5*<1|1>
*     If the weight of <1|2> or <1|1> is less than 0.3
      istop = 0
      icase1 = 1
      icase2 = 2
      icase3 = 4
      Do kRoot = 1,nRoots
        S1jk = S1(iTemp(kRoot),kRoot)
        S2jk = S2(iTemp(kRoot),kRoot)
        If ( S1jk.lt.0.5d0*S2jk ) istop = IOR(istop,icase1)
        If ( Sqrt(S1jk).lt.0.316d0 ) istop = IOR(istop,icase2)
        If ( Sqrt(S2jk).lt.0.3d0 ) istop = IOR(istop,icase3)
      End Do

*     If the stop flag has been set write an approriate message
*     and to stop execution change the iteration counter.
      If ( istop.ge.1 ) then
        Write (6,*)
        Write (6,'(6X,120A1)') ('=',i=1,120)
        If ( IAND(istop,icase1).ne.0 ) then
          Write (6,'(6X,A)')
     &    'The projection of the CI vector(s) onto the model vector(s)'
          Write (6,'(6X,A)')
     &    'is smaller than half the norm of the subspace.'
        End If
        If ( IAND(istop,icase2).ne.0 ) then
          Write (6,'(6X,A)')
     &    'The overlap of the projected CI vector(s) and the model '//
     &    'vector(s) is smaller than 0.1'
        End If
        If ( IAND(istop,icase3).ne.0 ) then
          Write (6,'(6X,A)')
     &    'The weight(s) of the subspace is(are) smaller than 30% '//
     &    'of the total wave function(s)'
        End If
        Write (6,'(6X,A)')
     &  'Please, check your model space'
        Write (6,'(6X,A)')
     &  'The program stops after the next iteration'
        Write (6,'(6X,120A1)') ('=',i=1,120)
        Write (6,*)
        Write (6,*)
        MAXIT=ITER
      Else
        Call iCopy(nRoots,iTemp,1,iRoot,1)
      End If

      Return
      End
