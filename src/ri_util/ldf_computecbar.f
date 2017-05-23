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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine LDF_ComputeCBar(iAtomPair,ip_CBar,l_CBar,ip_Z,l_Z,irc)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Input:
C        iAtomPair - index of atom pair in atom pair list
C
C     Output:
C        ip_CBar, l_CBar - pointer and dimension of array
C                          containing CBar
C        ip_Z, l_Z - pointer and dimension of array containing Z
C        irc - return code
C
C     Return codes:
C        irc=-1: input error (iAtomPair out of bounds)
C        irc=0:  all is fine
C        irc=1:  internal error
C
C     Purpose:
C
C     Compute
C
C        CBar[uv,J] = { (uv|J) - sum(I=1,J-1) CBar[uv,I]*Z[J,I] }/Z[J,J]
C
C     where
C
C        sum(K) Z[I,K]*Z[J,K] = G[I,J] = (I|J)
C
C     and u and v are basis functions centered on atom iAtom and jAtom,
C     respectively. (The atom pair iAtomPair implies iAtom and jAtom.)
C     The auxiliary basis functions are taken from the
C     union of the sets on atom iAtom and atom jAtom plus two-center
C     functions (if any, and if iAtom != jAtom). Z vectors are computed
C     locally and discarded. If iAtom=jAtom, the Cbar array is packed
C     according to u>=v (triangular storage), else it is rectangular.
C
C     If linear dependence in the auxiliary basis is detected, the
C     linearly dependent functions are removed (atom pair info will be
C     updated to reflect the changes).
C
      Implicit None
      Integer iAtomPair
      Integer ip_CBar, l_CBar
      Integer ip_Z, l_Z
      Integer irc
#include "ldf_atom_pair_info.fh"
#include "WrkSpc.fh"

      Character*15 SecNam
      Parameter (SecNam='LDF_ComputeCBar')

      Integer iAtom, jAtom
      Integer nBas_iAtom, nBas_jAtom
      Integer M
      Integer ip_G, l_G
      Integer nuv, ip_this_vec

      Real*8 ZJJ_inv

      Integer  LDF_nBas_Atom, LDF_nBasAux_Pair
      External LDF_nBas_Atom, LDF_nBasAux_Pair

      Integer i, j
      Integer iTri, AP_Atoms
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      ! Initializations
      irc=0
      ip_CBar=0
      l_CBar=0
      ip_Z=0
      l_Z=0
#if defined (_DEBUG_)
      If (iAtomPair.lt.1 .or. iAtomPair.gt.NumberOfAtomPairs) Then
         Write(6,'(A,A,I8)')
     &   SecNam,': iAtomPair out of bounds:',iAtomPair
         irc=-1
         Return
      End If
#endif

      ! Get atoms corresponding to atom pair iAtomPair
      iAtom=AP_Atoms(1,iAtomPair)
      jAtom=AP_Atoms(2,iAtomPair)

      ! Get number of auxiliary basis functions
      M=LDF_nBasAux_Pair(iAtomPair)

      ! Get number of basis functions on each atom
      nBas_iAtom=LDF_nBas_Atom(iAtom)
      nBas_jAtom=LDF_nBas_Atom(jAtom)
      nuv=nBas_iAtom*nBas_jAtom

      ! Allocate CBar array
      l_CBar=nuv*M
      Call GetMem('CBar','Allo','Real',ip_CBar,l_CBar)

      ! Allocate G matrix
      l_G=M*M
      Call GetMem('GMatrix','Allo','Real',ip_G,l_G)

      ! Set index arrays (IndxG, IndxG2, 2CList)
      Call LDF_SetIndxG(iAtomPair)

      ! Compute G matrix
      Call LDF_ComputeGMat(iAtomPair,M,Work(ip_G))

      ! Compute Z vectors by CD of G
      ! (CBar may be used as work space)
      ! Linear dependence is handled here.
      ! Z vectors are stored as lower triangle.
      Call LDF_ComputeZVec(iAtomPair,ip_CBar,l_Cbar,ip_G,l_G,ip_Z,l_Z,
     &                     irc)
      If (irc.ne.0) Then
         Write(6,'(A,A,I8)')
     &   SecNam,': LDF_ComputeZVec returned code',irc
         irc=1
         Call LDF_UnsetIndxG()
         Call GetMem('GMatrix','Free','Real',ip_G,l_G)
         Call GetMem('CBar','Free','Real',ip_CBar,l_CBar)
         Return
      End If

      ! Unset index arrays for G
      Call LDF_UnsetIndxG()

      ! Deallocate G matrix
      Call GetMem('GMatrix','Free','Real',ip_G,l_G)

      ! Resize CBar array if needed
      If (LDF_nBasAux_Pair(iAtomPair).lt.M) Then
         Call GetMem('CBar','Free','Real',ip_CBar,l_CBar)
         M=LDF_nBasAux_Pair(iAtomPair)
         l_CBar=nuv*M
         Call GetMem('CBar','Allo','Real',ip_CBar,l_CBar)
      End If

      ! Set index arrays for G, including newly detected lin dep
      Call LDF_SetIndxG(iAtomPair)

      ! Compute integrals (uv|J), store in CBar
      Call LDF_ComputeIntegrals_uvJ(iAtomPair,l_Cbar,Work(ip_CBar))

      ! Unset index arrays for G
      Call LDF_UnsetIndxG()

      ! Compute CBar:
      ! CBar[uv,J] = { (uv|J) - sum(I=1,J-1) CBar[uv,I]*Z[J,I] }/Z[J,J]
      Do J=1,M
         ! Compute CBar[uv,J]
         ZJJ_inv=1.0d0/Work(ip_Z-1+iTri(J,J))
         ip_this_vec=ip_CBar+nuv*(J-1)
         Call dScal_(nuv,ZJJ_inv,Work(ip_this_vec),1)
         ! Subtract from all subsequent integrals
         Do I=J+1,M
            Call dAXPY_(nuv,-Work(ip_Z-1+iTri(I,J)),
     &                 Work(ip_this_vec),1,
     &                 Work(ip_CBar+nuv*(I-1)),1)
         End Do
      End Do

      End
