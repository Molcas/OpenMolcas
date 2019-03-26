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
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
*  DefinePairDomain
*
*> @brief
*>   Define pair domains
*> @author Thomas Bondo Pedersen
*>
*> @details
*> Pair domains are defined by the union of individual orbital
*> domains, \p iDomain(0:nAtom,nOcc), see SubRoutine ::DefineDomain.
*>
*> On exit, the contents of \p iPairDomain array are:
*>
*> - \p iPairDomain(0,ij): number of atoms in pair domain \c ij.
*> - \p iPairDomain(n,ij): id of atom \c n (\c n = ``1``, ``2``, ..., \p iPairDomain(0,ij)) in pair domain \c ij.
*>
*> The contents of \p Rmin array are:
*>
*> - \p Rmin(ij): minimum distance between atoms in pair \p ij.
*>
*> The contents of \p iClass array are (for \c i = ``1``, ``2`` , ``3``, ..., \p nRThr):
*>
*> - \p iClass(ij) = \c i-1 implies \p Rmin(ij) &le; \p RThr(i) and
*> - \p iClass(ij) = \p nRThr implies \p Rmin(ij) > \p RThr(nRThr)
*>
*> Note that if classification is not wanted, simply specify
*> \p nRThr = ``0`` (in which case arrays \p iClass and \p RThr are not referenced,
*> may be dummy arguments in the call to this routine).
*> The \p iDomain array should be set by SubRoutine ::DefineDomain
*> before calling this routine.
*>
*> Return codes:
*>
*> - \p irc = ``0``: all OK.
*>
*> (at the moment, this is the only possible return code, but that
*> might change.)
*>
*> @param[out] irc         Return code
*> @param[out] iPairDomain Pair domain definition
*> @param[out] iClass      Classification
*> @param[out] Rmin        Minimum interatomic distance in each pair
*> @param[in]  iDomain     Orbital domain definition
*> @param[in]  RThr        Thresholds for classification
*> @param[in]  Coord       Nuclear coordinates
*> @param[in]  nAtom       Number of atoms
*> @param[in]  nOcc        Number of orbitals for which domains are defined
*> @param[in]  nRThr       Number of thresholds for classification
************************************************************************
      SubRoutine DefinePairDomain(irc,iPairDomain,iClass,Rmin,iDomain,
     &                            RThr,Coord,nAtom,nOcc,nRThr)
      Implicit Real*8 (a-h,o-z)
      Integer iPairDomain(0:nAtom,*), iClass(*)
      Real*8  Rmin(*)
      Integer iDomain(0:nAtom,nOcc)
      Real*8  RThr(*), Coord(3,nAtom)
#include "WrkSpc.fh"

C     Set return code.
C     ----------------

      irc = 0
      If (nOcc .lt. 2) Return

C     Set pair domains as union of individual domains: [ij]=[i]U[j] for
C     i>=j.
C     -----------------------------------------------------------------

      nnOcc = nOcc*(nOcc+1)/2
      lT = (nAtom+1)*nnOcc
      Call iCopy(lT,[0],0,iPairDomain,1)

      l_Union = nAtom*nOcc
      Call GetMem('Union','Allo','Inte',ip_Union,l_Union)

      Call iCopy(l_Union,[0],0,iWork(ip_Union),1)
      kOff = ip_Union - 1
      Do i = 1,nOcc
         iOff = kOff + nAtom*(i-1)
         Do iA = 1,iDomain(0,i)
            iAtom = iDomain(iA,i)
            iWork(iOff+iAtom) = 1
         End Do
      End Do

      kOff = ip_Union - 1
      ij = 0
      Do j = 1,nOcc
         ij = ij + 1
         lD = iDomain(0,j) + 1
         Call iCopy(lD,iDomain(0,j),1,iPairDomain(0,ij),1) ! case i=j
         jOff = kOff + nAtom*(j-1)
         Do i = j+1,nOcc
            iOff = kOff + nAtom*(i-1)
            iCount = 0
            ij = ij + 1
            Do kAtom = 1,nAtom
               isThere = iWork(jOff+kAtom) + iWork(iOff+kAtom)
               If (isThere .gt. 0) Then
                  iCount = iCount + 1
                  iPairDomain(iCount,ij) = kAtom
               End If
            End Do
            iPairDomain(0,ij) = iCount
         End Do
      End Do

      Call GetMem('Union','Free','Inte',ip_Union,l_Union)

C     Set min. distance between any two atoms in pairs of domains.
C     ------------------------------------------------------------

      ij = 0
      Do j = 1,nOcc
         ij = ij + 1
         Rmin(ij) = 0.0d0 ! case i=j
         Do i = j+1,nOcc
            ij = ij + 1
            Rmin(ij) = 1.0d15
            Do jA = 1,iDomain(0,j)
               jAtom = iDomain(jA,j)
               Do iA = 1,iDomain(0,i)
                  iAtom = iDomain(iA,i)
                  R = sqrt((Coord(1,iAtom)-Coord(1,jAtom))**2
     &                    +(Coord(2,iAtom)-Coord(2,jAtom))**2
     &                    +(Coord(3,iAtom)-Coord(3,jAtom))**2)
                  Rmin(ij) = min(Rmin(ij),R)
               End Do
            End Do
         End Do
      End Do

C     Set classification for each pair domain.
C     Note that the caller must have ordered the thresholds according to
C     RThr(1) < RThr(2) < RThr(3) < ... < RThr(nRThr)
C     ------------------------------------------------------------------

      If (nRThr .gt. 0) Then
         Call iCopy(nnOcc,[nRThr],0,iClass,1)
         Do ij = 1,nnOcc
            i = 0
            Do While (i .lt. nRThr)
               i = i + 1
               If (Rmin(ij) .le. RThr(i)) Then
                  iClass(ij) = i - 1
                  i = nRThr ! break while loop
               End If
            End Do
         End Do
      End If

      End
