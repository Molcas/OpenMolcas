.. _UG\:sec\:structure_of_the_ecp_libraries:

Structure of the ECP libraries
..............................

.. compound::

  The start of a given basis set and AIMP is identified by the line ::

    /label

  where "label" is defined below,
  in the input description to :program:`seward`.
  Then, comment lines, effective charge, and basis set follow,
  with the same structure that the all-electron Basis Set Library
  (see items 1. to 4. in Sec. :ref:`UG:sec:structure_of_the_ae_basis_set_libraries`.)
  Next, the AIMP/ECP/PP is specified as follows:

5. The pseudo potential approach :cite:`Kahn:68,Christiansen:79,Durand:75`,
   see eqs. (3) and (4) in Ref. :cite:`Skylaris:98`,
   with the following lines:

   #. The keyword :kword:`PP`
      On the same line follows the atomic symbol of the element, the number of core electrons (:math:`N_c`) and
      :math:`L`, where :math:`L-1` is the largest angular momentum orbital belonging to the core. This line
      is followed by :math:`L+1` identical sections. The first of these sections is the so-called :math:`L` potential
      and the subsequent sections corresponds to the S-\ :math:`L`, P-\ :math:`L`, D-\ :math:`L`, etc. potentials.
      Each sections start with a line specifying the number of Gaussian terms in the potential.
      This line is then followed by a single line for each Gaussian specifying
      the powers (:math:`n_{kl}`), the Gaussian exponent (:math:`\zeta_{kl}`), and the associated coefficient
      (:math:`d_{kl}`).

      Note that the pseudo potential input is mutually exclusive to the :kword:`M1`, :kword:`M2`,
      :kword:`COREREP`, and :kword:`PROJOP` keywords!

#. The Coulomb local model potential, eq.(6) in Ref. :cite:`Huzinaga:87`
   with the following lines:

   #. The keyword :kword:`M1`,
      which identifies the terms with :math:`n_k=0`.
   #. The number of terms.
      If greater than 0, lines :ref:`3 <it:alphak0>` and :ref:`4 <it:Ak0>` are read.

      .. _it\:alphak0:

   #. The exponents :math:`\alpha_k`.

      .. _it\:Ak0:

   #. The coefficients :math:`A_k`.
      (divided by the negative of the effective charge).
   #. The keyword :kword:`M2`,
      which identifies the terms with :math:`n_k=1`.
   #. The number of terms.
      If greater than 0, lines :ref:`7 <it:alphak1>` and :ref:`8 <it:Ak1>` are read.

      .. _it\:alphak1:

   #. The exponents :math:`\alpha_k`.

      .. _it\:Ak1:

   #. The coefficients :math:`A_k`.
      (divided by the negative of the effective charge).

#. A line with the keyword :kword:`COREREP`
   followed by another one with a real constant.
   This is not used now but it is reserved for future use.
#. The projection operator, eq.(3) in Ref. :cite:`Huzinaga:87`
   with the following lines:

   #. The keyword :kword:`PROJOP`.
   #. The maximum angular momentum (l)
      of the frozen core (embedding) orbitals.
      Lines :ref:`3 <it:nPrim>` to :ref:`6 <it:coeffs>` are repeated for
      each angular momentum l.

      .. _it\:nPrim:

   #. The number of primitives and the number of orbitals
      (more properly, degenerate sets of orbitals or l-shells)
      for angular momentum l.
      As an option, these two integers can be followed by the occupation numbers
      of the l-shells; default values are 2 for l=0, 6 for l=1, etc.
   #. The projection constants, :math:`-2\varepsilon_c`.
   #. The exponents of the primitive functions.

      .. _it\:coeffs:

   #. The coefficients of the orbitals, one per column,
      using general contraction format.

#. The spectral representation operator,
   eq.(7) in Ref. :cite:`Huzinaga:87` for NR-AIMP,
   eq.(3) in Ref. :cite:`Barandiaran:90` for relativistic CG-AIMP,
   and eqs.(1) and (7) in Ref. :cite:`Rakowitz:99b` for relativistic NP-AIMP,
   with the following lines:

   #. The keyword :kword:`Spectral Representation Operator`.
   #. One of the keywords
      :kword:`Valence`, :kword:`Core`, or :kword:`External`.
      :kword:`Valence` indicates that the set of primitive functions specified
      in the basis set data will be used for the spectral representation
      operator;
      this is the standard for ab initio *core* model potentials.
      :kword:`Core` means that the set of primitives specified
      in the :kword:`PROJOP` section will be used instead;
      this is the standard for complete-ion ab initio *embedding* model
      potentials.
      :kword:`External` means that a set of primitives specific for the
      spectral representation operator will be provided in the next lines.
      In this case the format is
      one line in which an integer number specifies the highest angular momentum
      of the external basis sets; then, for each angular momentum the input is
      formated as for lines
      :ref:`2 <it:ae-2>`, :ref:`3 <it:ae-3>`, and :ref:`4 <it:ae-4>`
      in Sec. :ref:`UG:sec:structure_of_the_ae_basis_set_libraries`.
   #. The keyword :kword:`Exchange`.
   #. For relativistic AIMPs one of the keywords :kword:`NoPair` or :kword:`1stOrder Relativistic Correction`.
      :kword:`NoPair` indicates that scalar relativistic no-pair Douglas-Kroll AIMP integrals are to be calculated.
      :kword:`1stOrder Relativistic Correction` means that Cowan-Griffin-based scalar relativistic AIMP, CG-AIMP's, are used.
      In the latter case, in the next line a *keyword* follows which, in the library :file:`QRPLIB`,
      identifies
      the starting of the numerical mass-velocity plus Darwin potentials
      (eq.(2) in Ref. :cite:`Barandiaran:90`).
      (In :file:`QRPLIB` a line with
      "*keyword* :kword:`mv&dw potentials start`" must
      exist, followed by the number of points in the radial logarithmic
      grid, the values of the radial coordinate r, and,
      for each valence orbital, its label (2S, 4P, etc.),
      and the values of the mass-velocity plus Darwin potentials
      at the corresponding values of r; these data must end up with
      a line
      "*keyword* :kword:`mv&dw potentials end`".)
   #. The keyword :kword:`End of Spectral Representation Operator`.

Below is an example of an entry in the :file:`ECP` library for an AIMP. ::

  /S.ECP.Barandiaran.7s6p1d.1s1p1d.6e-CG-AIMP.        -- label (note that type is ECP)
  Z.Barandiaran and L.Seijo, Can.J.Chem. 70(1992)409. -- 1st ref. line
  core[Ne] val[3s,3p]  (61/411/1*)=2s3p1d recommended -- 2nd ref. line
  *SQR-SP(7/6/1)                 (61/411/1)           -- comment line
    6.000000         2                                -- eff. charge & highest ang.mom.
                                                      -- blank line
      7    1                                          -- 7s -> 1s
     1421.989530                                      -- s-exponent
     211.0266560                                      -- s-exponent
     46.72165060                                      -- s-exponent
     4.310564040                                      -- s-exponent
     1.966475840                                      -- s-exponent
     .4015383790                                      -- s-exponent
     .1453058790                                      -- s-exponent
     .004499703540                                    -- contr. coeff.
     .030157124800                                    -- contr. coeff.
     .089332590700                                    -- contr. coeff.
    -.288438151000                                    -- contr. coeff.
    -.279252515000                                    -- contr. coeff.
     .700286615000                                    -- contr. coeff.
     .482409523000                                    -- contr. coeff.
      6    1                                          -- 6p -> 1p
     78.08932440                                      -- p-exponent
     17.68304310                                      -- p-exponent
     4.966340810                                      -- p-exponent
     .5611646780                                      -- p-exponent
     .2130782690                                      -- p-exponent
     .8172415400E-01                                  -- p-exponent
    -.015853278200                                    -- contr. coeff.
    -.084808963800                                    -- contr. coeff.
    -.172934245000                                    -- contr. coeff.
     .420961662000                                    -- contr. coeff.
     .506647309000                                    -- contr. coeff.
     .200082121000                                    -- contr. coeff.
      1    1                                          -- 1d -> 1d
     .4210000000                                      -- d-exponent
    1.000000000000                                    -- contr. coeff.
  *                                                   -- comment line
  * Core AIMP: SQR-2P                                 -- comment line
  *                                                   -- comment line
  * Local Potential Parameters : (ECP convention)    -- comment line
  *                            A(AIMP)=-Zeff*A(ECP)   -- comment line
  M1                                                  -- M1 operator
      9                                               -- number of M1 terms
     237485.0100                                      -- M1 exponent
     24909.63500                                      -- M1 exponent
     4519.833100                                      -- M1 exponent
     1082.854700                                      -- M1 exponent
     310.5610000                                      -- M1 exponent
     96.91851000                                      -- M1 exponent
     26.63059000                                      -- M1 exponent
     9.762505000                                      -- M1 exponent
     4.014487500                                      -- M1 exponent
                                                      -- blank line
     .019335998333                                    -- M1 coeff.
     .031229360000                                    -- M1 coeff.
     .061638463333                                    -- M1 coeff.
     .114969451667                                    -- M1 coeff.
     .190198283333                                    -- M1 coeff.
     .211928633333                                    -- M1 coeff.
     .336340950000                                    -- M1 coeff.
     .538432350000                                    -- M1 coeff.
     .162593178333                                    -- M1 coeff.
  M2                                                  -- M2 operator
      0                                               -- number of M2 terms
  COREREP                                             -- CoreRep operator
     1.0                                              -- CoreRep constant
  PROJOP                                              -- Projection operator
      1                                               -- highest ang. mom.
      8    2                                          -- 8s -> 2s
    184.666320      18.1126960                        -- 1s,2s proj. op. constants
     3459.000000                                      -- s-exponent
     620.3000000                                      -- s-exponent
     171.4000000                                      -- s-exponent
     58.53000000                                      -- s-exponent
     22.44000000                                      -- s-exponent
     6.553000000                                      -- s-exponent
     2.777000000                                      -- s-exponent
     1.155000000                                      -- s-exponent
     .018538249000   .005054826900                    -- contr. coeffs.
     .094569248000   .028197248000                    -- contr. coeffs.
     .283859290000   .088959130000                    -- contr. coeffs.
     .454711270000   .199724180000                    -- contr. coeffs.
     .279041370000   .158375340000                    -- contr. coeffs.
     .025985763000  -.381198090000                    -- contr. coeffs.
    -.005481472900  -.621887210000                    -- contr. coeffs.
     .001288714400  -.151789890000                    -- contr. coeffs.
      7    1                                          -- 7p -> 1p
    13.3703160                                        -- 2p proj. op. constant
     274.0000000                                      -- p-exponent
     70.57000000                                      -- p-exponent
     24.74000000                                      -- p-exponent
     9.995000000                                      -- p-exponent
     4.330000000                                      -- p-exponent
     1.946000000                                      -- p-exponent
     .8179000000                                      -- p-exponent
     .008300916100                                    -- cont. coeff.
     .048924254000                                    -- cont. coeff.
     .162411660000                                    -- cont. coeff.
     .327163550000                                    -- cont. coeff.
     .398615170000                                    -- cont. coeff.
     .232548200000                                    -- cont. coeff.
     .034091088000                                    -- cont. coeff.
  *                                                   -- comment line
  Spectral Representation Operator                    -- SR operator
  Valence primitive basis                             -- SR basis specification
  Exchange                                            -- Exchange operator
  1stOrder Relativistic Correction                    -- mass-vel + Darwin oper.
  SQR-2P                                              -- label in QRPLIB
  End of Spectral Representation Operator             -- end of SR operator

Below is an example of an entry in the :file:`ECP` library for a pseudo potential. ::

  /Hg.ECP.Dolg.4s4p2d.2s2p1d.2e-MWB                   -- label (note the type ECP)
  W. Kuechle, M. Dolg, H. Stoll, H. Preuss, Mol. Phys.-- ref. line 1
  74, 1245 (1991); J. Chem. Phys. 94, 3011 (1991).    -- ref. line 2
      2.00000    2                                    -- eff. charge & highest ang.mom.
  *s functions                                        -- comment line
    4  2                                              -- 4s -> 2s
    0.13548420E+01                                    -- s-exponent
    0.82889200E+00                                    -- s-exponent
    0.13393200E+00                                    -- s-exponent
    0.51017000E-01                                    -- s-exponent
    0.23649400E+00  0.00000000E+00                    -- contr. coeff.
   -0.59962800E+00  0.00000000E+00                    -- contr. coeff.
    0.84630500E+00  0.00000000E+00                    -- contr. coeff.
    0.00000000E+00  0.10000000E+01                    -- contr. coeff.
  *p functions                                        -- comment line
    4  2                                              -- 4p -> 2p
    0.10001460E+01                                    -- p-exponent
    0.86645300E+00                                    -- p-exponent
    0.11820600E+00                                    -- p-exponent
    0.35155000E-01                                    -- p-exponent
    0.14495400E+00  0.00000000E+00                    -- contr. coeff.
   -0.20497100E+00  0.00000000E+00                    -- contr. coeff.
    0.49030100E+00  0.00000000E+00                    -- contr. coeff.
    0.00000000E+00  0.10000000E+01                    -- contr. coeff.
  *d functions                                        -- comment line
    1  1                                              -- 1d -> 1d
    0.19000000E+00                                    -- d-exponent
    0.10000000E+01                                    -- contr. coeff.
  *                                                   -- comment line
  PP,Hg,78,5;                                         -- PP operator, label, # of core elec., L
  1; ! H POTENTIAL                                    -- # number of exponents in the H potential
  2, 1.00000000,.000000000;                           -- power, exponent and coeff.
  3; ! S-H POTENTIAL                                  -- # number of exponents in the S-H potential
  2,0.227210000,-.69617800;                           -- power, exponent and coeff.
  2, 1.65753000,27.7581050;                           -- power, exponent and coeff.
  2, 10.0002480,48.7804750;                           -- power, exponent and coeff.
  2; ! P-H POTENTIAL                                  -- # number of exponents in the P-H potential
  2,0.398377000,-2.7358110;                           -- power, exponent and coeff.
  2,0.647307000,8.57563700;                           -- power, exponent and coeff.
  2; ! D-H POTENTIAL                                  -- # number of exponents in the D-H potential
  2,0.217999000,-.01311800;                           -- power, exponent and coeff.
  2,0.386058000,2.79286200;                           -- power, exponent and coeff.
  1; ! F-H POTENTIAL                                  -- # number of exponents in the F-H potential
  2,0.500000000,-2.6351640;                           -- power, exponent and coeff.
  1; ! G-H POTENTIAL                                  -- # number of exponents in the G-H potential
  2,0.800756000,-13.393716;                           -- power, exponent and coeff.
  *                                                   -- comment line
  Spectral Representation Operator                    -- SR operator
  End of Spectral Representation Operator             -- end of SR operator
