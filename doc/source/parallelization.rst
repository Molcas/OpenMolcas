Parallelization efforts for |molcas| modules
============================================

Presented below is a table of modules in |molcas| that *can* benifit from
parallel execution through distribution of work and/or resources. If a module
is not listed in this table, and the module-specific documentation does not
mention anything about parallelization, then you have to assume the module is
not (efficiently) parallelized. This means that even though it will get
executed in parallel, all processes will perform the same serial calculation!
Be aware that for parallel modules with serial components, the use of the
serial components (indirectly or through the use of a keyword) might adversely
affect CPU and memory usage for large calculations. In that case, you might
have to increase the runtime or memory, or avoid/use keywords that
activate/deactivate the serial components.

.. table:: Modules in |molcas| which benefit from parallel processing.
   :name: tab:mpp_effort:

   ================= ============================================== =============================
   Module            Parallel speed-up expected for                 Notable non-parallel parts
   ================= ============================================== =============================
   :program:`seward` | conventional 2-el integrals                  | 1-el integrals
                     | Cholesky vectors                             | Douglas--Kroll--Hess
                                                                    | properties
   :program:`scf`    | orbital optimization                         | properties
   :program:`rasscf` | orbital optimization                         | CI optimization
                                                                    | properties
   :program:`mbpt2`
   :program:`caspt2` | Cholesky vectors                             | conventional 2-el integrals
                                                                    | properties
                                                                    | multi-state interaction
   :program:`alaska` | displacements (if using numerical gradients)
   :program:`geo`    | displacements
   ================= ============================================== =============================
