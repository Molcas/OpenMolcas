.. table:: Symmetries available in |molcas|, including generators and symmetry elements.
   :name: tab:symmetry_symbol

   ============== ============== ================ =========== ========= ============== ================ ================= =========== =================== =================== ===================
   Group          Generators                                  Elements
   -------------- ------------------------------------------- -----------------------------------------------------------------------------------------------------------------------------------
   |zws|          :math:`g_1`    :math:`g_2`      :math:`g_3` :math:`E` :math:`g_1`    :math:`g_2`      :math:`g_1g_2`    :math:`g_3` :math:`g_1g_3`      :math:`g_2g_3`      :math:`g_1g_2g_3`
   ============== ============== ================ =========== ========= ============== ================ ================= =========== =================== =================== ===================
   :math:`C_1`                                                :math:`E`
   :math:`C_2`    :math:`C_2`                                 :math:`E` :math:`C_2`
   :math:`C_s`    :math:`\sigma`                              :math:`E` :math:`\sigma`
   :math:`C_i`    :math:`i`                                   :math:`E` :math:`i`
   :math:`C_{2v}` :math:`C_2`    :math:`\sigma_v`             :math:`E` :math:`C_2`    :math:`\sigma_v` :math:`\sigma_v'`
   :math:`C_{2h}` :math:`C_2`    :math:`i`                    :math:`E` :math:`C_2`    :math:`i`        :math:`\sigma_h`
   :math:`D_2`    :math:`C_2^z`  :math:`C_2^y`                :math:`E` :math:`C_2^z`  :math:`C_2^y`    :math:`C_2^x`
   :math:`D_{2h}` :math:`C_2^z`  :math:`C_2^y`    :math:`i`   :math:`E` :math:`C_2^z`  :math:`C_2^y`    :math:`C_2^x`     :math:`i`   :math:`\sigma^{xy}` :math:`\sigma^{xz}` :math:`\sigma^{yz}`
   ============== ============== ================ =========== ========= ============== ================ ================= =========== =================== =================== ===================

.. table:: Symmetries available in |molcas|, including |molcas| keywords and symmetry elements.
   :name: tab:symmetry_input

   ============== ============ ============ ============ ========= ============== ================ ================= =========== =================== =================== ===================
   Group          Generators                             Elements
   -------------- -------------------------------------- -----------------------------------------------------------------------------------------------------------------------------------
   |zws|          :math:`g_1`  :math:`g_2`  :math:`g_3`  :math:`E` :math:`g_1`    :math:`g_2`      :math:`g_1g_2`    :math:`g_3` :math:`g_1g_3`      :math:`g_2g_3`      :math:`g_1g_2g_3`
   ============== ============ ============ ============ ========= ============== ================ ================= =========== =================== =================== ===================
   :math:`C_1`                                           :math:`E`
   :math:`C_2`    :kword:`xy`                            :math:`E` :math:`C_2`
   :math:`C_s`    :kword:`x`                             :math:`E` :math:`\sigma`
   :math:`C_i`    :kword:`xyz`                           :math:`E` :math:`i`
   :math:`C_{2v}` :kword:`xy`  :kword:`y`                :math:`E` :math:`C_2`    :math:`\sigma_v` :math:`\sigma_v'`
   :math:`C_{2h}` :kword:`xy`  :kword:`xyz`              :math:`E` :math:`C_2`    :math:`i`        :math:`\sigma_h`
   :math:`D_2`    :kword:`xy`  :kword:`xz`               :math:`E` :math:`C_2^z`  :math:`C_2^y`    :math:`C_2^x`
   :math:`D_{2h}` :kword:`xy`  :kword:`xz`  :kword:`xyz` :math:`E` :math:`C_2^z`  :math:`C_2^y`    :math:`C_2^x`     :math:`i`   :math:`\sigma^{xy}` :math:`\sigma^{xz}` :math:`\sigma^{yz}`
   ============== ============ ============ ============ ========= ============== ================ ================= =========== =================== =================== ===================
