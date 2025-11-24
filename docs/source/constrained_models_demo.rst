Constrained Models Demo
=======================

This demo generates omics‑constrained flux models using SCOOTI’s CFR/DFA framework, with a ready‑to‑run config and wrapper.

Direct run
----------

Run SCOOTI with the provided constrained demo config:

.. code-block:: bash

   bash scooti/run_flux.sh examples/constrained_demo/constrained_demo_config.json

Convenience wrapper
-------------------

Alternatively, use the wrapper script which prints output locations and a few example files:

.. code-block:: bash

   bash examples/constrained_demo/run_constrained.sh

Output
------

Results are saved under ``examples/constrained_demo/out/constrained_models/``.

Notes
-----

- Edit the JSON to adjust parameters such as ``rho``, ``kappa``, labels (``uplabel``/``dwlabel``), and ``data_dir`` for your significant gene lists.
- Switch to DFA by setting ``DFA_kappa > 0`` and providing metabolomics series.
- Use ``paraLen``/``jj`` and a small bash loop to scan parameters if needed.


Timing
------
- Typical run: 2–8 hours depending on solver (Gurobi/CPLEX), model size, and number of samples/conditions.
- Parameter scans (kappa/rho): add proportionally to above.

Notes
-----
- Ensure ``medium"`` is set (e.g., ``DMEMF12``) and consistent with demo configs.
- Demo config: ``examples/constrained_demo/constrained_demo_config.json`` (uses Shen2019 GEM).
