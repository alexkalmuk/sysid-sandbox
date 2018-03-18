### 1. MPC + LSCR for ARX systems

The problem is to construct robust control in
probabilistic sense for all possible values from some predefined
set. In that case it is a good way to shrink the set over time
in order to get better sets. The proposed algorithm is based
on the modified LSCR (Leave-out Sign-dominant Correlation
Regions) and applied to state space model used in MPC.

#### How to launch

The inputs for algorithm are located in mpc_lscr_arx/MPC_LSCR_Test_2D_Simple.m. You can read more [here](https://github.com/alexkalmuk/sysid-sandbox/wiki/Launching-MPC-LSRC-example)
