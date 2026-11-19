# Matched four-method experiment

Twenty paired realizations use simulation seeds `seed+10000+r`, optimizer seeds `seed+20000+r`, the trusted held-site manifest, four multistarts, and the confirm-profile optimization budget.

Compared wall time includes fitting, any method-specific rank selection, and reconstruction. Data simulation, the shared numerical-correlation table, and shared basis-bank construction are excluded.

KL tolerance: 0.01. Initial screened rank: 121.

## Four-method summary

- numerical_inversion_reference: latent RMSE 0.112745 [0.111949, 0.113551]; observed RMSE 0.250753 [0.249819, 0.251772]; time 149.268 s [138.323, 160.230]; convergence 1.00.
- fixed_hsgp_M121: latent RMSE 0.113926 [0.113131, 0.114689]; observed RMSE 0.251304 [0.250340, 0.252357]; time 44.796 s [40.739, 49.171]; convergence 0.95.
- fixed_high_rank_hsgp: latent RMSE 0.112762 [0.111979, 0.113576]; observed RMSE 0.250758 [0.249786, 0.251762]; time 229.951 s [210.067, 248.950]; convergence 1.00.
- frequency_adaptive_hsgp: latent RMSE 0.114000 [0.113241, 0.114733]; observed RMSE 0.251286 [0.250342, 0.252362]; time 76.588 s [71.973, 81.648]; convergence 1.00.

## Paired RMSE differences

- fixed_hsgp_M121 minus numerical_inversion_reference: latent 0.001181 [0.001007, 0.001343]; observed 0.000551 [0.000480, 0.000617].
- fixed_high_rank_hsgp minus numerical_inversion_reference: latent 0.000016 [-0.000027, 0.000063]; observed 0.000005 [-0.000023, 0.000032].
- frequency_adaptive_hsgp minus numerical_inversion_reference: latent 0.001255 [0.001131, 0.001398]; observed 0.000534 [0.000440, 0.000637].
- frequency_adaptive_hsgp minus fixed_hsgp_M121: latent 0.000074 [-0.000058, 0.000216]; observed -0.000018 [-0.000126, 0.000091].
- frequency_adaptive_hsgp minus fixed_high_rank_hsgp: latent 0.001238 [0.001101, 0.001391]; observed 0.000528 [0.000430, 0.000636].

## Adaptive timing contrasts

- Versus numerical_inversion_reference: mean saved 72.680 s [63.491, 81.928]; mean saved fraction 0.4786 [0.4442, 0.5132].
- Versus fixed_hsgp_M121: mean saved -31.792 s [-33.635, -30.106]; mean saved fraction -0.7393 [-0.8112, -0.6711].
- Versus fixed_high_rank_hsgp: mean saved 153.364 s [132.837, 171.247]; mean saved fraction 0.6530 [0.6035, 0.6890].
