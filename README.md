## Kink-Antikink Collisions in the $\phi^4$ Model

Simulating collisions between topological solitons, specifically kink-antikink collisions in the $\phi^4$ model, using a finite difference simulation scheme. Identified velocity thresholds for interesting observable behaviour consistent with published results.

## Key Results
- Critical velocity threshold. For initial velocity greater than $v = 0.2598$, the pair of waves always undergo a single inelastic collision before escaping to spatial infinity. We call this a one-bounce interaction. Some energy is lost in the collision, as is observed by the lower final velocity than initial.
- Lower trapping threshold. For $v \le 0.193$, the pair collide and are trapped in a permanant oscillatory state. This is due to the lack of sufficient energy to escape the mutual attraction between them.
- More complex behaviour in the intermediate window. Specifically, two-bounce windows observed.

## Result Diagrams

[Kink-antikink collision, v = 0.35](./images/phi-4.png)

*One-bounce interaction at v = 0.35 (t = 0, 10, 17, 19, 20, 22, 27, 45): the kink and
antikink approach, collide, and reflect apart to infinity.*

[One and two-bounce interaction comparison](./images/velocity_comparison.png)

*Displays the center of the domain ($x = 0$) across time steps for different starting velocities. Dips in $u(0, t)$ indicate when collisions occur.*

Can also be found in the report.

## Full Report
The complete write-up, including derivations and full results, is available
[here](./Kink_Antikink_Collisions_in_Phi4_Model.pdf).

## Code
The MATLAB code, which can be altered to reproduce the results in the report, is available [here](./kinkantikinkcollisionsphi4.m).
