# Description
HotSpot is a pre-RTL thermal simulator intended for use early in the design process. HotSpot supports simulation of traditional 2D Integrated Circuits (2D ICs) and [3D ICs](https://en.wikipedia.org/wiki/Three-dimensional_integrated_circuit) as well as  [microfluidic cooling](https://en.wikipedia.org/wiki/Microfluidics). If this is your first time seeing HotSpot, check out our [Getting Started](https://github.com/uvahotspot/HotSpot/wiki/Getting-Started) page. To see HotSpot in action, check out our simulation examples in the `examples` directory!

## ThermoDSE cooling-control extension

The `hotspot7-cooling-krylov` branch adds a four-branch shared-pump model for
cooling-service experiments. Microchannel inlet tokens may carry a branch ID,
for example `2:0`. The cooling contract is controlled by
`cooling_branch_count`, `valve_resistance_0` through
`valve_resistance_3`, `manifold_inlet_resistance`,
`pump_curve_resistance`, and `pump_efficiency`. Branch flows remain hydraulic
state variables; they are not assigned independently by the optimizer.

Set `HOTSPOT_G7_HYDRAULIC_REPORT` to record branch flows, pressure, pump power,
mass-conservation error, and pump-curve residual. The reported flow unit is
`m^3/s`, so pump power is `total_flow * pump_pressure / pump_efficiency`.

The extension also exposes a HotSpot-native matrix-free thermal action.
`HOTSPOT_G7_SOLVER=auto` admits PCG only for the no-flow SPD regime and retains
SuperLU for a flowing, nonsymmetric coupled operator. `partitioned` enables an
experimental full-system FGMRES path with a solid/fluid block preconditioner;
it is not a production default. `HOTSPOT_G7_OPERATOR_PREFIX` writes assembled
and matrix-free audit data for matched operator-validation runs.
