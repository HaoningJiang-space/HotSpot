/* Linux: cc -O2 -ffunction-sections -fdata-sections -DSUPERLU=0
 * scripts/test_hydraulic_edges.c -Wl,--gc-sections -lm -o /external/test-edges
 * Include implementation to exercise the exact assembly/flow helper. */
#include <assert.h>
#include <stdarg.h>
#include "../microchannel.c"

void fatal(char *message) { fprintf(stderr, "%s", message); exit(1); }

static int near(double a, double b) {
  return fabs(a-b) <= 1e-12 * fmax(fabs(a), fabs(b));
}

int main(void) {
  microchannel_config_t c = default_microchannel_config();
  double gx, gy;
  c.cell_width = 4e-4; c.cell_height = 2e-4;
  c.cell_thickness = 1e-4; c.coolant_visc = 1e-3;
  gx = edge_hydro_conductance(&c, 1);
  gy = edge_hydro_conductance(&c, 0);
  /* Rectangular-duct approximation: (1-.63 h/w) h^3 w/(12 mu L). */
  assert(near(gx, 2.854166666666667e-11));
  assert(near(gy, 1.404166666666667e-10));
  c.cell_width = 2e-4; c.cell_height = 4e-4;
  assert(near(edge_hydro_conductance(&c, 1), gy));
  assert(near(edge_hydro_conductance(&c, 0), gx));
  c.cell_width = c.cell_height = 2e-4;
  assert(near(edge_hydro_conductance(&c, 1), edge_hydro_conductance(&c, 0)));
  {
    int maprow[] = {0,1}; int *mapping[] = {maprow};
    double pressures[] = {2000,1000};
    c.mapping = mapping; c.b = pressures;
    assert(near(flow_rate(&c,0,0,0,1),
                1000*edge_hydro_conductance(&c,1)));
    assert(near(flow_rate(&c,0,0,0,1), -flow_rate(&c,0,1,0,0)));
  }
  {
    int cells[] = {INLET, FLUID, OUTLET}; int *types[] = {cells};
    int k;
    c = default_microchannel_config();
    c.num_rows = 1; c.num_columns = 3; c.cell_types = types;
    c.cell_width = 4e-4; c.cell_height = 2e-4;
    c.cell_thickness = 1e-4; c.coolant_visc = 1e-3;
    build_pressure_matrix(&c);
    assert(near(c.A[1][0], -gx));
    assert(near(c.A[1][2], -gx));
    assert(near(c.A[1][1], 2*gx));
    c.b[0] = 2000; c.b[1] = 1000; c.b[2] = 0;
    assert(near(flow_rate(&c,0,0,0,1), 1000*gx));
    assert(near(flow_rate(&c,0,0,0,1), flow_rate(&c,0,1,0,2)));
    for(k=0;k<3;k++) free(c.A[k]);
    free(c.A); free(c.b); free(c.mapping[0]); free(c.mapping);
  }
  puts("PASS: rectangular edges, rotation, square parity, reverse flow");
  return 0;
}
