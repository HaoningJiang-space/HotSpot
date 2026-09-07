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
  {
    FILE *report = tmpfile();
    char buffer[4096], *field;
    double bank_g[] = {1e-9, 1e-9, 1e-9, 1e-9};
    double value;
    assert(report);
    c.cooling_branch_count = 4;
    c.solved_pump_pressure = 100.;
    c.pump_efficiency = .5;
    c.closed_branch_mask = 1;
    c.branch_flow[0] = 0.;
    c.branch_flow[1] = 5e-8;
    c.valve_resistance[1] = 1e9;
    write_hydraulic_report(report, &c, bank_g);
    rewind(report);
    assert(fgets(buffer, sizeof(buffer), report));
    field = strstr(buffer, "branch_0_valve_drop_pa=");
    assert(field && sscanf(field, "branch_0_valve_drop_pa=%lf", &value) == 1 && value == 100.);
    field = strstr(buffer, "branch_1_pressure_residual_pa=");
    assert(field && sscanf(field, "branch_1_pressure_residual_pa=%lf", &value) == 1 && fabs(value) < 1e-12);
    field = strstr(buffer, "branch_1_pump_power_w=");
    assert(field && sscanf(field, "branch_1_pump_power_w=%lf", &value) == 1 && near(value, 1e-5));
    fclose(report);
    c = default_microchannel_config();
  }
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
  {
    int refinement;
    double reference_flow = 0.;
    for(refinement = 1; refinement <= 3; refinement++) {
      int row, col;
      microchannel_config_t *p = malloc(sizeof(*p));
      *p = default_microchannel_config();
      strcpy(p->network_file, "synthetic-fixed-duct-test");
      p->num_rows = 8 * refinement;
      p->num_columns = 4 * refinement;
      p->cell_width = .01 / p->num_columns;
      p->cell_height = .004 / p->num_rows;
      p->cell_thickness = .0003;
      p->cooling_branch_count = 4;
      p->manifold_inlet_resistance = 1e9;
      p->pump_curve_resistance = 3e10;
      p->pump_efficiency = .7;
      p->pumping_pressure = 52000.;
      p->cell_types = calloc(p->num_rows, sizeof(int *));
      p->branch_ids = calloc(p->num_rows, sizeof(int *));
      for(row = 0; row < 4; row++) p->valve_resistance[row] = 5e10;
      for(row = 0; row < p->num_rows; row++) {
        p->cell_types[row] = calloc(p->num_columns, sizeof(int));
        p->branch_ids[row] = calloc(p->num_columns, sizeof(int));
        if((row / refinement) % 2) {
          for(col = 0; col < p->num_columns; col++) p->cell_types[row][col] = FLUID;
          p->cell_types[row][0] = OUTLET;
          p->cell_types[row][p->num_columns - 1] = INLET;
          p->branch_ids[row][p->num_columns - 1] = row / (2 * refinement);
        }
      }
      solve_physical_straight_ducts(p);
      if(refinement == 1) reference_flow = p->total_flow;
      assert(near(p->total_flow, reference_flow));
      assert(fabs(p->pump_curve_residual) < 1e-9);
      assert(p->hydraulic_conservation_error < 1e-18);
      for(row = refinement; row < 2 * refinement; row++) {
        double q = flow_rate(p, row, 1, row, 0);
        assert(near(q * refinement, p->branch_flow[0]));
        assert(near(q, -flow_rate(p, row, 0, row, 1)));
        if(row + 1 < 2 * refinement) assert(flow_rate(p, row, 1, row + 1, 1) == 0.);
      }
      free(p->physical_row_flow);
      p->physical_row_flow = NULL;
      p->closed_branch_mask = 1;
      solve_physical_straight_ducts(p);
      assert(p->branch_flow[0] == 0.);
      assert(p->branch_flow[1] > 0.);
      free(p->physical_row_flow);
      p->physical_row_flow = NULL;
      p->closed_branch_mask = 15;
      solve_physical_straight_ducts(p);
      assert(p->total_flow == 0. && p->pump_power == 0.);
      assert(p->solved_pump_pressure == p->pumping_pressure);
      free_microchannel(p);
    }
    puts("PASS: physical duct flow invariant under 1x/2x/3x refinement");
  }
  return 0;
}
