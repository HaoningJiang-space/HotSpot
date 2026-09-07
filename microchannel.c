#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "microchannel.h"

#if SUPERLU > 0
#include "slu_ddefs.h"
#endif

#define MAX_LINE_SIZE 4096
#define DEBUG 0

static int uses_shared_pump(const microchannel_config_t *config)
{
  return config->cooling_branch_count > 0;
}

static int pressure_node_count(const microchannel_config_t *config)
{
  if(uses_shared_pump(config))
    return config->n_fluid_cells + config->cooling_branch_count + 1;
  return config->n_fluid_cells +
    (config->pump_internal_res != 0.0 ? EXTRA_PRESSURE_NODES : 0);
}

static void record_hydraulic_operating_point(microchannel_config_t *config);

// default microchannel configuration parameters
microchannel_config_t default_microchannel_config(void)
{
  microchannel_config_t config;
  int branch;

  config.cell_width        = 100e-6;
  config.cell_height       = 100e-6;
  config.cell_thickness    = 100e-6;
  config.pumping_pressure  = 5000;
  config.pump_internal_res = 0;            // ideal pump
  config.cooling_branch_count = 0;
  for(branch = 0; branch < MAX_COOLING_BRANCHES; branch++) {
    config.valve_resistance[branch] = 0.0;
    config.branch_flow[branch] = 0.0;
  }
  config.pump_curve_resistance = 0.0;
  config.manifold_inlet_resistance = 0.0;
  config.pump_efficiency       = 1.0;
  config.inlet_temperature = 300;
  config.coolant_capac     = 4172638;      // water
  config.coolant_res       = 1.647717911;  // water
  config.coolant_visc      = 0.000889;     // water
  config.wall_capac        = 1635660;      // silicon
  config.wall_res          = 0.0076923077; // silicon
  config.htc               = 27132;
  config.num_rows          = -1;
  config.num_columns       = -1;
  config.n_fluid_cells     = -1;
  config.cell_types        = NULL;
  config.branch_ids        = NULL;
  config.mapping           = NULL;
  config.A                 = NULL;
  config.b                 = NULL;
  config.nnz               = 0;
  config.total_flow        = 0.0;
  config.solved_pump_pressure = 0.0;
  config.pump_power        = 0.0;
  config.hydraulic_conservation_error = 0.0;
  config.pump_curve_residual = 0.0;

  return config;
}

/*
 * parse a table of name-value string pairs and add the configuration
 * parameters to 'config'
 */
void microchannel_config_add_from_strs(microchannel_config_t *config, materials_list_t *materials_list, str_pair *table, int size)
{
  int idx, branch;
  char option[STR_SIZE];

  if ((idx = get_str_index(table, size, "pumping_pressure")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->pumping_pressure) != 1)
      fatal("invalid format for configuration  parameter pumping_pressure\n");
  if ((idx = get_str_index(table, size, "pump_internal_res")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->pump_internal_res) != 1)
      fatal("invalid format for configuration  parameter pumping internal resistance\n");
  if ((idx = get_str_index(table, size, "cooling_branch_count")) >= 0)
    if(sscanf(table[idx].value, "%d", &config->cooling_branch_count) != 1)
      fatal("invalid format for configuration parameter cooling_branch_count\n");
  for(branch = 0; branch < MAX_COOLING_BRANCHES; branch++) {
    sprintf(option, "valve_resistance_%d", branch);
    if ((idx = get_str_index(table, size, option)) >= 0)
      if(sscanf(table[idx].value, "%lf", &config->valve_resistance[branch]) != 1)
        fatal("invalid format for branch valve resistance\n");
  }
  if ((idx = get_str_index(table, size, "pump_curve_resistance")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->pump_curve_resistance) != 1)
      fatal("invalid format for configuration parameter pump_curve_resistance\n");
  if ((idx = get_str_index(table, size, "manifold_inlet_resistance")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->manifold_inlet_resistance) != 1)
      fatal("invalid format for configuration parameter manifold_inlet_resistance\n");
  if ((idx = get_str_index(table, size, "pump_efficiency")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->pump_efficiency) != 1)
      fatal("invalid format for configuration parameter pump_efficiency\n");
  if ((idx = get_str_index(table, size, "inlet_temperature")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->inlet_temperature) != 1)
      fatal("invalid format for configuration  parameter inlet_temperature\n");
  if ((idx = get_str_index(table, size, "coolant_capac")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->coolant_capac) != 1)
      fatal("invalid format for configuration  parameter coolant_capacity\n");
  if ((idx = get_str_index(table, size, "coolant_res")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->coolant_res) != 1)
      fatal("invalid format for configuration  parameter coolant resistivity\n");
  if ((idx = get_str_index(table, size, "coolant_visc")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->coolant_visc) != 1)
      fatal("invalid format for configuration  parameter coolant dynamic viscosity\n");
  if ((idx = get_str_index(table, size, "wall_capac")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->wall_capac) != 1)
      fatal("invalid format for configuration  parameter wall_capacity\n");
  if ((idx = get_str_index(table, size, "wall_res")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->wall_res) != 1)
      fatal("invalid format for configuration  parameter wall_resistivity\n");
  if ((idx = get_str_index(table, size, "htc")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->htc) != 1)
      fatal("invalid format for configuration  parameter heat transfer coefficient\n");
  if ((idx = get_str_index(table, size, "network_file")) >= 0)
    if(sscanf(table[idx].value, "%s", config->network_file) != 1)
      fatal("invalid format for configuration  parameter network_file\n");

  if ((idx = get_str_index(table, size, "wall_material")) >= 0) {
    char material_name[STR_SIZE];
    if(sscanf(table[idx].value, "%s", material_name) != 1)
      fatal("invalid format for configuration parameter wall_material\n");

    config->wall_res = 1.0 / get_material_thermal_conductivity(materials_list, material_name);
    config->wall_capac = get_material_volumetric_heat_capacity(materials_list, material_name);

    if(config->wall_res < 0 || config->wall_capac < 0)
      fatal("material name specified in configuration parameter wall_material not found\n");
  }

  if ((idx = get_str_index(table, size, "coolant_material")) >= 0) {
    char material_name[STR_SIZE];
    if(sscanf(table[idx].value, "%s", material_name) != 1)
      fatal("invalid format for configuration parameter coolant_material\n");

    config->coolant_res = 1.0 / get_material_thermal_conductivity(materials_list, material_name);
    config->coolant_capac = get_material_volumetric_heat_capacity(materials_list, material_name);
    config->coolant_visc = get_material_dynamic_viscosity(materials_list, material_name);

    if(config->coolant_res < 0 || config->coolant_capac < 0 || config->coolant_visc < 0)
      fatal("material name specified in configuration parameter coolant_material not found\n");
  }

  if(config->cooling_branch_count < 0 ||
     config->cooling_branch_count > MAX_COOLING_BRANCHES)
    fatal("cooling_branch_count must be between zero and four\n");
  if(config->pump_curve_resistance < 0.0)
    fatal("pump_curve_resistance must be nonnegative\n");
  if(config->pump_efficiency <= 0.0 || config->pump_efficiency > 1.0)
    fatal("pump_efficiency must be in (0, 1]\n");
  for(branch = 0; branch < config->cooling_branch_count; branch++)
    if(config->valve_resistance[branch] <= 0.0)
      fatal("active branch valve resistance must be positive\n");
  if(uses_shared_pump(config) && config->manifold_inlet_resistance <= 0.0)
    fatal("shared-pump manifold inlet resistance must be positive\n");
  if(uses_shared_pump(config) && config->pump_internal_res != 0.0)
    fatal("shared-pump controls cannot be combined with legacy pump_internal_res\n");
}

/*
 * convert config into a table of name-value pairs. returns the no.
 * of parameters converted
 */
int microchannel_config_to_strs(microchannel_config_t *config, str_pair *table, int max_entries)
{
  if (max_entries < 25)
    fatal("not enough entries in table\n");

  sprintf(table[0].name, "cell_width");
  sprintf(table[1].name, "cell_height");
  sprintf(table[2].name, "cell_thickness");
  sprintf(table[3].name, "pumping_pressure");
  sprintf(table[4].name, "pump_internal_res");
  sprintf(table[5].name, "inlet_temperature");
  sprintf(table[6].name, "coolant_capac");
  sprintf(table[7].name, "coolant_res");
  sprintf(table[8].name, "coolant_visc");
  sprintf(table[9].name, "wall_capac");
  sprintf(table[10].name, "wall_res");
  sprintf(table[11].name, "htc");
  sprintf(table[12].name, "network_file");
  sprintf(table[13].name, "floorplan_file");
  sprintf(table[14].name, "num_rows");
  sprintf(table[15].name, "num_columns");
  sprintf(table[16].name, "n_fluid_cells");
  sprintf(table[17].name, "cooling_branch_count");
  sprintf(table[18].name, "valve_resistance_0");
  sprintf(table[19].name, "valve_resistance_1");
  sprintf(table[20].name, "valve_resistance_2");
  sprintf(table[21].name, "valve_resistance_3");
  sprintf(table[22].name, "pump_curve_resistance");
  sprintf(table[23].name, "pump_efficiency");
  sprintf(table[24].name, "manifold_inlet_resistance");

  sprintf(table[0].value, "%e", config->cell_width);
  sprintf(table[1].value, "%e", config->cell_height);
  sprintf(table[2].value, "%e", config->cell_thickness);
  sprintf(table[3].value, "%e", config->pumping_pressure);
  sprintf(table[4].value, "%e", config->pump_internal_res);
  sprintf(table[5].value, "%e", config->inlet_temperature);
  sprintf(table[6].value, "%e", config->coolant_capac);
  sprintf(table[7].value, "%e", config->coolant_res);
  sprintf(table[8].value, "%e", config->coolant_visc);
  sprintf(table[9].value, "%e", config->wall_capac);
  sprintf(table[10].value, "%e", config->wall_res);
  sprintf(table[11].value, "%e", config->htc);
  sprintf(table[12].value, "%s", config->network_file);
  sprintf(table[13].value, "%s", config->floorplan_file);
  sprintf(table[14].value, "%d", config->num_rows);
  sprintf(table[15].value, "%d", config->num_columns);
  sprintf(table[16].value, "%d", config->n_fluid_cells);
  sprintf(table[17].value, "%d", config->cooling_branch_count);
  sprintf(table[18].value, "%e", config->valve_resistance[0]);
  sprintf(table[19].value, "%e", config->valve_resistance[1]);
  sprintf(table[20].value, "%e", config->valve_resistance[2]);
  sprintf(table[21].value, "%e", config->valve_resistance[3]);
  sprintf(table[22].value, "%e", config->pump_curve_resistance);
  sprintf(table[23].value, "%e", config->pump_efficiency);
  sprintf(table[24].value, "%e", config->manifold_inlet_resistance);

  return 25;
}

void solve_pressure_circuit(microchannel_config_t *config) {
#if SUPERLU > 0
  SuperMatrix A, L, U, B;
  double *a, *rhs;
  int *asub, *xa;
  int *perm_r;
  int *perm_c;
  int nnz, nrhs, info, i, m, n, perc_spec;
  int j;
  superlu_options_t options;
  SuperLUStat_t stat;

  int node_count = pressure_node_count(config);
  m = n = node_count;
  nnz = config->nnz;
  if( !(a = doubleMalloc(nnz)) ) ABORT("malloc failed");
  if( !(asub = intMalloc(nnz)) ) ABORT("malloc failed");
  if( !(xa = intMalloc(n + 1)) ) ABORT("malloc failed");

  int v = 0;
  int x = 1;
  xa[0] = 0;
  for(i = 0; i < m; i++) {
    for(j = 0; j < n; j++) {
      if(config->A[i][j] != 0) {
        a[v] = config->A[i][j];
        asub[v++] = j;
      }
    }
    xa[x++] = v;
  }
  if(v != nnz)
    fatal("Pressure-network nonzero count does not match assembled matrix\n");

  if(DEBUG) {
    fprintf(stderr, "\n");
    for(i = 0; i < nnz; i++)
      fprintf(stderr, "a[%d] = %.15lf\n", i, a[i]);

    for(i = 0; i < nnz; i++)
      fprintf(stderr, "asub[%d] = %d\n", i, asub[i]);

    for(i = 0; i < n+1; i++)
      fprintf(stderr, "xa[%d] = %d\n", i, xa[i]);
  }

  dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NR, SLU_D, SLU_GE);

  nrhs = 1;
  if( !(rhs = doubleMalloc(m * nrhs)) ) ABORT("malloc failed");

  for(i = 0; i < m; i++) {
    rhs[i] = config->b[i];
  }

  dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);

  if( !(perm_r = intMalloc(m)) ) ABORT("malloc failed");
  if( !(perm_c = intMalloc(n)) ) ABORT("malloc failed");

  set_default_options(&options);
  options.ColPerm = NATURAL;

  StatInit(&stat);

  dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
  if(info != 0)
    fatal("Unable to solve microchannel pressure network\n");

  if(DEBUG) {
    //dPrint_CompCol_Matrix("A", &A);
    //dPrint_Dense_Matrix("B", &B);
  }

  DNformat *Astore = (DNformat *) B.Store;
  double *dp = (double *) Astore->nzval;

  for(i = 0; i < n; i++) {
      config->b[i] = dp[i];
      if(DEBUG)
        fprintf(stderr, "config->b[%d] = %e\n", i, config->b[i]);
    }

  record_hydraulic_operating_point(config);

  SUPERLU_FREE(rhs);
  SUPERLU_FREE(perm_r);
  SUPERLU_FREE(perm_c);
  Destroy_CompCol_Matrix(&A);
  Destroy_SuperMatrix_Store(&B);
  Destroy_SuperNode_Matrix(&L);
  Destroy_CompCol_Matrix(&U);
  StatFree(&stat);
#else

  int node_count = pressure_node_count(config);
  gaussj(config->A, node_count, config->b);
  record_hydraulic_operating_point(config);

  if(DEBUG) {
    int i;
    for(i = 0; i < node_count; i++)
      fprintf(stderr, "config->b[%d] = %e\n", i, config->b[i]);
  }

#endif
}

// Parse config's CSV file to build internal array of microchannel network
void microchannel_build_network(microchannel_config_t *config) {
  char line[MAX_LINE_SIZE], str[STR_SIZE];
  char *cell;
  int cell_type, branch_id, i = 0, j = 0;
  FILE *fp = fopen(config->network_file, "r");

  if(DEBUG)
    fprintf(stderr, "network_file = %s\n", config->network_file);

  if(!fp) {
    strncpy(str, "Unable to open ", STR_SIZE);
    strncat(str, config->network_file, STR_SIZE - strlen(str));
    strncat(str, "\n", STR_SIZE - strlen(str));
    fatal(str);
  }

  int nr = config->num_rows;
  int nc = config->num_columns;

  if(DEBUG)
    fprintf(stderr, "num_rows: %d, num_cols: %d\n", nr, nc);

  config->cell_types = calloc(nr, sizeof(int *));
  config->branch_ids = calloc(nr, sizeof(int *));

  if (config->cell_types == NULL || config->branch_ids == NULL)
    fatal("Unable to allocate microchannel cell metadata\n");

  for(i = 0; i < nr; i++) {
    config->cell_types[i] = calloc(nc, sizeof(int));
    config->branch_ids[i] = malloc(nc * sizeof(int));

    if (config->cell_types[i] == NULL || config->branch_ids[i] == NULL)
      fatal("Unable to allocate microchannel row metadata\n");
    for(j = 0; j < nc; j++)
      config->branch_ids[i][j] = -1;
  }

  // parse network file to build cell_types array
  i = j = 0;
  fgets(line, MAX_LINE_SIZE, fp);
  while(!feof(fp)) {
    if(i >= nr) {
      fclose(fp);
      sprintf(str, "Microchannel has more rows than num_rows(%d)\n", nr);
      fatal(str);
    }
    cell = strtok(line, " ,\n");
    while(cell != NULL) {
      if(j >= nc) {
        fclose(fp);
        sprintf(str, "Microchannel row %d has more cells than num_columns(%d)\n", i, nc);
        fatal(str);
      }
      branch_id = -1;
      if(sscanf(cell, "%d:%d", &cell_type, &branch_id) < 1) {
        fclose(fp);
        fatal("Invalid microchannel cell token\n");
      }
      if(branch_id < -1 || branch_id >= config->cooling_branch_count) {
        fclose(fp);
        fatal("Microchannel branch identifier is outside cooling_branch_count\n");
      }
      config->cell_types[i][j] = cell_type;
      config->branch_ids[i][j] = branch_id;
      if(DEBUG)
        fprintf(stderr, "Adding cell %d, %d; type = %d\n", i, j, config->cell_types[i][j]);
      j++;


      cell = strtok(NULL, " ,\n");
    }

    if(j != nc) {
      fclose(fp);
      sprintf(str, "Microchannel row %d has %d cells, expected %d\n", i, j, nc);
      fatal(str);
    }

    i++; j = 0;

    fgets(line, MAX_LINE_SIZE, fp);
  }

  fclose(fp);
  if(i != nr)
    fatal("Microchannel has fewer rows than num_rows\n");

  if(uses_shared_pump(config)) {
    for(i = 0; i < nr; i++)
      for(j = 0; j < nc; j++)
        if(IS_INLET_CELL(config, i, j) && config->branch_ids[i][j] < 0)
          fatal("Every inlet requires a branch identifier in shared-pump mode\n");
  }

  // Create floorplan file for microchannel
  strcpy(config->floorplan_file, config->network_file);
  char *ext = strstr(config->floorplan_file, NETWORK_EXTENSION);
  strcpy(ext, FLOORPLAN_EXTENSION);
  fp = fopen(config->floorplan_file, "w");

  fprintf(fp, "# Name\tw\th\tx\ty\tc_v\tp\n");
  for(i = 0; i < nr; i++) {
    for(j = 0; j < nc; j++) {
      if(IS_FLUID_CELL(config, i, j)) {
        // For floorplan files, x = y = 0 is the bottom left corner, but for our parsing i = j = 0 is
        // the top left cell, so we have to account for that
        fprintf(fp, "Cell_%d_%d\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n", i, j, config->cell_width, config->cell_height,
                j*config->cell_width, (config->num_rows - i - 1)*config->cell_height, config->coolant_capac,
                config->coolant_res);
      }
      else {
        fprintf(fp, "Cell_%d_%d\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n", i, j, config->cell_width, config->cell_height,
                j*config->cell_width, (config->num_rows - i - 1)*config->cell_height, config->wall_capac,
                config->wall_res);
      }
    }
  }

  fclose(fp);

  printf("Creating pressure circuit...\n");
  build_pressure_matrix(config);
  printf("Solving pressure circuit...\n");
  solve_pressure_circuit(config);
}

static double edge_hydro_conductance(microchannel_config_t *config,
                                    int horizontal) {
  double h = config->cell_thickness;
  /* Cross-section is transverse to the edge; length is along the edge.
   * This cell-network model does not represent transverse refinement of a
   * single physical duct. Such refinement requires a separate duct model. */
  double w = horizontal ? config->cell_height : config->cell_width;
  double L = horizontal ? config->cell_width : config->cell_height;
  double viscosity = config->coolant_visc;
  double ret_val;

  if(!isfinite(h) || !isfinite(w) || !isfinite(L) || !isfinite(viscosity) ||
     h <= 0 || w <= 0 || L <= 0 || viscosity <= 0)
    fatal("Invalid hydraulic edge dimensions or viscosity\n");

  if(h == w)
    ret_val = (0.42229 * pow(h, 4)) / (12 * viscosity * L);
  else if(h > w)
    ret_val = ((1 - 0.63*(w / h)) * (pow(w, 3)) * (h)) / (12 * viscosity * L);
  else
    ret_val = ((1 - 0.63*(h / w)) * (pow(h, 3)) * (w)) / (12 * viscosity * L);

  return ret_val;
}

void build_pressure_matrix(microchannel_config_t *config) {
  int i, j, branch;
  int nr = config->num_rows;
  int nc = config->num_columns;
  int **mapping;
  int n = 0;
  config->nnz = 0;

  mapping = calloc(nr, sizeof(int *));

  if(mapping == NULL)
    fatal("Unable to allocate pressure circuit mapping\n");

  for(i = 0; i < nr; i++) {
    mapping[i] = calloc(nc, sizeof(int));

    if(mapping[i] == NULL) {
      fatal("Unable to allocate pressure circuit mapping\n");
    }
  }

  // Assign unique number to each fluid cell
  for(i = 0; i < nr; i++) {
    for(j = 0; j < nc; j++) {
      if(IS_FLUID_CELL(config, i, j)) {
        mapping[i][j] = n++;
      }
      else
        mapping[i][j] = -1;
    }
  }
  config->n_fluid_cells = n;
  config->mapping = mapping;

  // If we're modeling a non-ideal pump, include one extra pressure node.
  // Keep the dimension local to this config: multiple microchannel layers may
  // legitimately use different pump models.
  int node_count = pressure_node_count(config);
  int shared_pump = uses_shared_pump(config);
  int branch_node_base = config->n_fluid_cells;
  int pump_node = config->n_fluid_cells + config->cooling_branch_count;
  int branch_inlets[MAX_COOLING_BRANCHES] = {0, 0, 0, 0};

  config->A = calloc(node_count, sizeof(double *));

  if(!config->A)
    fatal("Unable to allocate matrix A for pressure circuit\n");
  for(i = 0; i < node_count; i++) {
    config->A[i] = calloc(node_count, sizeof(double));

    if(!config->A[i])
      fatal("Unable to allocate matrix A for pressure circuit\n");
  }

  config->b = calloc(node_count, sizeof(double));

  if(!config->b)
    fatal("Unable to allocate matrix b for pressure circuit\n");


  if(DEBUG) {
    fprintf(stderr, "Mapping number: %d\n", n);
    for(i = 0; i < nr; i++) {
      for(j = 0; j < nc; j++) {
        fprintf(stderr, "mapping[%d][%d] = %d\n", i, j, mapping[i][j]);
      }
    }
  }

  // Iterate through all cells
  double diagonal_val = 0;
  double hydro_conductance;
  for(i = 0; i < nr; i++) {
    for(j = 0; j < nc; j++) {
      if(config->cell_types[i][j] == FLUID ||
         (config->cell_types[i][j] == INLET &&
          (shared_pump || config->pump_internal_res != 0))) {
        // northern cell
          if(i > 0 && IS_FLUID_CELL(config, i-1, j)) {
            hydro_conductance = -edge_hydro_conductance(config, FALSE);
            if(DEBUG)
              fprintf(stderr, "[%d, %d]: Northern cell. Setting A[%d][%d] = %.15lf\n", i, j, mapping[i][j], mapping[i-1][j], hydro_conductance);

            config->A[mapping[i][j]][mapping[i-1][j]] = hydro_conductance;
            diagonal_val += hydro_conductance;
            config->nnz++;
          }

        // southern cell
          if(i < nr - 1 && IS_FLUID_CELL(config, i+1, j)) {
            hydro_conductance = -edge_hydro_conductance(config, FALSE);
            if(DEBUG)
              fprintf(stderr, "[%d, %d]: Southern cell. Setting A[%d][%d] = %.15lf\n", i, j, mapping[i][j], mapping[i+1][j], hydro_conductance);

            config->A[mapping[i][j]][mapping[i+1][j]] = hydro_conductance;
            diagonal_val += hydro_conductance;
            config->nnz++;
          }

        // western cell
          if(j > 0 && IS_FLUID_CELL(config, i, j-1)) {
            hydro_conductance = -edge_hydro_conductance(config, TRUE);
            if(DEBUG)
              fprintf(stderr, "[%d, %d]: Western Cell. Setting A[%d][%d] = %.15lf\n", i, j, mapping[i][j], mapping[i][j-1], hydro_conductance);

            config->A[mapping[i][j]][mapping[i][j-1]] = hydro_conductance;
            diagonal_val += hydro_conductance;
            config->nnz++;
          }

        // eastern cell
          if(j < nc - 1 && IS_FLUID_CELL(config, i, j+1)) {
            hydro_conductance = -edge_hydro_conductance(config, TRUE);
            if(DEBUG)
              fprintf(stderr, "[%d, %d]: Eastern Cell. Setting A[%d][%d] = %.15lf\n", i, j, mapping[i][j], mapping[i][j+1], hydro_conductance);

            config->A[mapping[i][j]][mapping[i][j+1]] = hydro_conductance;
            diagonal_val += hydro_conductance;
            config->nnz++;
          }

        // diagonal
        if(DEBUG)
          fprintf(stderr, "[%d, %d]: Diagonal. Setting A[%d][%d] = %.15lf\n", i, j, mapping[i][j], mapping[i][j], -diagonal_val);

        config->A[mapping[i][j]][mapping[i][j]] = -diagonal_val;
        config->nnz++;
        diagonal_val = 0;
      }

      if(IS_INLET_CELL(config, i, j)) {
        if(shared_pump) {
          double manifold_conductance =
            1.0 / config->manifold_inlet_resistance;
          int branch_node;
          branch = config->branch_ids[i][j];
          branch_node = branch_node_base + branch;
          config->A[mapping[i][j]][branch_node] = -manifold_conductance;
          config->A[mapping[i][j]][mapping[i][j]] += manifold_conductance;
          config->A[branch_node][mapping[i][j]] = -manifold_conductance;
          config->A[branch_node][branch_node] += manifold_conductance;
          branch_inlets[branch]++;
          config->nnz += 2;
        }
        // Non-ideal pump
        else if(config->pump_internal_res != 0) {
          if(DEBUG) {
            fprintf(stderr, "[%d, %d]: Inlet. Setting A[%d][%d] = %e\n", i, j, mapping[i][j], config->n_fluid_cells, -1.0 / config->pump_internal_res);
            fprintf(stderr, "[%d, %d]: Inlet. Setting A[%d][%d] = %e\n", i, j, mapping[i][j], mapping[i][j], 1.0 / config->pump_internal_res);
          }

          // Inlet cells are connected to the pump through the pump's internal
          // resistance
          config->A[mapping[i][j]][config->n_fluid_cells] = -1.0 / config->pump_internal_res;
          config->A[mapping[i][j]][mapping[i][j]] += 1.0 / config->pump_internal_res;
          config->nnz++; // diagonal val already exists, so we're only adding one nonzero value
        }

        // Ideal pump
        else {
          if(DEBUG) {
            fprintf(stderr, "[%d, %d]: Inlet. Setting A[%d][%d] = %e\n", i, j, mapping[i][j], mapping[i][j], 1.0);
            fprintf(stderr, "[%d, %d]: Inlet. Setting b[%d] = %e\n", i, j, mapping[i][j], config->pumping_pressure);
          }
          config->A[mapping[i][j]][mapping[i][j]] = 1.0;
          config->b[mapping[i][j]] = config->pumping_pressure;
          config->nnz++;
        }
      }

      else if(IS_OUTLET_CELL(config, i, j)) {
        if(DEBUG)
          fprintf(stderr, "[%d, %d]: OUTLET. Setting A[%d][%d] = 1\n", i, j, mapping[i][j], mapping[i][j]);

        config->A[mapping[i][j]][mapping[i][j]] = 1;
        config->nnz++;
      }
    }
  }

  if(shared_pump) {
    double pump_conductance_sum = 0.0;
    for(branch = 0; branch < config->cooling_branch_count; branch++) {
      int branch_node = branch_node_base + branch;
      double valve_conductance = 1.0 / config->valve_resistance[branch];
      if(branch_inlets[branch] == 0)
        fatal("Every active cooling branch requires at least one inlet\n");
      config->A[branch_node][branch_node] += valve_conductance;
      config->A[branch_node][pump_node] = -valve_conductance;
      pump_conductance_sum += valve_conductance;
      config->nnz += 2;
    }
    if(config->pump_curve_resistance > 0.0) {
      double pump_ground_conductance =
        1.0 / config->pump_curve_resistance;
      config->A[pump_node][pump_node] =
        pump_conductance_sum + pump_ground_conductance;
      config->b[pump_node] =
        pump_ground_conductance * config->pumping_pressure;
      config->nnz++;
      for(branch = 0; branch < config->cooling_branch_count; branch++) {
        config->A[pump_node][branch_node_base + branch] =
          -1.0 / config->valve_resistance[branch];
        config->nnz++;
      }
    }
    else {
      config->A[pump_node][pump_node] = 1.0;
      config->b[pump_node] = config->pumping_pressure;
      config->nnz++;
    }
  }
  else if(config->pump_internal_res != 0) {
    if(DEBUG) {
      fprintf(stderr, "Pump Node. Setting A[%d][%d] = %e\n", config->n_fluid_cells, config->n_fluid_cells, 1.0);
      fprintf(stderr, "Pump Node. Setting b[%d] = %e\n", config->n_fluid_cells, config->pumping_pressure);
    }

    // Handle extra node representing pump
    config->A[config->n_fluid_cells][config->n_fluid_cells] = 1;
    config->b[config->n_fluid_cells] = config->pumping_pressure;
    config->nnz++;
  }

  if(DEBUG) {
    fprintf(stderr, "Nonzero values (%d total):\n", config->nnz);
     for(i = 0; i < node_count; i++) {
      for(j = 0; j < node_count; j++) {
        if(config->A[i][j] != 0)
          fprintf(stderr, "A[%d][%d] = %e\n", i, j, config->A[i][j]);
      }
    }

    if(DEBUG)
      for(i = 0; i < node_count; i++)
        fprintf(stderr, "b[%d] = %e\n", i, config->b[i]);

  }
}

double flow_rate(microchannel_config_t * config, int cell1_i, int cell1_j, int cell2_i, int cell2_j) {
  double *pressure = config->b;
  int **mapping = config->mapping;
  if(abs(cell1_i - cell2_i) + abs(cell1_j - cell2_j) != 1)
    fatal("Hydraulic flow requires adjacent cells\n");
  return (pressure[mapping[cell1_i][cell1_j]] - pressure[mapping[cell2_i][cell2_j]]) *
    edge_hydro_conductance(config, cell1_i == cell2_i);
}

static double inlet_network_flow(microchannel_config_t *config, int row,
                                 int column)
{
  double flow = 0.0;
  if(row > 0 && IS_FLUID_CELL(config, row - 1, column))
    flow += flow_rate(config, row, column, row - 1, column);
  if(row < config->num_rows - 1 &&
     IS_FLUID_CELL(config, row + 1, column))
    flow += flow_rate(config, row, column, row + 1, column);
  if(column > 0 && IS_FLUID_CELL(config, row, column - 1))
    flow += flow_rate(config, row, column, row, column - 1);
  if(column < config->num_columns - 1 &&
     IS_FLUID_CELL(config, row, column + 1))
    flow += flow_rate(config, row, column, row, column + 1);
  return flow;
}

static void write_hydraulic_report(FILE *stream,
                                   microchannel_config_t *config)
{
  int branch;
  fprintf(stream,
          "HotSpot 7 hydraulic state: network=%s branches=%d "
          "pump_pressure_pa=%.17g total_flow_m3_s=%.17g "
          "pump_power_w=%.17g conservation_error_m3_s=%.17g "
          "pump_curve_residual_pa=%.17g",
          config->network_file, config->cooling_branch_count,
          config->solved_pump_pressure, config->total_flow,
          config->pump_power, config->hydraulic_conservation_error,
          config->pump_curve_residual);
  for(branch = 0; branch < config->cooling_branch_count; branch++)
    fprintf(stream, " branch_%d_flow_m3_s=%.17g", branch,
            config->branch_flow[branch]);
  fprintf(stream, "\n");
}

static void record_hydraulic_operating_point(microchannel_config_t *config)
{
  int row, column, branch;
  double network_flow = 0.0;
  double valve_flow = 0.0;
  double branch_network_flow[MAX_COOLING_BRANCHES] = {0.0, 0.0, 0.0, 0.0};
  const char *report_path;
  FILE *report;

  for(branch = 0; branch < MAX_COOLING_BRANCHES; branch++)
    config->branch_flow[branch] = 0.0;

  if(uses_shared_pump(config))
    config->solved_pump_pressure =
      config->b[config->n_fluid_cells + config->cooling_branch_count];
  else if(config->pump_internal_res != 0.0)
    config->solved_pump_pressure = config->b[config->n_fluid_cells];
  else
    config->solved_pump_pressure = config->pumping_pressure;

  for(row = 0; row < config->num_rows; row++) {
    for(column = 0; column < config->num_columns; column++) {
      if(IS_INLET_CELL(config, row, column)) {
        double inlet_flow = inlet_network_flow(config, row, column);
        network_flow += inlet_flow;
        if(uses_shared_pump(config)) {
          branch = config->branch_ids[row][column];
          branch_network_flow[branch] += inlet_flow;
        }
      }
    }
  }

  if(uses_shared_pump(config)) {
    for(branch = 0; branch < config->cooling_branch_count; branch++) {
      int branch_node = config->n_fluid_cells + branch;
      config->branch_flow[branch] =
        (config->solved_pump_pressure - config->b[branch_node]) /
        config->valve_resistance[branch];
      valve_flow += config->branch_flow[branch];
    }
  }

  config->total_flow = uses_shared_pump(config) ? valve_flow : network_flow;
  config->hydraulic_conservation_error = 0.0;
  if(uses_shared_pump(config))
    for(branch = 0; branch < config->cooling_branch_count; branch++)
      config->hydraulic_conservation_error +=
        fabs(config->branch_flow[branch] - branch_network_flow[branch]);
  config->pump_curve_residual = uses_shared_pump(config) ?
    config->solved_pump_pressure +
      config->pump_curve_resistance * config->total_flow -
      config->pumping_pressure : 0.0;
  config->pump_power = config->total_flow * config->solved_pump_pressure /
    config->pump_efficiency;

  write_hydraulic_report(stdout, config);
  report_path = getenv("HOTSPOT_G7_HYDRAULIC_REPORT");
  if(report_path && report_path[0]) {
    report = fopen(report_path, "a");
    if(!report)
      fatal("Unable to open HotSpot 7 hydraulic report\n");
    write_hydraulic_report(report, config);
    fclose(report);
  }
}

// Copy user-defined parameters from one microchannel config to another
void copy_microchannel(microchannel_config_t *dst, microchannel_config_t *src) {
  int branch;
  dst->pumping_pressure  = src->pumping_pressure;
  dst->pump_internal_res = src->pump_internal_res;
  dst->cooling_branch_count = src->cooling_branch_count;
  for(branch = 0; branch < MAX_COOLING_BRANCHES; branch++)
    dst->valve_resistance[branch] = src->valve_resistance[branch];
  dst->pump_curve_resistance = src->pump_curve_resistance;
  dst->manifold_inlet_resistance = src->manifold_inlet_resistance;
  dst->pump_efficiency = src->pump_efficiency;
  dst->inlet_temperature = src->inlet_temperature;
  dst->coolant_capac     = src->coolant_capac;
  dst->coolant_res       = src->coolant_res;
  dst->coolant_visc      = src->coolant_visc;
  dst->wall_capac        = src->wall_capac;
  dst->wall_res          = src->wall_res;
  dst->htc               = src->htc;
}

void free_microchannel(microchannel_config_t *config) {
  int i;
  if(config) {
    if(config->cell_types) {
      for(i = 0; i < config->num_rows; i++) {
        free(config->cell_types[i]);
      }
      free(config->cell_types);
    }

    if(config->branch_ids) {
      for(i = 0; i < config->num_rows; i++)
        free(config->branch_ids[i]);
      free(config->branch_ids);
    }

    if(config->A) {
      int node_count = pressure_node_count(config);
      for(i = 0; i < node_count; i++) {
        free(config->A[i]);
      }
      free(config->A);
    }

    if(config->b) {
      free(config->b);
    }

    if(config->mapping) {
      for(i = 0; i < config->num_rows; i++) {
        free(config->mapping[i]);
      }
      free(config->mapping);
    }

    free(config);
  }
}
