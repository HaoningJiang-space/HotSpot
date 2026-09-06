#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>
#include <time.h>
#ifdef _MSC_VER
#define strcasecmp    _stricmp
#define strncasecmp   _strnicmp
#else
#include <strings.h>
#endif
#include <math.h>

#include "temperature_grid.h"
#include "flp.h"
#include "util.h"

#ifndef GATE_GS0
#define GATE_GS0 0
#endif

#ifndef GATE_GS1_ARM
#define GATE_GS1_ARM 0
#endif

#if GATE_GS0 > 0 && GATE_GS1_ARM > 0
#error GATE_GS0 and GATE_GS1_ARM cannot be enabled together
#endif

#define GS1_ARM_GS_NATIVE 1
#define GS1_ARM_SUPERLU 3
#define GS1_ARM_MFPCG 4
#define GS1_ARM_CSR_PCG 5
#define GS1_ARM_STRUCTURED_PCG 6

/* DSE-level candidate stopping contract. Builds may override it explicitly. */
#ifndef PCG_RELATIVE_RESIDUAL
#define PCG_RELATIVE_RESIDUAL 1.0e-4
#endif

#if GATE_GS1_ARM > 0
static double gs1_program_start_s;

static double gs1_now(void)
{
  struct timespec value;
  if (clock_gettime(CLOCK_MONOTONIC, &value) != 0)
    fatal("unable to read G-S1 monotonic clock\n");
  return (double)value.tv_sec + (double)value.tv_nsec * 1.0e-9;
}

__attribute__((constructor)) static void gs1_capture_program_start(void)
{
  struct timespec value;
  if (clock_gettime(CLOCK_MONOTONIC, &value) == 0)
    gs1_program_start_s = (double)value.tv_sec + (double)value.tv_nsec * 1.0e-9;
}

static const char *gs1_arm_name(void);
static const char *gs1_output_prefix(void);
static void gs1_write_state(grid_model_t *model, const double *block_state,
                            const char *prefix);
static void gs1_write_metadata(grid_model_t *model, const char *prefix);
#endif

#if SUPERLU < 1 || GATE_GS1_ARM == GS1_ARM_CSR_PCG || \
    GATE_GS1_ARM == GS1_ARM_STRUCTURED_PCG
static void jacobi_pcg_steady_grid(grid_model_t *model,
                                   grid_model_vector_t *power,
                                   grid_model_vector_t *temp);
#if GATE_GS1_ARM == GS1_ARM_GS_NATIVE
static void gs1_gauss_seidel_steady_grid(grid_model_t *model,
                                         grid_model_vector_t *power,
                                         grid_model_vector_t *temp);
#endif
#endif

static double relative_difference(double a, double b);

#if SUPERLU > 0
/* Lib for SuperLU */
#include "slu_ddefs.h"
#endif

/*BU_3D: We modified R-C computations to return resistance or capacitance for a specified grid cell.
 * If the grid cell does not have a unique value assigned in the flp file, we return the default values 
 * for that layer.the find_res_3D function will return the rx,ry,rz value of the grid.  
 * This function is used with the macros later.
 * It will calculate the joint resistance only if there are values defined within the array or rx,ry,rz values. */
double find_res_3D(int n, int i, int j, grid_model_t *model,int choice)
{
  int hasRes = model->layers[n].b2gmap[i][j]->hasRes;
  //Returns the rx of the grid cell
  if(choice==1){
      if(!hasRes)			
        return model->layers[n].rx;
      else
        return model->layers[n].b2gmap[i][j]->rx;
  }

  //Returns the ry of the grid cell
  else if(choice==2){	
      if(!hasRes)			
        return model->layers[n].ry;
      else
        return model->layers[n].b2gmap[i][j]->ry;
  }

  //Returns the rz of the grid cell
  else if(choice==3){	
      if(!hasRes)			
        return model->layers[n].rz;
      else
        return model->layers[n].b2gmap[i][j]->rz;
  }

  return 0;
}//end->BU_3D

/* BU_3D: finds capacitance of 3D cell.*/
double find_cap_3D(int n, int i, int j, grid_model_t *model)
{
  if (model->layers[n].b2gmap[i][j]->lock == TRUE) {
      // Return the capacitance of the unit that meets the occupancy threshold
      return model->layers[n].b2gmap[i][j]->capacitance;
  } else {
      // No unit that meets the occupancy threshold.
      // Return the layer's (default) capacitance instead
      return model->layers[n].c;
  }
}//end->BU_3D

enum {
  CONDUCTANCE_NORTH,
  CONDUCTANCE_SOUTH,
  CONDUCTANCE_EAST,
  CONDUCTANCE_WEST,
  CONDUCTANCE_ABOVE,
  CONDUCTANCE_BELOW,
  CONDUCTANCE_COUNT
};

#define STEADY_CONDUCTANCE(model, n, i, j, direction) \
  ((model)->steady_conductance[((((n) * (model)->rows + (i)) * \
                                 (model)->cols + (j)) * CONDUCTANCE_COUNT) + \
                               (direction)])

static void populate_steady_conductance(grid_model_t *model)
{
  int n, i, j;
  int nl = model->n_layers;
  int nr = model->rows;
  int nc = model->cols;

  free_dvector(model->steady_conductance);
  model->steady_conductance = dvector(nl * nr * nc * CONDUCTANCE_COUNT);
  for (n = 0; n < nl; n++)
    for (i = 0; i < nr; i++)
      for (j = 0; j < nc; j++) {
        STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_NORTH) =
          (i > 0) ? 1.0 / find_res_3D(n, i-1, j, model, 2) : 0.0;
        STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_SOUTH) =
          (i < nr-1) ? 1.0 / find_res_3D(n, i+1, j, model, 2) : 0.0;
        STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_EAST) =
          (j < nc-1) ? 1.0 / find_res_3D(n, i, j+1, model, 1) : 0.0;
        STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_WEST) =
          (j > 0) ? 1.0 / find_res_3D(n, i, j-1, model, 1) : 0.0;
        STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_ABOVE) =
          (n > 0) ? 1.0 / find_res_3D(n-1, i, j, model, 3) : 0.0;
        STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_BELOW) =
          (n < nl-1) ? 1.0 / find_res_3D(n, i, j, model, 3) : 0.0;
      }
}

/* constructors	*/
/*BU_3D: 
 * - Added parameter do_detailed_3D to this function
 * - Assign grid specific resistivty and capacitance if unit occupying the grid have resistivity & capacitance values
 * - If occupancy is > OCCUPANCY_THRESHOLD% then use the values obtained from the resistivity & capacitance of the occupying unit */
blist_t *new_blist(int idx, double occupancy, double res, double specificHeat,int first,int do_detailed_3D, double cw, double ch, double thickness)
{
  blist_t *ptr = (blist_t *) calloc (1, sizeof(blist_t));
  if (!ptr)
    fatal("memory allocation error\n");
  ptr->idx = idx;
  ptr->occupancy = occupancy;
  ptr->next = NULL;
  /*BU_3D: 
   * - If occupancy is greater than OCCUPANCY_THRESHOLD% lock in thermal resistance values 
   * - Else assume 50/50 occupancy*/
  if(first && do_detailed_3D){
      if(occupancy >= OCCUPANCY_THRESHOLD){
          ptr->lock=TRUE;
          ptr->rx =  getr(1/res, cw, ch * thickness);	
          ptr->ry =  getr(1/res, ch, cw * thickness);
          ptr->rz =  getr(1/res, thickness, cw * ch);
          ptr->capacitance = getcap(specificHeat, thickness, cw * ch);	
      }
      else{
          ptr->lock=FALSE;
          ptr->rx = 1 / ((1 / getr(1 / res, cw, ch * thickness)) * occupancy);
          ptr->ry = 1 / ((1 / getr(1 / res, ch, cw * thickness)) * occupancy);
          ptr->rz = 1 / ((1 / getr(1 / res, thickness, cw * ch)) * occupancy);
          ptr->capacitance = getcap(specificHeat, thickness, cw * ch);	
      }
      return ptr;
  }
  /*end->BU_3D*/
  return ptr;
}

blist_t ***new_b2gmap(int rows, int cols)
{
  int i;
  blist_t ***b2gmap;

  b2gmap = (blist_t ***) calloc (rows, sizeof(blist_t **));
  b2gmap[0] = (blist_t **) calloc (rows * cols, sizeof(blist_t *));
  if (!b2gmap || !b2gmap[0])
    fatal("memory allocation error\n");

  for(i=1; i < rows; i++) 
    b2gmap[i] = b2gmap[0] + cols * i;

  return b2gmap;	
}

/* destructor	*/
void delete_b2gmap(blist_t ***b2gmap, int rows, int cols)
{
  int i, j;
  blist_t *ptr, *temp;

  /* free the linked list	*/
  for(i=0; i < rows; i++)
    for(j=0; j < cols; j++) {
        ptr = b2gmap[i][j];
        while(ptr) {
            temp = ptr->next;
            free(ptr);
            ptr = temp;
        }
    }

  /* free the array space	*/
  free(b2gmap[0]);
  free(b2gmap);
}

/* re-initialize */
void reset_b2gmap(grid_model_t *model, layer_t *layer)
{
  int i, j;
  blist_t *ptr, *temp;

  /* free the linked list	*/
  for(i=0; i < model->rows; i++)
    for(j=0; j < model->cols; j++) {
        ptr = layer->b2gmap[i][j];
        while(ptr) {
            temp = ptr->next;
            free(ptr);
            ptr = temp;
        }
        layer->b2gmap[i][j] = NULL;
    }
}

/* create a linked list node and append it at the end
 * BU_3D: The parameter int do_detailed_3D is added to this function*/
void blist_append(blist_t *head, int idx, double occupancy,double res, double specificHeat,int first, int do_detailed_3D,double cw, double ch, double thickness)
{
  /*BU_3D: 
   * - If a block occupies a grid cell by OCCUPANCY_THRESHOLD% or more, the grid cell gets that blocks thermal resistance values
   * 	otherwise, 50/50 sharing is assumed with whichever block is occupying the other portion.
   * - Add resistances in parallel */
  if(do_detailed_3D && (head->lock!=TRUE)){
      if(occupancy >= OCCUPANCY_THRESHOLD){
          head->lock=TRUE;
          head->rx =  getr(1/res, cw, ch * thickness);	
          head->ry =  getr(1/res, ch, cw * thickness);
          head->rz =  getr(1/res, thickness, cw * ch);
          head->capacitance = getcap(specificHeat, thickness, cw * ch);	
      }
      else{
          head->rx = 1/((1/head->rx) + ((1 / getr(1 / res, cw, ch * thickness)) * occupancy));
          head->ry = 1/((1/head->ry) + ((1 / getr(1 / res, ch, cw * thickness)) * occupancy));
          head->rz = 1/((1/head->rz) + ((1 / getr(1 / res, thickness, cw * ch)) * occupancy));
          head->lock=FALSE;
      }
  }
  /*end->BU_3D*/
  blist_t *tail = NULL;

  if(!head)
    fatal("blist_append called with empty list\n");

  /* traverse till the end	*/
  for(; head; head = head->next)
    tail = head;

  /* append
   * BU_3D: added do_detailed_3D to this function call*/
  tail->next = new_blist(idx, occupancy, res, specificHeat, first, do_detailed_3D, cw, ch, thickness);
}

/* compute the power/temperature average weighted by occupancies	*/
double blist_avg(blist_t *ptr, flp_t *flp, double *v, int type)
{
  double  val = 0.0;

  for(; ptr; ptr = ptr->next) {
      if (type == V_POWER)
        val += ptr->occupancy * v[ptr->idx] / (flp->units[ptr->idx].width * 
                                               flp->units[ptr->idx].height);
      else if (type == V_TEMP)		   
        val += ptr->occupancy * v[ptr->idx];
      else
        fatal("unknown vector type\n");
  }		

  return val;		   
}

/* setup the block and grid mapping data structures	*/
void set_bgmap(grid_model_t *model, layer_t *layer)
{
  /* i1, i2, j1 and j2 are indices of the boundary grid cells	*/
  int i, j, u, i1, i2, j1, j2;

  /* shortcuts for cell width(cw) and cell height(ch)	*/
  double cw = model->width / model->cols;
  double ch = model->height / model->rows;
  /* shortcut for unit resistivity & specific heat*/
  double sh,res;
  /* initialize	*/
  reset_b2gmap(model, layer);

  /* for each functional unit	*/
  for(u=0; u < layer->flp->n_units; u++) {
      /* shortcuts for unit boundaries	*/
      double lu = layer->flp->units[u].leftx;
      double ru = lu + layer->flp->units[u].width;
      double bu = layer->flp->units[u].bottomy;
      double tu = bu + layer->flp->units[u].height;

      /* top index (lesser row) = rows - ceil (topy / cell height)	*/
      i1 = model->rows - tolerant_ceil(tu/ch);
      /* bottom index (greater row) = rows - floor (bottomy / cell height)	*/
      i2 = model->rows - tolerant_floor(bu/ch);
      /* left index = floor (leftx / cell width)	*/
      j1 = tolerant_floor(lu/cw);
      /* right index = ceil (rightx / cell width)	*/
      j2 = tolerant_ceil(ru/cw);
      /* sanity check	*/
      if((i1 < 0) || (j1 < 0))
        fatal("negative grid cell start index!\n");
      if((i2 > model->rows) || (j2 > model->cols))
        fatal("grid cell end index out of bounds!\n");
      if((i1 >= i2) || (j1 >= j2))
        fatal("invalid floorplan spec or grid resolution\n");

      /* setup g2bmap	*/
      layer->g2bmap[u].i1 = i1;
      layer->g2bmap[u].i2 = i2;
      layer->g2bmap[u].j1 = j1;
      layer->g2bmap[u].j2 = j2;

      /* setup b2gmap	*/
      /* for each grid cell in this unit	*/
      for(i=i1; i < i2; i++) {
          for(j=j1; j < j2; j++) {
              /* grid cells fully overlapped by this unit	*/

              /*BU_3D*/
              // - Load values from floorplan into each grid
              if(layer->flp->units[u].hasRes && model->config.detailed_3D_used)
                res = layer->flp->units[u].resistivity;
              else
                res = 1/layer->k;
              if(layer->flp->units[u].hasSh && model->config.detailed_3D_used)
                sh = layer->flp->units[u].specificheat;
              else
                sh = layer->sp;
              /*end->BU_3D*/

              if ((i > i1) && (i < i2-1) && (j > j1) && (j < j2-1)) {
                  /* first unit in the list	*/
                  if (!layer->b2gmap[i][j]){
                      layer->b2gmap[i][j] = new_blist(u, 1.0,res,sh,1,model->config.detailed_3D_used,cw,ch,layer->thickness);
                      layer->b2gmap[i][j]->hasRes = TRUE;//BU_3D assign hasRes to the b2gdata structure
                      layer->b2gmap[i][j]->hasCap = layer->flp->units[u].hasSh;//BU_3D assign hasSh to the b2gdata structure
                  }
                  else {
                      /* this should not occur since the grid cell is 
                       * fully covered and hence, no other unit should 
                       * be sharing it */
                      blist_append(layer->b2gmap[i][j], u, 1.0,res,sh,1,model->config.detailed_3D_used,cw,ch,layer->thickness);
                      warning("overlap of functional blocks?\n");
                  }
                  /* boundary grid cells partially overlapped by this unit	*/
              } else {
                  /* shortcuts for cell boundaries	*/
                  double lc = j * cw, rc = (j+1) * cw;
                  double tc = model->height - i * ch;
                  double bc = model->height - (i+1) * ch;

                  /* shortcuts for overlap width and height	*/
                  double oh = (MIN(tu, tc) - MAX(bu, bc));
                  double ow = (MIN(ru, rc) - MAX(lu, lc));
                  double occupancy;

                  /* overlap tolerance	*/
                  if (eq(oh/ch, 0))
                    oh = 0;
                  else if (eq(oh/ch, 1))
                    oh = ch;

                  if (eq(ow/cw, 0))
                    ow = 0;
                  else if (eq(ow/cw, 1))
                    ow = cw;

                  occupancy = (oh * ow) / (ch * cw);
                  if (oh < 0 || ow < 0)
                    fatal("negative overlap!\n");

                  /* first unit in the list	*/
                  if (!layer->b2gmap[i][j]){
                      layer->b2gmap[i][j] = new_blist(u, occupancy,res,sh,1,model->config.detailed_3D_used,cw,ch,layer->thickness);
                      layer->b2gmap[i][j]->hasRes = TRUE;//BU_3D assign hasRes to the b2gdata structure
                      layer->b2gmap[i][j]->hasCap = layer->flp->units[u].hasSh;//BU_3D assign hasSh to the b2gdata structure
                  }
                  else
                    /* append at the end */
                    blist_append(layer->b2gmap[i][j], u, occupancy,res,sh,0,model->config.detailed_3D_used,cw,ch,layer->thickness);
              }
          }
      }
  }
}

/* populate default set of layers	*/ 
void populate_default_layers(grid_model_t *model, flp_t *flp_default)
{
  int silidx, intidx, metalidx, c4idx;

  if (!model->config.model_secondary) {
      silidx = LAYER_SI;
      intidx = LAYER_INT;
  } else {
      silidx = SEC_PACK_LAYERS + SEC_CHIP_LAYERS + LAYER_SI;
      intidx = SEC_PACK_LAYERS + SEC_CHIP_LAYERS + LAYER_INT;
      c4idx  = SEC_PACK_LAYERS + LAYER_C4;
      metalidx = SEC_PACK_LAYERS + LAYER_METAL;		
  }

  /* silicon */
  model->layers[silidx].no = silidx;
  model->layers[silidx].has_lateral = TRUE;
  model->layers[silidx].has_power = TRUE;
  model->layers[silidx].k = model->config.k_chip;
  model->layers[silidx].thickness = model->config.t_chip;
  model->layers[silidx].sp = model->config.p_chip;
  model->layers[silidx].flp = flp_default;
  model->layers[silidx].b2gmap = new_b2gmap(model->rows, model->cols);
  model->layers[silidx].g2bmap = (glist_t *) calloc(flp_default->n_units, sizeof(glist_t));
  if (!model->layers[silidx].g2bmap)
    fatal("memory allocation error\n");

  /* interface material	*/
  model->layers[intidx].no = intidx;
  model->layers[intidx].has_lateral = TRUE;
  model->layers[intidx].has_power = FALSE;
  model->layers[intidx].k = model->config.k_interface;
  model->layers[intidx].thickness = model->config.t_interface;
  model->layers[intidx].sp = model->config.p_interface;
  model->layers[intidx].flp = flp_default;
  model->layers[intidx].b2gmap = model->layers[silidx].b2gmap;
  model->layers[intidx].g2bmap = model->layers[silidx].g2bmap;

  if (model->config.model_secondary) {
      /* metal layer	*/
      model->layers[metalidx].no = metalidx;
      model->layers[metalidx].has_lateral = TRUE;
      model->layers[metalidx].has_power = FALSE;
      model->layers[metalidx].k = K_METAL;
      model->layers[metalidx].thickness = model->config.t_metal;
      model->layers[metalidx].sp = SPEC_HEAT_METAL;
      model->layers[metalidx].flp = flp_default;
      model->layers[metalidx].b2gmap = model->layers[silidx].b2gmap;
      model->layers[metalidx].g2bmap = model->layers[silidx].g2bmap;

      /* C4/underfill layer*/
      model->layers[c4idx].no = c4idx;
      model->layers[c4idx].has_lateral = TRUE;
      model->layers[c4idx].has_power = FALSE;
      model->layers[c4idx].k = K_C4;
      model->layers[c4idx].thickness = model->config.t_c4;
      model->layers[c4idx].sp = SPEC_HEAT_C4;
      model->layers[c4idx].flp = flp_default;
      model->layers[c4idx].b2gmap = model->layers[silidx].b2gmap;
      model->layers[c4idx].g2bmap = model->layers[silidx].g2bmap;
  }
}

/* populate the package layers	*/ 
void append_package_layers(grid_model_t *model)
{
  /* shortcut	*/
  int nl = model->n_layers;
  int silidx, spidx, hsidx, subidx, solderidx, pcbidx;

  spidx = nl - DEFAULT_PACK_LAYERS + LAYER_SP;
  hsidx = nl - DEFAULT_PACK_LAYERS + LAYER_SINK;

  /* spreader	*/
  model->layers[spidx].no = spidx;
  model->layers[spidx].has_lateral = TRUE;
  model->layers[spidx].has_power = FALSE;
  model->layers[spidx].k = model->config.k_spreader;
  model->layers[spidx].thickness = model->config.t_spreader;
  model->layers[spidx].sp = model->config.p_spreader;
  model->layers[spidx].flp = model->layers[spidx-1].flp;
  model->layers[spidx].b2gmap = model->layers[spidx-1].b2gmap;
  model->layers[spidx].g2bmap = model->layers[spidx-1].g2bmap;

  /* heatsink	*/
  model->layers[hsidx].no = hsidx;
  model->layers[hsidx].has_lateral = TRUE;
  model->layers[hsidx].has_power = FALSE;
  model->layers[hsidx].k = model->config.k_sink;
  model->layers[hsidx].thickness = model->config.t_sink;
  model->layers[hsidx].sp = model->config.p_sink;
  model->layers[hsidx].flp = model->layers[spidx-1].flp;
  model->layers[hsidx].b2gmap = model->layers[spidx-1].b2gmap;
  model->layers[hsidx].g2bmap = model->layers[spidx-1].g2bmap;

  /* BU_3D: If detailed_3d option is on, the hasRes & hasCap variables must be set false in the b2gmap*/
  if(model->config.detailed_3D_used){
      model->layers[spidx].b2gmap = new_b2gmap(model->rows, model->cols);
      model->layers[hsidx].b2gmap = model->layers[spidx].b2gmap;
  }
  /*end->BU_3D*/

  if (model->config.model_secondary) {
      if(model->has_lcf)
        silidx = SEC_PACK_LAYERS;
      else
        silidx = SEC_PACK_LAYERS + SEC_CHIP_LAYERS + LAYER_SI;
      subidx = LAYER_SUB;
      solderidx = LAYER_SOLDER;
      pcbidx = LAYER_PCB;	

      /* package substrate	*/
      model->layers[subidx].no = subidx;
      model->layers[subidx].has_lateral = TRUE;
      model->layers[subidx].has_power = FALSE;
      model->layers[subidx].k = K_SUB;
      model->layers[subidx].thickness = model->config.t_sub;
      model->layers[subidx].sp = SPEC_HEAT_SUB;
      model->layers[subidx].flp = model->layers[silidx].flp;
      model->layers[subidx].b2gmap = model->layers[silidx].b2gmap;
      model->layers[subidx].g2bmap = model->layers[silidx].g2bmap;

      /* solder balls	*/
      model->layers[solderidx].no = solderidx;
      model->layers[solderidx].has_lateral = TRUE;
      model->layers[solderidx].has_power = FALSE;
      model->layers[solderidx].k = K_SOLDER;
      model->layers[solderidx].thickness = model->config.t_solder;
      model->layers[solderidx].sp = SPEC_HEAT_SOLDER;
      model->layers[solderidx].flp = model->layers[silidx].flp;
      model->layers[solderidx].b2gmap = model->layers[silidx].b2gmap;
      model->layers[solderidx].g2bmap = model->layers[silidx].g2bmap;

      /* PCB	*/
      model->layers[pcbidx].no = pcbidx;
      model->layers[pcbidx].has_lateral = TRUE;
      model->layers[pcbidx].has_power = FALSE;
      model->layers[pcbidx].k = K_PCB;
      model->layers[pcbidx].thickness = model->config.t_pcb;
      model->layers[pcbidx].sp = SPEC_HEAT_PCB;
      model->layers[pcbidx].flp = model->layers[silidx].flp;
      model->layers[pcbidx].b2gmap = model->layers[silidx].b2gmap;
      model->layers[pcbidx].g2bmap = model->layers[silidx].g2bmap;

      if(model->config.detailed_3D_used){
          model->layers[subidx].b2gmap = new_b2gmap(model->rows, model->cols);
          model->layers[solderidx].b2gmap = model->layers[subidx].b2gmap;
          model->layers[pcbidx].b2gmap = model->layers[subidx].b2gmap;
      }
  }
}

/* parse the layer file open for reading	*/
void parse_layer_file(grid_model_t *model, FILE *fp)
{
  char line[LINE_SIZE], *ptr, cval;
  int count, i = 0, field = LCF_SNO, ival;
  double dval;
  int base, inner_layers;

  if (!model->config.model_secondary){
      base = 0;
      inner_layers = model->n_layers - DEFAULT_PACK_LAYERS;
  }
  else{
      base = SEC_PACK_LAYERS;
      inner_layers = model->n_layers - DEFAULT_PACK_LAYERS - SEC_PACK_LAYERS;
  }

  fseek(fp, 0, SEEK_SET);
  count = 0;
  while (!feof(fp) && count < (inner_layers * LCF_NPARAMS)) {
      fgets(line, LINE_SIZE, fp);
      if (feof(fp))
        break;

      /* ignore comments and empty lines	*/
      ptr = strtok(line, " \r\t\n");
      if (!ptr || ptr[0] == '#')
        continue;

      switch (field) 
        {
        case LCF_SNO:
          if (sscanf(ptr, "%d", &ival) != 1)
            fatal("invalid layer number\n");
          if(ival >= inner_layers || ival < 0)
            fatal("layer number must be >= 0 and < no. of layers specified\n");
          if (model->layers[base+ival].no != 0)
            fatal("layer numbers must be unique\n");
          i = base + ival;
          model->layers[i].no = base + ival;
          break;
        case LCF_LATERAL:
          if (sscanf(ptr, "%c", &cval) != 1)
            fatal("invalid layer heat flow indicator\n");
          if (cval == 'Y' || cval == 'y')
            model->layers[i].has_lateral = TRUE;
          else if (cval == 'N' || cval == 'n')			
            model->layers[i].has_lateral = FALSE;
          else			
            fatal("invalid layer heat flow indicator\n");
          break;
        case LCF_POWER:
          if (sscanf(ptr, "%c", &cval) != 1)
            fatal("invalid layer power dissipation indicator\n");
          if (cval == 'Y' || cval == 'y')
            model->layers[i].has_power = TRUE;
          else if (cval == 'N' || cval == 'n')			
            model->layers[i].has_power = FALSE;
          else			
            fatal("invalid layer power dissipation indicator\n");
          break;
        case LCF_SP:
          if (sscanf(ptr, "%lf", &dval) != 1)
            fatal("invalid specific heat\n");
          model->layers[i].sp = dval;
          break;
        case LCF_RHO:
          if (sscanf(ptr, "%lf", &dval) != 1)
            fatal("invalid resistivity\n");
          model->layers[i].k = 1.0 / dval;
          break;
        case LCF_THICK:
          if (sscanf(ptr, "%lf", &dval) != 1)
            fatal("invalid thickness\n");
          model->layers[i].thickness = dval;
          break;
        case LCF_FLP:
          model->layers[i].flp = read_flp(ptr, FALSE);
          /* first layer	*/
          if (count < LCF_NPARAMS) {
              model->width = get_total_width(model->layers[i].flp);
              model->height = get_total_height(model->layers[i].flp);
          } else if(!eq(model->width, get_total_width(model->layers[i].flp)) || 
                    !eq(model->height, get_total_height(model->layers[i].flp)))
            fatal("width and height differ across layers\n");
          break;
        default:
          fatal("invalid field id\n");
          break;
        }
      field = (field + 1) % LCF_NPARAMS;
      count++;
  }

  /* allocate the block-grid maps */
  for(i=base; i < (base+inner_layers); i++) {
      model->layers[i].b2gmap = new_b2gmap(model->rows, model->cols);
      model->layers[i].g2bmap = (glist_t *) calloc(model->layers[i].flp->n_units, 
                                                   sizeof(glist_t));
      if (!model->layers[i].g2bmap)
        fatal("memory allocation error\n");
  }
}

/* populate layer info either from the default floorplan or from
 * the layer configuration file (lcf)
 */
void populate_layers_grid(grid_model_t *model, flp_t *flp_default)
{
  char str[STR_SIZE];
  FILE *fp = NULL;

  /* lcf file specified	*/
  if (model->has_lcf) {
      if (!strcasecmp(model->config.grid_layer_file, "stdin"))
        fp = stdin;
      else
        fp = fopen (model->config.grid_layer_file, "r");
      if (!fp) {
          sprintf(str, "error opening file %s\n", model->config.grid_layer_file);
          fatal(str);
      }
  }

  /* compute the no. of layers	*/
  if (!model->config.model_secondary) {
      if (model->has_lcf) {
          model->n_layers = count_significant_lines(fp);
          if (model->n_layers % LCF_NPARAMS)
            fatal("wrong no. of lines in layer file\n");
          model->n_layers /= LCF_NPARAMS;
          /* default no. of layers when lcf file is not specified	*/	
      } else
        model->n_layers = DEFAULT_CHIP_LAYERS;
  } else {
      if (model->has_lcf) {
          model->n_layers = count_significant_lines(fp);
          if (model->n_layers % LCF_NPARAMS)
            fatal("wrong no. of lines in layer file\n");
          model->n_layers /= LCF_NPARAMS;
          /* default no. of layers when lcf file is not specified	*/	
      } else
        model->n_layers = DEFAULT_CHIP_LAYERS + SEC_CHIP_LAYERS;
  }

  /* allocate initial memory including package layers	*/
  if (!model->config.model_secondary) {
      model->n_layers += DEFAULT_PACK_LAYERS;
      model->layers = (layer_t *) calloc (model->n_layers, sizeof(layer_t));
      if (!model->layers)
        fatal("memory allocation error\n");
  } else {
      model->n_layers += DEFAULT_PACK_LAYERS + SEC_PACK_LAYERS;
      model->layers = (layer_t *) calloc (model->n_layers, sizeof(layer_t));
      if (!model->layers)
        fatal("memory allocation error\n");
  }

  /* read in values from the lcf when specified	*/
  if (model->has_lcf) {
      parse_layer_file(model, fp);
      warning("layer configuration file specified. overriding default floorplan with those in lcf file...\n");
  } else {
      /* default set of layers	*/
      populate_default_layers(model, flp_default);
  }

  /* append the package layers	*/
  append_package_layers(model);

  if (model->has_lcf && fp != stdin)
    fclose(fp);
}

/* constructor */
grid_model_t *alloc_grid_model(thermal_config_t *config, flp_t *flp_default, int do_detailed_3D)
{
  int i;
  grid_model_t *model;

#if SUPERLU < 1
  if (config->grid_rows & (config->grid_rows-1) ||
      config->grid_cols & (config->grid_cols-1))
    fatal("grid rows and columns should both be powers of two\n");
#endif

  model = (grid_model_t *) calloc (1, sizeof(grid_model_t));
  if (!model)
    fatal("memory allocation error\n");
#if GATE_GS1_ARM > 0
  model->gs1_process_start_s = gs1_program_start_s;
  model->gs1_relative_residual = NAN;
  model->gs1_native_delta = NAN;
  model->gs1_action_audit_max_abs_error = NAN;
  model->gs1_action_audit_max_relative_error = NAN;
#endif
  model->config = *config;
  model->rows = config->grid_rows;
  model->cols = config->grid_cols;
  if(do_detailed_3D) //BU_3D: check if heterogenous RC model is on
    model->config.detailed_3D_used = TRUE; 
  if(!strcasecmp(model->config.grid_map_mode, GRID_AVG_STR))
    model->map_mode = GRID_AVG;
  else if(!strcasecmp(model->config.grid_map_mode, GRID_MIN_STR))
    model->map_mode = GRID_MIN;
  else if(!strcasecmp(model->config.grid_map_mode, GRID_MAX_STR))
    model->map_mode = GRID_MAX;
  else if(!strcasecmp(model->config.grid_map_mode, GRID_CENTER_STR))
    model->map_mode = GRID_CENTER;
  else
    fatal("unknown mapping mode\n");

  /* layer configuration file specified?	*/
  if(strcmp(model->config.grid_layer_file, NULLFILE))
    model->has_lcf = TRUE;
  else {
      model->has_lcf = FALSE;
      model->base_n_units = flp_default->n_units;
  }

  /* get layer information	*/
  populate_layers_grid(model, flp_default);

  /* count the total no. of blocks */
  model->total_n_blocks = 0;
  for(i=0; i < model->n_layers; i++)
    model->total_n_blocks += model->layers[i].flp->n_units;

  /* allocate internal state	*/
  model->last_steady = new_grid_model_vector(model);
  model->last_trans = new_grid_model_vector(model);

  return model;
}
#if DEBUG3D > 0
void debug_print_layer_det3D(grid_model_t *model, layer_t *layer);
#endif 
void populate_R_model_grid(grid_model_t *model, flp_t *flp)
{
  int i, base;
  double cw, ch;
#if GATE_GS1_ARM > 0
  double gs1_started_s = gs1_now();
#endif

  int inner_layers;
  int silidx, hsidx, pcbidx;
  int model_secondary = model->config.model_secondary;
  int nl = model->n_layers;

  if(model_secondary){
      if(model->has_lcf){
          base = SEC_PACK_LAYERS;
          inner_layers = nl - DEFAULT_PACK_LAYERS - SEC_PACK_LAYERS;
          silidx = SEC_PACK_LAYERS + LAYER_SI;
      }
      else{
          base = SEC_CHIP_LAYERS + SEC_PACK_LAYERS;
          inner_layers = nl - DEFAULT_PACK_LAYERS - SEC_PACK_LAYERS - SEC_CHIP_LAYERS;
          silidx = SEC_PACK_LAYERS + SEC_CHIP_LAYERS + LAYER_SI;
      }
      pcbidx = LAYER_PCB;	
  }
  else{
      base = 0;
      inner_layers = nl - DEFAULT_PACK_LAYERS;
      silidx = LAYER_SI;
  }

  /* setup the block-grid maps; flp parameter is ignored */
  if(model->has_lcf)
    for(i=base; i < (base+inner_layers); i++){
        set_bgmap(model, &model->layers[i]);
        //BU_3D Print the thermal resistance of each grid cell in form
        // 	row,col,rx,ry,rz		output as CSV file
#if DEBUG3D > 0
        debug_print_layer_det3D(model,&model->layers[i]);
        fprintf(stdout,"\n");
#endif
    }
  /* only the silicon layer has allocated space for the maps. 
   * all the rest just point to it. so it is sufficient to
   * setup the block-grid map for the silicon layer alone.
   * further, for default layer configuration, the `flp' 
   * parameter should be the same as that of the silicon 
   * layer. finally, the chip width and height information 
   * need to be populated for default layer configuration
   */
  else {
      if (flp != model->layers[silidx].flp)
        fatal("mismatch between the floorplan and the thermal model\n");
      model->width = get_total_width(flp);
      model->height = get_total_height(flp);
      set_bgmap(model, &model->layers[silidx]);
  }

  /* BU_3D: If detailed 3D is used we need to set variables in the b2gmap of the spreader layer*/
  if(model->config.detailed_3D_used){
      int ii,jj;
      int spidx = nl - DEFAULT_PACK_LAYERS + LAYER_SP;
      set_bgmap(model, &model->layers[spidx]);
      for(ii=0; ii < model->rows; ii++)
        for(jj=0; jj < model->cols; jj++){
            /* Make hasRes & hasCap as false*/
            model->layers[spidx].b2gmap[ii][jj]->hasRes = FALSE;
            model->layers[spidx].b2gmap[ii][jj]->hasCap = FALSE;
            model->layers[spidx].b2gmap[ii][jj]->lock   = FALSE;
        }
      if(model_secondary){
          int subidx = LAYER_SUB;
          set_bgmap(model, &model->layers[subidx]);
          for(ii=0; ii < model->rows; ii++)
            for(jj=0; jj < model->cols; jj++){
                /* Make hasRes & hasCap as false*/
                model->layers[subidx].b2gmap[ii][jj]->hasRes = FALSE;
                model->layers[subidx].b2gmap[ii][jj]->hasCap = FALSE;
                model->layers[subidx].b2gmap[ii][jj]->lock   = FALSE;
            }
      }
  } //end->BU_3D

  /* sanity check on floorplan sizes	*/
  if (model->width > model->config.s_sink || 
      model->height > model->config.s_sink || 
      model->width > model->config.s_spreader || 
      model->height > model->config.s_spreader) {
      print_flp(model->layers[silidx].flp);
      print_flp_fig(model->layers[silidx].flp);
      fatal("inordinate floorplan size!\n");
  }

  /* shortcuts for cell width(cw) and cell height(ch)	*/
  cw = model->width / model->cols;
  ch = model->height / model->rows;

  /* package R's	*/
  populate_package_R(&model->pack, &model->config, model->width, model->height);

  hsidx = nl - DEFAULT_PACK_LAYERS + LAYER_SINK;
  /* layer specific resistances	*/
  for(i=0; i < model->n_layers; i++){
      if (model->layers[i].has_lateral) {
          model->layers[i].rx =  getr(model->layers[i].k, cw, ch * model->layers[i].thickness);
          model->layers[i].ry =  getr(model->layers[i].k, ch, cw * model->layers[i].thickness);
      } else {
          /* positive infinity	*/
          model->layers[i].rx = LARGENUM;
          model->layers[i].ry = LARGENUM;
      }
      model->layers[i].rz =  getr(model->layers[i].k, model->layers[i].thickness, cw * ch);

      /* heatsink	is connected to ambient. divide r_convec proportional to cell area */
      if (i == hsidx){
          model->layers[i].rz += model->config.r_convec * 
            (model->config.s_sink * model->config.s_sink) / (cw * ch);
      }
      //if (model_secondary && (i == pcbidx)){
      //    /* PCB	is connected to ambient. divide r_convec_sec proportional to cell area */
      //      model->layers[i].rz += model->config.r_convec_sec * 
      //        (model->config.s_pcb * model->config.s_pcb) / (cw * ch);
      //}
  }

  /* done	*/
  if (model->config.detailed_3D_used &&
      GATE_GS1_ARM != GS1_ARM_GS_NATIVE &&
      GATE_GS1_ARM != GS1_ARM_SUPERLU)
    populate_steady_conductance(model);
  model->r_ready = TRUE;
#if GATE_GS1_ARM > 0
  model->gs1_coefficient_s += gs1_now() - gs1_started_s;
#endif
}

void populate_C_model_grid(grid_model_t *model, flp_t *flp)
{
  int i;
  int silidx, hsidx, pcbidx;

  /* shortcuts */
  double cw = model->width / model->cols;
  double ch = model->height / model->rows;
  int nl = model->n_layers;
  int model_secondary = model->config.model_secondary;

  hsidx = nl - DEFAULT_PACK_LAYERS + LAYER_SINK;
  if(model_secondary){
      if(model->has_lcf)
        silidx = SEC_PACK_LAYERS + LAYER_SI;
      else
        silidx = SEC_PACK_LAYERS + SEC_CHIP_LAYERS + LAYER_SI;

      pcbidx = LAYER_PCB;	
  }
  else{
      silidx = LAYER_SI;
  }

  if (!model->r_ready)
    fatal("R model not ready\n");
  if (!model->has_lcf && flp != model->layers[silidx].flp)
    fatal("different floorplans for R and C models!\n");

  /* package C's	*/
  populate_package_C(&model->pack, &model->config, model->width, model->height);

  /* layer specific capacitances	*/
  for(i=0; i < nl; i++){
      model->layers[i].c =  getcap(model->layers[i].sp, model->layers[i].thickness, cw * ch);
  }

  /* last layer (heatsink) is connected to the ambient. 
   * divide c_convec proportional to cell area 
   */
  model->layers[hsidx].c += C_FACTOR * model->config.c_convec * (cw * ch) / \
                            (model->config.s_sink * model->config.s_sink);
  if (model_secondary) {		
      model->layers[pcbidx].c += C_FACTOR * model->config.c_convec_sec * (cw * ch) / \
                                 (model->config.s_pcb * model->config.s_pcb);
  }

  /* done	*/	
  model->c_ready = TRUE;
}

/* destructor	*/
void delete_grid_model(grid_model_t *model)
{
  int i, base;

  int inner_layers;
  int silidx;
  int model_secondary = model->config.model_secondary;
  int nl = model->n_layers;

  if (model_secondary){
      if(model->has_lcf){
          base = SEC_PACK_LAYERS;
          inner_layers = nl - DEFAULT_PACK_LAYERS - SEC_PACK_LAYERS;
          silidx = SEC_PACK_LAYERS + LAYER_SI;
      }
      else{
          base = SEC_CHIP_LAYERS + SEC_PACK_LAYERS;
          inner_layers = nl - DEFAULT_PACK_LAYERS - SEC_PACK_LAYERS - SEC_CHIP_LAYERS;
          silidx = SEC_PACK_LAYERS + SEC_CHIP_LAYERS + LAYER_SI;
      }
  }
  else{
      base = 0;
      inner_layers = nl - DEFAULT_PACK_LAYERS;
      silidx = LAYER_SI;
  }

  if (model->has_lcf) 
    for(i=base; i < (base+inner_layers); i++){
        delete_b2gmap(model->layers[i].b2gmap, model->rows, model->cols);
        free(model->layers[i].g2bmap);
        free_flp(model->layers[i].flp, FALSE);
    }
  /* only the silicon layer has allocated space for the maps. 
   * all the rest just point to it. also, its floorplan was
   * allocated elsewhere. so, we don't need to deallocate those.
   */
  else {
      delete_b2gmap(model->layers[silidx].b2gmap, model->rows, model->cols);
      free(model->layers[silidx].g2bmap);
  }
  /*BU_3D: If detailed_3d is used then perform delete the b2gmap 
   * 	 that was allocated for the package */
  if(model->config.detailed_3D_used){
      int spidx = nl - DEFAULT_PACK_LAYERS + LAYER_SP;
      delete_b2gmap(model->layers[spidx].b2gmap, model->rows, model->cols);
      if(model_secondary){
          int subidx = LAYER_SUB;
          delete_b2gmap(model->layers[subidx].b2gmap, model->rows, model->cols);
      }
  }//end->BU_3D

  free_grid_model_vector(model->last_steady);
  free_grid_model_vector(model->last_trans);
  free_dvector(model->steady_conductance);
  free(model->layers);
  free(model);
}

/* differs from 'dvector()' in that memory for internal nodes is also allocated	*/
double *hotspot_vector_grid(grid_model_t *model)
{
  double total_nodes;
  if (model->total_n_blocks <= 0)
    fatal("total_n_blocks is not greater than zero\n");

  if (model->config.model_secondary)
    total_nodes = model->total_n_blocks + EXTRA + EXTRA_SEC;
  else
    total_nodes = model->total_n_blocks + EXTRA;

  return dvector(total_nodes);
}

/* copy 'src' to 'dst' except for a window of 'size'
 * elements starting at 'at'. useful in floorplan
 * compaction. can be used only with default layer
 * configuration as all layers should have the same
 * floorplan. incompatible with layer configuration
 * file.
 */
void trim_hotspot_vector_grid(grid_model_t *model, double *dst, double *src, 
                              int at, int size)
{
  int i;

  double total_nodes;
  if (model->config.model_secondary)
    total_nodes = model->total_n_blocks + EXTRA + EXTRA_SEC;
  else
    total_nodes = model->total_n_blocks + EXTRA;

  if (model->has_lcf)
    fatal("trim_hotspot_vector_grid called with lcf file\n");
  for (i=0; i < at && i < total_nodes; i++)
    dst[i] = src[i];
  for(i=at+size; i < total_nodes; i++)
    dst[i-size] = src[i];
}

/* update the model corresponding to floorplan compaction	*/						 
void resize_thermal_model_grid(grid_model_t *model, int n_units)
{
  int i;

  if (model->has_lcf)
    fatal("resize_thermal_model_grid called with lcf file\n");
  if (n_units > model->base_n_units)
    fatal("resizing grid model to more than the allocated space\n");

  /* count the total no. of blocks again */
  model->total_n_blocks = 0;
  for(i=0; i < model->n_layers; i++)
    model->total_n_blocks += model->layers[i].flp->n_units;

  /* nothing more needs to be done because the only data structure
   * that is dependent on flp->n_units is g2bmap (others are
   * dependent on 'grid size' which does not change because
   * of resizing). g2bmap is a 1-d array and needs no reallocation
   */
}

/* sets the temperature of a vector 'temp' allocated using 'hotspot_vector'	*/
void set_temp_grid(grid_model_t *model, double *temp, double val)
{
  int i;

  double total_nodes;
  if (model->config.model_secondary)
    total_nodes = model->total_n_blocks + EXTRA + EXTRA_SEC;
  else
    total_nodes = model->total_n_blocks + EXTRA;

  if (model->total_n_blocks <= 0)
    fatal("total_n_blocks is not greater than zero\n");
  for(i=0; i < total_nodes; i++)
    temp[i] = val;
}

/* dump the steady state grid temperatures of the top layer onto 'file'	*/
void dump_top_layer_temp_grid (grid_model_t *model, char *file, grid_model_vector_t *temp)
{
  int i, j;
  char str[STR_SIZE];
  FILE *fp;
  int silidx;

  if (!model->config.model_secondary) {
      silidx = LAYER_SI;
  } else {
      if(model->has_lcf)
        silidx = SEC_PACK_LAYERS;
      else
        silidx = SEC_PACK_LAYERS + SEC_CHIP_LAYERS;
  }

  if (!model->r_ready)
    fatal("R model not ready\n");

  if (!strcasecmp(file, "stdout"))
    fp = stdout;
  else if (!strcasecmp(file, "stderr"))
    fp = stderr;
  else 	
    fp = fopen (file, "w");

  if (!fp) {
      sprintf (str,"error: %s could not be opened for writing\n", file);
      fatal(str);
  }

  for(i=0; i < model->rows; i++){
      for(j=0; j < model->cols; j++){
          /* top layer of the most-recently computed 
           * steady state temperature	
           */
  //        fprintf(fp, "%d\t%d\t%.2f\n", j, model->rows-i-1, 
  //                model->last_steady->cuboid[silidx][i][j]); 
          fprintf(fp, "%d\t%.2f\n", i*model->cols+j, 
                  model->last_steady->cuboid[silidx][i][j]); 
      }
      fprintf(fp, "\n");
  }

  //int l;
  //for(l=0;  l < model->n_layers; l++){
  //    for(i=0; i < model->rows; i++){
  //        for(j=0; j < model->cols; j++){
  //            fprintf(fp, "%d\t%.2f\n", i*model->cols+j, 
  //                    model->last_steady->cuboid[l][i][j]); 
  //        }
  //    }
  //    fprintf(fp, "\n");
  //}
 
  //rzdebug
  //int l;
  //double max_temp = 0;
  //for(l=0; l<model->n_layers; l++)
  //  for(i=0; i<model->rows; i++)
  //    for(j=0; j<model->cols; j++)
  //      if(model->last_steady->cuboid[l][i][j] > max_temp)
  //        max_temp = model->last_steady->cuboid[l][i][j];
  //printf("%lf\n", max_temp);

  if(fp != stdout && fp != stderr)
    fclose(fp);	
}

/* dump the steady state grid temperatures of the top layer onto 'file'	*/
void dump_steady_temp_grid (grid_model_t *model, char *file)
{
  /* top layer of the most-recently computed steady state temperature	*/
  dump_top_layer_temp_grid(model, file, model->last_steady);
}

/* dump temperature vector alloced using 'hotspot_vector' to 'file' */ 
void dump_temp_grid(grid_model_t *model, double *temp, char *file)
{
  int i, n, base = 0;
  char str[STR_SIZE];
  FILE *fp;

  int extra_nodes;
  int model_secondary = model->config.model_secondary;
  int nl = model->n_layers;
  int spidx, hsidx, silidx, intidx, c4idx, metalidx, subidx, solderidx, pcbidx;

  if (model_secondary)
    extra_nodes = EXTRA + EXTRA_SEC;
  else
    extra_nodes = EXTRA;

  if (!strcasecmp(file, "stdout"))
    fp = stdout;
  else if (!strcasecmp(file, "stderr"))
    fp = stderr;
  else 	
    fp = fopen (file, "w");

  if (!fp) {
      sprintf (str,"error: %s could not be opened for writing\n", file);
      fatal(str);
  }

  spidx = nl - DEFAULT_PACK_LAYERS + LAYER_SP;
  hsidx = nl - DEFAULT_PACK_LAYERS + LAYER_SINK;
  if (model_secondary) {
      subidx = LAYER_SUB;
      solderidx = LAYER_SOLDER;
      pcbidx = LAYER_PCB;	
      if(!model->has_lcf){
          silidx = SEC_PACK_LAYERS + SEC_CHIP_LAYERS + LAYER_SI;
          intidx = SEC_PACK_LAYERS + SEC_CHIP_LAYERS + LAYER_INT;
          c4idx  = SEC_PACK_LAYERS + LAYER_C4;
          metalidx = SEC_PACK_LAYERS + LAYER_METAL;		
      }
  }
  else{
      silidx = LAYER_SI;
      intidx = LAYER_INT;
  }

  /* layer temperatures	*/
  for(n=0; n < model->n_layers; n++) {
      if (!model_secondary) {
          /* default set of layers	*/
          if (!model->has_lcf) {
              if(n == silidx)
                strcpy(str,"");
              else if(n == intidx)
                strcpy(str,"iface_");
              else if(n == spidx)
                strcpy(str,"hsp_");
              else if(n == hsidx)
                strcpy(str,"hsink_");
              else
                fatal("unknown layer\n");
          } else {
              if (n == spidx)
                strcpy(str, "hsp_");	/* spreader layer	*/
              else if (n == hsidx)
                strcpy(str, "hsink_");	/* heatsink layer	*/
              else	/* other layers	*/
                sprintf(str,"layer_%d_", n);
          }
      } else {
          /* default set of layers	*/
          if (!model->has_lcf) {
              if(n == silidx)
                strcpy(str,"");
              else if(n == intidx)
                strcpy(str,"iface_");
              else if(n == spidx)
                strcpy(str,"hsp_");
              else if(n == hsidx)
                strcpy(str,"hsink_");
              else if(n == metalidx)
                strcpy(str,"metal_");
              else if(n == c4idx)
                strcpy(str,"c4_");
              else if(n == subidx)
                strcpy(str,"sub_");
              else if(n == solderidx)
                strcpy(str,"solder_");
              else if(n == pcbidx)
                strcpy(str,"pcb_");
              else
                fatal("unknown layer\n");
              /* layer configuration file	*/
          } else {
              if (n == spidx)
                strcpy(str, "hsp_");	/* spreader layer	*/
              else if (n == hsidx)
                strcpy(str, "hsink_");	/* heatsink layer	*/
              else if (n == subidx)
                strcpy(str, "sub_");	/* package substrate layer	*/
              else if (n == solderidx)
                strcpy(str, "solder_");	/* solder layer	*/
              else if (n == pcbidx)
                strcpy(str, "pcb_");	/* pcb layer	*/
              else	/* other layers	*/
                sprintf(str,"layer_%d_", n);
          }
      }

      for(i=0; i < model->layers[n].flp->n_units; i++)
        fprintf(fp, "%s%s\t%.2f\n", str, 
                model->layers[n].flp->units[i].name, temp[base+i]);
      base += model->layers[n].flp->n_units;	
  }

  if (base != model->total_n_blocks)
    fatal("total_n_blocks failed to tally\n");

  /* internal node temperatures	*/
  for (i=0; i < extra_nodes; i++) {
      sprintf(str, "inode_%d", i);
      fprintf(fp, "%s\t%.2f\n", str, temp[base+i]);
  }
  if(fp != stdout && fp != stderr)
    fclose(fp);	
}

void copy_temp_grid(grid_model_t *model, double *dst, double *src)
{
  if (!model->config.model_secondary)
    copy_dvector(dst, src, model->total_n_blocks + EXTRA);
  else
    copy_dvector(dst, src, model->total_n_blocks + EXTRA + EXTRA_SEC);
}

/* 
 * read temperature vector alloced using 'hotspot_vector' from 'file'
 * which was dumped using 'dump_temp'. values are clipped to thermal
 * threshold based on 'clip'
 */ 
void read_temp_grid(grid_model_t *model, double *temp, char *file, int clip)
{
  int i, n, idx, base = 0;
  double max=0, val;
  char *ptr, str1[LINE_SIZE], str2[LINE_SIZE];
  char name[STR_SIZE], format[STR_SIZE];
  FILE *fp;

  int model_secondary = model->config.model_secondary;
  int extra_nodes;
  int nl = model->n_layers;
  int spidx, hsidx, silidx, intidx, c4idx, metalidx, subidx, solderidx, pcbidx;

  if (model->config.model_secondary)
    extra_nodes = EXTRA + EXTRA_SEC;
  else
    extra_nodes = EXTRA;

  if (!strcasecmp(file, "stdin"))
    fp = stdin;
  else
    fp = fopen (file, "r");

  if (!fp) {
      sprintf (str1,"error: %s could not be opened for reading\n", file);
      fatal(str1);
  }	

  spidx = nl - DEFAULT_PACK_LAYERS + LAYER_SP;
  hsidx = nl - DEFAULT_PACK_LAYERS + LAYER_SINK;
  if (model_secondary) {
      subidx = LAYER_SUB;
      solderidx = LAYER_SOLDER;
      pcbidx = LAYER_PCB;	
      if(!model->has_lcf){
          silidx = SEC_PACK_LAYERS + SEC_CHIP_LAYERS + LAYER_SI;
          intidx = SEC_PACK_LAYERS + SEC_CHIP_LAYERS + LAYER_INT;
          c4idx  = SEC_PACK_LAYERS + LAYER_C4;
          metalidx = SEC_PACK_LAYERS + LAYER_METAL;		
      }
  }
  else{
      silidx = LAYER_SI;
      intidx = LAYER_INT;
  }

  /* temperatures of the different layers	*/
  for (n=0; n < model->n_layers; n++) {
      if (!model_secondary) {
          /* default set of layers	*/
          if (!model->has_lcf) {
              if(n == silidx)
                strcpy(format,"%s%lf");
              else if(n == intidx)
                strcpy(format,"iface_%s%lf");
              else if(n == spidx)
                strcpy(format,"hsp_%s%lf");
              else if(n == hsidx)
                strcpy(format,"hsink_%s%lf");
              else
                fatal("unknown layer\n");
          } else {
              /* layer configuration file	*/
              if (n == spidx)
                strcpy(format, "hsp_%s%lf");	/* spreader layer	*/
              else if (n == hsidx)
                strcpy(format, "hsink_%s%lf");	/* heatsink layer	*/
              else	/* other layers	*/
                sprintf(format,"layer_%d_%%s%%lf", n);
          }
      } else {
          /* default set of layers	*/
          if (!model->has_lcf) {
              if(n == silidx)
                strcpy(format,"%s%lf");
              else if(n == intidx)
                strcpy(format,"iface_%s%lf");
              else if(n == spidx)
                strcpy(format,"hsp_%s%lf");
              else if(n == hsidx)
                strcpy(format,"hsink_%s%lf");
              else if(n == metalidx)
                strcpy(format,"metal_%s%lf");
              else if(n == c4idx)
                strcpy(format,"c4_%s%lf");
              else if(n == subidx)
                strcpy(format,"sub_%s%lf");
              else if(n == solderidx)
                strcpy(format,"solder_%s%lf");
              else if(n == pcbidx)
                strcpy(format,"pcb_%s%lf");
              else
                fatal("unknown layer\n");
              /* layer configuration file	*/
          } else {
              if (n == spidx)
                strcpy(format, "hsp_%s%lf");	/* spreader layer	*/
              else if (n == hsidx)
                strcpy(format, "hsink_%s%lf");	/* heatsink layer	*/
              else if (n == subidx)
                strcpy(format, "sub_%s%lf");	/* package substrate layer	*/
              else if (n == solderidx)
                strcpy(format, "solder_%s%lf");	/* solder ball layer	*/
              else if (n == pcbidx)
                strcpy(format, "pcb_%s%lf");	/* PCB layer	*/
              else	/* other layers	*/
                sprintf(format,"layer_%d_%%s%%lf", n);
          }
      }

      for (i=0; i < model->layers[n].flp->n_units; i++) {
          fgets(str1, LINE_SIZE, fp);
          if (feof(fp))
            fatal("not enough lines in temperature file\n");
          strcpy(str2, str1);
          /* ignore comments and empty lines	*/
          ptr = strtok(str1, " \r\t\n");
          if (!ptr || ptr[0] == '#') {
              i--;
              continue;
          }
          if (sscanf(str2, format, name, &val) != 2)
            fatal("invalid temperature file format\n");
          idx = get_blk_index(model->layers[n].flp, name);
          if (idx >= 0)
            temp[base+idx] = val;
          else	/* since get_blk_index calls fatal, the line below cannot be reached	*/
            fatal ("unit in temperature file not found in floorplan\n");

          /* find max temp on the top layer 
           * (silicon for the default set of layers)
           */
          if (n == 0 && temp[idx] > max)
            max = temp[idx];
      }
      base += model->layers[n].flp->n_units;
  }

  if (base != model->total_n_blocks)
    fatal("total_n_blocks failed to tally\n");

  /* internal node temperatures	*/
  for (i=0; i < extra_nodes; i++) {
      fgets(str1, LINE_SIZE, fp);
      if (feof(fp))
        fatal("not enough lines in temperature file\n");
      strcpy(str2, str1);
      /* ignore comments and empty lines	*/
      ptr = strtok(str1, " \r\t\n");
      if (!ptr || ptr[0] == '#') {
          i--;
          continue;
      }
      if (sscanf(str2, "%s%lf", name, &val) != 2)
        fatal("invalid temperature file format\n");
      sprintf(str1, "inode_%d", i);
      if (strcasecmp(str1, name))
        fatal("invalid temperature file format\n");
      temp[base+i] = val;	
  }

  fgets(str1, LINE_SIZE, fp);
  if (!feof(fp))
    fatal("too many lines in temperature file\n");

  if(fp != stdin)
    fclose(fp);	

  /* clipping	*/
  if (clip && (max > model->config.thermal_threshold)) {
      /* if max has to be brought down to thermal_threshold, 
       * (w.r.t the ambient) what is the scale down factor?
       */
      double factor = (model->config.thermal_threshold - model->config.ambient) / 
        (max - model->config.ambient);

      /* scale down all temperature differences (from ambient) by the same factor	*/
      for (i=0; i < model->total_n_blocks + extra_nodes; i++)
        temp[i] = (temp[i]-model->config.ambient)*factor + model->config.ambient;
  }
}

/* dump power numbers to file	*/
void dump_power_grid(grid_model_t *model, double *power, char *file)
{
  int i, n, base = 0;
  char str[STR_SIZE];
  FILE *fp;

  if (!strcasecmp(file, "stdout"))
    fp = stdout;
  else if (!strcasecmp(file, "stderr"))
    fp = stderr;
  else 	
    fp = fopen (file, "w");
  if (!fp) {
      sprintf (str,"error: %s could not be opened for writing\n", file);
      fatal(str);
  }

  /* dump values only for the layers dissipating power	*/
  for(n=0; n < model->n_layers; n++) {
      if (model->layers[n].has_power) {
          for(i=0; i < model->layers[n].flp->n_units; i++)
            if (model->has_lcf)
              fprintf(fp, "layer_%d_%s\t%.6f\n", n, 
                      model->layers[n].flp->units[i].name, power[base+i]);
            else 
              fprintf(fp, "%s\t%.6f\n", 
                      model->layers[n].flp->units[i].name, power[base+i]);
      }				
      base += model->layers[n].flp->n_units;		
  }

  if(fp != stdout && fp != stderr)
    fclose(fp);	
}

/* 
 * read power vector alloced using 'hotspot_vector' from 'file'
 * which was dumped using 'dump_power'. 
 */ 
void read_power_grid(grid_model_t *model, double *power, char *file)
{
  int i, idx, n, base = 0;
  double val;
  char *ptr, str1[LINE_SIZE], str2[LINE_SIZE]; 
  char name[STR_SIZE], format[STR_SIZE];
  FILE *fp;

  if (!strcasecmp(file, "stdin"))
    fp = stdin;
  else
    fp = fopen (file, "r");
  if (!fp) {
      sprintf (str1,"error: %s could not be opened for reading\n", file);
      fatal(str1);
  }

  /* lcf file could potentially specify more than one power dissipating 
   * layer. hence, units with zero power within a layer cannot be left
   * out in the power file.
   */
  if (model->has_lcf) {
      for(n=0; n < model->n_layers; n++) {
          if (model->layers[n].has_power)
            for(i=0; i < model->layers[n].flp->n_units; i++) {
                fgets(str1, LINE_SIZE, fp);
                if (feof(fp))
                  fatal("not enough lines in power file\n");
                strcpy(str2, str1);

                /* ignore comments and empty lines	*/
                ptr = strtok(str1, " \r\t\n");
                if (!ptr || ptr[0] == '#') {
                    i--;
                    continue;
                }

                sprintf(format, "layer_%d_%%s%%lf", n);
                if (sscanf(str2, format, name, &val) != 2)
                  fatal("invalid power file format\n");
                idx = get_blk_index(model->layers[n].flp, name);
                if (idx >= 0)
                  power[base+idx] = val;
                /* since get_blk_index calls fatal, the line below cannot be reached	*/
                else
                  fatal ("unit in power file not found in floorplan\n");
            }
          base += model->layers[n].flp->n_units;	
      }
      fgets(str1, LINE_SIZE, fp);
      if (!feof(fp))
        fatal("too many lines in power file\n");
      /* default layer configuration. so only one layer
       * has power dissipation. units with zero power 
       * can be omitted in the power file
       */
  } else {
      while(!feof(fp)) {
          fgets(str1, LINE_SIZE, fp);
          if (feof(fp))
            break;
          strcpy(str2, str1);

          /* ignore comments and empty lines	*/
          ptr = strtok(str1, " \r\t\n");
          if (!ptr || ptr[0] == '#')
            continue;

          if (sscanf(str2, "%s%lf", name, &val) != 2)
            fatal("invalid power file format\n");
          idx = get_blk_index(model->layers[LAYER_SI].flp, name);
          if (idx >= 0)
            power[idx] = val;
          else	/* since get_blk_index calls fatal, the line below cannot be reached	*/
            fatal ("unit in power file not found in floorplan\n");
      }
  }

  if(fp != stdin)
    fclose(fp);
}

double find_max_temp_grid(grid_model_t *model, double *temp)
{
  int i;
  double max = 0.0;
  int silidx;

  if (!model->config.model_secondary) {
      silidx = LAYER_SI;
  } else {
      if(model->has_lcf)
        silidx = SEC_PACK_LAYERS;
      else
        silidx = SEC_PACK_LAYERS + SEC_CHIP_LAYERS;
  }

  /* max temperature occurs on the top-most layer	*/
  for(i=0; i < model->layers[silidx].flp->n_units; i++) {
      if (temp[i] < 0)
        fatal("negative temperature!\n");
      else if (max < temp[i])
        max = temp[i];
  }

  return max;
}

double find_avg_temp_grid(grid_model_t *model, double *temp)
{
  int i, n, base = 0, count = 0; 
  double sum = 0.0;
  /* average temperature of all the power dissipating blocks	*/
  for(n=0; n < model->n_layers; n++) {
      if (model->layers[n].has_power) {
          for(i=0; i < model->layers[n].flp->n_units; i++) {
              if (temp[base+i] < 0)
                fatal("negative temperature!\n");
              else 
                sum += temp[base+i];
          }
          count += model->layers[n].flp->n_units;
      }	
      base += model->layers[n].flp->n_units;
  }

  if (!count)
    fatal("no power dissipating units?!\n");
  return (sum / count);
}

/* calculate average heatsink temperature for natural convection package */
double calc_sink_temp_grid(grid_model_t *model, double *temp, thermal_config_t *config)
{
  int i, n, base = 0; 
  int hsidx = model->n_layers - DEFAULT_PACK_LAYERS + LAYER_SINK;
  double sum = 0.0;
  double spr_size = config->s_spreader*config->s_spreader;
  double sink_size = config->s_sink*config->s_sink;

  /* heat sink core	*/
  for(n=0; n < hsidx; n++)
    base += model->layers[n].flp->n_units;

  for(i=base; i < base+model->layers[n].flp->n_units; i++)
    if (temp[i] < 0)
      fatal("negative temperature!\n");
    else /* area-weighted average */
      sum += temp[i]*(model->layers[n].flp->units[i-base].width*model->layers[n].flp->units[i-base].height);

  /* heat sink periphery	*/
  base = model->total_n_blocks;

  for(i=SINK_C_W; i <= SINK_C_E; i++)
    if (temp[i+base] < 0)
      fatal("negative temperature!\n");
    else
      sum += temp[i+base]*0.25*(config->s_spreader+model->height)*(config->s_spreader-model->width);

  for(i=SINK_C_N; i <= SINK_C_S; i++)
    if (temp[i+base] < 0)
      fatal("negative temperature!\n");
    else
      sum += temp[i+base]*0.25*(config->s_spreader+model->width)*(config->s_spreader-model->height);

  for(i=SINK_W; i <= SINK_S; i++)
    if (temp[i+base] < 0)
      fatal("negative temperature!\n");
    else
      sum += temp[i+base]*0.25*(sink_size-spr_size);

  return (sum / sink_size);
}

/* grid_model_vector routines	*/

/* constructor	*/
grid_model_vector_t *new_grid_model_vector(grid_model_t *model)
{
  grid_model_vector_t *v;

  int extra_nodes;
  if (model->config.model_secondary)
    extra_nodes = EXTRA + EXTRA_SEC;
  else
    extra_nodes = EXTRA;

  v = (grid_model_vector_t *) calloc (1, sizeof(grid_model_vector_t));
  if (!v)
    fatal("memory allocation error\n");

  v->cuboid = dcuboid_tail(model->rows, model->cols, model->n_layers, extra_nodes);
  v->extra = v->cuboid[0][0] + model->rows * model->cols * model->n_layers;
  return v;
}

/* destructor	*/
void free_grid_model_vector(grid_model_vector_t *v)
{
  free_dcuboid(v->cuboid);
  free(v);
}

/* translate power/temperature between block and grid vectors	*/
void xlate_vector_b2g(grid_model_t *model, double *b, grid_model_vector_t *g, int type)
{
  int i, j, n, base = 0;
  double area;

  int extra_nodes;
  if (model->config.model_secondary)
    extra_nodes = EXTRA + EXTRA_SEC;
  else
    extra_nodes = EXTRA;

  /* area of a single grid cell	*/
  area = (model->width * model->height) / (model->cols * model->rows);

  for(n=0; n < model->n_layers; n++) {
      for(i=0; i < model->rows; i++)
        for(j=0; j < model->cols; j++) {
            /* for each grid cell, the power density / temperature are 
             * the average of the power densities / temperatures of the 
             * blocks in it weighted by their occupancies
             */
            /* convert power density to power	*/ 
            if (type == V_POWER)
              g->cuboid[n][i][j] = blist_avg(model->layers[n].b2gmap[i][j], 
                                             model->layers[n].flp, &b[base], type) * area;
            /* no conversion necessary for temperature	*/ 
            else if (type == V_TEMP)
              g->cuboid[n][i][j] = blist_avg(model->layers[n].b2gmap[i][j], 
                                             model->layers[n].flp, &b[base], type);
            else
              fatal("unknown vector type\n");
        }
      /* keep track of the beginning address of this layer in the 
       * block power vector
       */
      base += model->layers[n].flp->n_units;							 
  }

  /* extra spreader and sink nodes	*/
  for(i=0; i < extra_nodes; i++)
    g->extra[i] = b[base+i];
}

/* translate temperature between grid and block vectors	*/
void xlate_temp_g2b(grid_model_t *model, double *b, grid_model_vector_t *g)
{
  int i, j, n, u, base = 0, count;
  int i1, j1, i2, j2, ci1, cj1, ci2, cj2;
  double min, max, avg;

  int extra_nodes;
  if (model->config.model_secondary)
    extra_nodes = EXTRA + EXTRA_SEC;
  else
    extra_nodes = EXTRA;

  for(n=0; n < model->n_layers; n++) {
      for(u=0; u < model->layers[n].flp->n_units; u++) {
          /* extent of this unit in grid cell units	*/
          i1 = model->layers[n].g2bmap[u].i1;
          j1 = model->layers[n].g2bmap[u].j1;
          i2 = model->layers[n].g2bmap[u].i2;
          j2 = model->layers[n].g2bmap[u].j2;

          /* map the center grid cell's temperature to the block	*/
          if (model->map_mode == GRID_CENTER) {
              /* center co-ordinates	*/	
              ci1 = (i1 + i2) / 2;
              cj1 = (j1 + j2) / 2;
              /* in case of even no. of cells, center 
               * is the average of two central cells
               */
              /* ci2 = ci1-1 when even, ci1 otherwise	*/  
              ci2 = ci1 - !((i2-i1) % 2);
              /* cj2 = cj1-1 when even, cj1 otherwise	*/  
              cj2 = cj1 - !((j2-j1) % 2);

              b[base+u] = (g->cuboid[n][ci1][cj1] + g->cuboid[n][ci2][cj1] + 
                           g->cuboid[n][ci1][cj2] + g->cuboid[n][ci2][cj2]) / 4;
              continue;
          }

          /* find the min/max/avg temperatures of the 
           * grid cells in this block
           */
          avg = 0.0;
          count = 0;
          min = max = g->cuboid[n][i1][j1];
          for(i=i1; i < i2; i++)
            for(j=j1; j < j2; j++) {
                avg += g->cuboid[n][i][j];
                if (g->cuboid[n][i][j] < min)
                  min = g->cuboid[n][i][j];
                if (g->cuboid[n][i][j] > max)
                  max = g->cuboid[n][i][j];
                count++;
            }

          /* map to output accordingly	*/
          switch (model->map_mode)
            {
            case GRID_AVG:
              b[base+u] = avg / count;
              break;
            case GRID_MIN:
              b[base+u] = min;
              break;
            case GRID_MAX:
              b[base+u] = max;
              break;
              /* taken care of already	*/	
            case GRID_CENTER:
              break;
            default:
              fatal("unknown mapping mode\n");
              break;
            }
      }
      /* keep track of the beginning address of this layer in the 
       * block power vector
       */
      base += model->layers[n].flp->n_units;							 
  }

  /* extra spreader and sink nodes	*/
  for(i=0; i < extra_nodes; i++)
    b[base+i] = g->extra[i];
}

/* setting package nodes' power numbers	*/
void set_internal_power_grid(grid_model_t *model, double *power)
{
  if (!model->config.model_secondary)
    zero_dvector(&power[model->total_n_blocks], EXTRA);
  else
    zero_dvector(&power[model->total_n_blocks], EXTRA+EXTRA_SEC);
}

/* set up initial temperatures for the steady state solution
 * heuristically (ignoring the lateral resistances)
 */
void set_heuristic_temp(grid_model_t *model, grid_model_vector_t *power, 
                        grid_model_vector_t *temp)
{
  int n, i, j, nl, nr, nc;
  double **sum;

  /* shortcuts	*/
  nl = model->n_layers;
  nr = model->rows;
  nc = model->cols;

  /* package temperatures	*/
  /* if all lateral resistances are considered infinity, all peripheral 
   * package nodes are at the ambient temperature
   */
  temp->extra[SINK_N] = temp->extra[SINK_S] = 
    temp->extra[SINK_E] = temp->extra[SINK_W] = 
    temp->extra[SINK_C_N] = temp->extra[SINK_C_S] = 
    temp->extra[SINK_C_E] = temp->extra[SINK_C_W] = 
    temp->extra[SP_N] = temp->extra[SP_S] = 
    temp->extra[SP_E] = temp->extra[SP_W] =
    model->config.ambient;

  if (model->config.model_secondary) {
      temp->extra[PCB_N] = temp->extra[PCB_S] = 
        temp->extra[PCB_E] = temp->extra[PCB_W] = 
        temp->extra[PCB_C_N] = temp->extra[PCB_C_S] = 
        temp->extra[PCB_C_E] = temp->extra[PCB_C_W] =
        temp->extra[SOLDER_N] = temp->extra[SOLDER_S] =
        temp->extra[SOLDER_E] = temp->extra[SOLDER_W] = 
        temp->extra[SUB_N] = temp->extra[SUB_S] = 
        temp->extra[SUB_E] = temp->extra[SUB_W] =
        model->config.ambient;
  }

  /* layer temperatures	*/
  /* add up power for each grid cell across all layers */
  sum = dmatrix(nr, nc);
  for(n=0; n < nl; n++)
    scaleadd_dvector(sum[0], sum[0], power->cuboid[n][0], nr*nc, 1.0);

  /* last layer	*/
  for(i=0; i < nr; i++)
    for(j=0; j < nc; j++)
      temp->cuboid[nl-1][i][j] = model->config.ambient + sum[i][j] * 
        model->layers[nl-1].rz;
  /* subtract away the layer's power	*/			
  scaleadd_dvector(sum[0], sum[0], power->cuboid[nl-1][0], nr*nc, -1.0);

  /* go from last-1 to first	*/
  for(n=nl-2; n >= 0; n--) {
      /* nth layer temp is n+1th temp + cumul_power * rz of the nth layer	*/
      scaleadd_dvector(temp->cuboid[n][0], temp->cuboid[n+1][0], sum[0],
                       nr*nc, model->layers[n].rz);
      /* subtract away the layer's power	*/			
      scaleadd_dvector(sum[0], sum[0], power->cuboid[n][0], nr*nc, -1.0);
  } 
  free_dmatrix(sum); 
}

/* single steady state iteration of grid solver - package part */
double single_iteration_steady_pack(grid_model_t *model, grid_model_vector_t *power,
                                    grid_model_vector_t *temp)
{
  int i, j;

  double delta[EXTRA+EXTRA_SEC], max = 0;
  /* sum of the conductances	*/
  double csum;
  /* weighted sum of temperatures	*/
  double wsum;

  /* shortcuts	*/
  double *v = temp->extra;
  package_RC_t *pk = &model->pack;
  thermal_config_t *c = &model->config;
  layer_t *l = model->layers;
  int nl = model->n_layers;
  int nr = model->rows;
  int nc = model->cols;
  int spidx, hsidx;
  int subidx = LAYER_SUB;
  int solderidx = LAYER_SOLDER;
  int pcbidx = LAYER_PCB;
  spidx = nl - DEFAULT_PACK_LAYERS + LAYER_SP;
  hsidx = nl - DEFAULT_PACK_LAYERS + LAYER_SINK;

  /* sink outer north/south	*/
  csum = 1.0/(pk->r_hs_per + pk->r_amb_per) + 1.0/(pk->r_hs2_y + pk->r_hs);
  wsum = c->ambient/(pk->r_hs_per + pk->r_amb_per) + v[SINK_C_N]/(pk->r_hs2_y + pk->r_hs);
  delta[SINK_N] = fabs(v[SINK_N] - wsum / csum);
  v[SINK_N] = wsum / csum;
  wsum = c->ambient/(pk->r_hs_per + pk->r_amb_per) + v[SINK_C_S]/(pk->r_hs2_y + pk->r_hs);
  delta[SINK_S] = fabs(v[SINK_S] - wsum / csum);
  v[SINK_S] = wsum / csum;

  /* sink outer west/east	*/
  csum = 1.0/(pk->r_hs_per + pk->r_amb_per) + 1.0/(pk->r_hs2_x + pk->r_hs);
  wsum = c->ambient/(pk->r_hs_per + pk->r_amb_per) + v[SINK_C_W]/(pk->r_hs2_x + pk->r_hs);
  delta[SINK_W] = fabs(v[SINK_W] - wsum / csum);
  v[SINK_W] = wsum / csum;
  wsum = c->ambient/(pk->r_hs_per + pk->r_amb_per) + v[SINK_C_E]/(pk->r_hs2_x + pk->r_hs);
  delta[SINK_E] = fabs(v[SINK_E] - wsum / csum);
  v[SINK_E] = wsum / csum;

  /* sink inner north/south	*/
  /* partition r_hs1_y among all the nc grid cells. edge cell has half the ry */
  csum = nc / (l[hsidx].ry / 2.0 + nc * pk->r_hs1_y);
  csum += 1.0/(pk->r_hs_c_per_y + pk->r_amb_c_per_y) + 
    1.0/pk->r_sp_per_y + 1.0/(pk->r_hs2_y + pk->r_hs);

  wsum = 0.0;
  for(j=0; j < nc; j++)
    wsum += temp->cuboid[hsidx][0][j];
  wsum /= (l[hsidx].ry / 2.0 + nc * pk->r_hs1_y);
  wsum += c->ambient/(pk->r_hs_c_per_y + pk->r_amb_c_per_y) + 
    v[SP_N]/pk->r_sp_per_y + v[SINK_N]/(pk->r_hs2_y + pk->r_hs);
  delta[SINK_C_N] = fabs(v[SINK_C_N] - wsum / csum);
  v[SINK_C_N] = wsum / csum;

  wsum = 0.0;
  for(j=0; j < nc; j++)
    wsum += temp->cuboid[hsidx][nr-1][j];
  wsum /= (l[hsidx].ry / 2.0 + nc * pk->r_hs1_y);
  wsum += c->ambient/(pk->r_hs_c_per_y + pk->r_amb_c_per_y) + 
    v[SP_S]/pk->r_sp_per_y + v[SINK_S]/(pk->r_hs2_y + pk->r_hs);
  delta[SINK_C_S] = fabs(v[SINK_C_S] - wsum / csum);
  v[SINK_C_S] = wsum / csum;

  /* sink inner west/east	*/
  /* partition r_hs1_x among all the nr grid cells. edge cell has half the rx */
  csum = nr / (l[hsidx].rx / 2.0 + nr * pk->r_hs1_x);
  csum += 1.0/(pk->r_hs_c_per_x + pk->r_amb_c_per_x) + 
    1.0/pk->r_sp_per_x + 1.0/(pk->r_hs2_x + pk->r_hs);

  wsum = 0.0;
  for(i=0; i < nr; i++)
    wsum += temp->cuboid[hsidx][i][0];
  wsum /= (l[hsidx].rx / 2.0 + nr * pk->r_hs1_x);
  wsum += c->ambient/(pk->r_hs_c_per_x + pk->r_amb_c_per_x) + 
    v[SP_W]/pk->r_sp_per_x + v[SINK_W]/(pk->r_hs2_x + pk->r_hs);
  delta[SINK_C_W] = fabs(v[SINK_C_W] - wsum / csum);
  v[SINK_C_W] = wsum / csum;

  wsum = 0.0;
  for(i=0; i < nr; i++)
    wsum += temp->cuboid[hsidx][i][nc-1];
  wsum /= (l[hsidx].rx / 2.0 + nr * pk->r_hs1_x);
  wsum += c->ambient/(pk->r_hs_c_per_x + pk->r_amb_c_per_x) + 
    v[SP_E]/pk->r_sp_per_x + v[SINK_E]/(pk->r_hs2_x + pk->r_hs);
  delta[SINK_C_E] = fabs(v[SINK_C_E] - wsum / csum);
  v[SINK_C_E] = wsum / csum;

  /* spreader north/south	*/
  /* partition r_sp1_y among all the nc grid cells. edge cell has half the ry */
  csum = nc / (l[spidx].ry / 2.0 + nc * pk->r_sp1_y);
  csum += 1.0/pk->r_sp_per_y;

  wsum = 0.0;
  for(j=0; j < nc; j++)
    wsum += temp->cuboid[spidx][0][j];
  wsum /= (l[spidx].ry / 2.0 + nc * pk->r_sp1_y);
  wsum += v[SINK_C_N]/pk->r_sp_per_y;
  delta[SP_N] = fabs(v[SP_N] - wsum / csum);
  v[SP_N] = wsum / csum;

  wsum = 0.0;
  for(j=0; j < nc; j++)
    wsum += temp->cuboid[spidx][nr-1][j];
  wsum /= (l[spidx].ry / 2.0 + nc * pk->r_sp1_y);
  wsum += v[SINK_C_S]/pk->r_sp_per_y;
  delta[SP_S] = fabs(v[SP_S] - wsum / csum);
  v[SP_S] = wsum / csum;

  /* spreader west/east	*/
  /* partition r_sp1_x among all the nr grid cells. edge cell has half the rx */
  csum = nr / (l[spidx].rx / 2.0 + nr * pk->r_sp1_x);
  csum += 1.0/pk->r_sp_per_x;

  wsum = 0.0;
  for(i=0; i < nr; i++)
    wsum += temp->cuboid[spidx][i][0];
  wsum /= (l[spidx].rx / 2.0 + nr * pk->r_sp1_x);
  wsum += v[SINK_C_W]/pk->r_sp_per_x;
  delta[SP_W] = fabs(v[SP_W] - wsum / csum);
  v[SP_W] = wsum / csum;

  wsum = 0.0;
  for(i=0; i < nr; i++)
    wsum += temp->cuboid[spidx][i][nc-1];
  wsum /= (l[spidx].rx / 2.0 + nr * pk->r_sp1_x);
  wsum += v[SINK_C_E]/pk->r_sp_per_x;
  delta[SP_E] = fabs(v[SP_E] - wsum / csum);
  v[SP_E] = wsum / csum;

  if (model->config.model_secondary) {
      /* secondary path package nodes */
      /* PCB outer north/south	*/
      csum = 1.0/(pk->r_amb_sec_per) + 1.0/(pk->r_pcb2_y + pk->r_pcb);
      wsum = c->ambient/(pk->r_amb_sec_per) + v[PCB_C_N]/(pk->r_pcb2_y + pk->r_pcb);
      delta[PCB_N] = fabs(v[PCB_N] - wsum / csum);
      v[PCB_N] = wsum / csum;
      wsum = c->ambient/(pk->r_amb_sec_per) + v[PCB_C_S]/(pk->r_pcb2_y + pk->r_pcb);
      delta[PCB_S] = fabs(v[PCB_S] - wsum / csum);
      v[PCB_S] = wsum / csum;

      /* PCB outer west/east	*/
      csum = 1.0/(pk->r_amb_sec_per) + 1.0/(pk->r_pcb2_x + pk->r_pcb);
      wsum = c->ambient/(pk->r_amb_sec_per) + v[PCB_C_W]/(pk->r_pcb2_x + pk->r_pcb);
      delta[PCB_W] = fabs(v[PCB_W] - wsum / csum);
      v[PCB_W] = wsum / csum;
      wsum = c->ambient/(pk->r_amb_sec_per) + v[PCB_C_E]/(pk->r_pcb2_x + pk->r_pcb);
      delta[PCB_E] = fabs(v[PCB_E] - wsum / csum);
      v[PCB_E] = wsum / csum;

      /* PCB inner north/south	*/
      /* partition r_pcb1_y among all the nc grid cells. edge cell has half the ry */
      csum = nc / (l[pcbidx].ry / 2.0 + nc * pk->r_pcb1_y);
      csum += 1.0/(pk->r_amb_sec_c_per_y) + 
        1.0/pk->r_pcb_c_per_y + 1.0/(pk->r_pcb2_y + pk->r_pcb);

      wsum = 0.0;
      for(j=0; j < nc; j++)
        wsum += temp->cuboid[pcbidx][0][j];
      wsum /= (l[pcbidx].ry / 2.0 + nc * pk->r_pcb1_y);
      wsum += c->ambient/(pk->r_amb_sec_c_per_y) + 
        v[SOLDER_N]/pk->r_pcb_c_per_y + v[PCB_N]/(pk->r_pcb2_y + pk->r_pcb);
      delta[PCB_C_N] = fabs(v[PCB_C_N] - wsum / csum);
      v[PCB_C_N] = wsum / csum;

      wsum = 0.0;
      for(j=0; j < nc; j++)
        wsum += temp->cuboid[pcbidx][nr-1][j];
      wsum /= (l[pcbidx].ry / 2.0 + nc * pk->r_pcb1_y);
      wsum += c->ambient/(pk->r_amb_sec_c_per_y) + 
        v[SOLDER_S]/pk->r_pcb_c_per_y + v[PCB_S]/(pk->r_pcb2_y + pk->r_pcb);
      delta[PCB_C_S] = fabs(v[PCB_C_S] - wsum / csum);
      v[PCB_C_S] = wsum / csum;

      /* PCB inner west/east	*/
      /* partition r_pcb1_x among all the nr grid cells. edge cell has half the rx */
      csum = nr / (l[pcbidx].rx / 2.0 + nr * pk->r_pcb1_x);
      csum += 1.0/(pk->r_amb_sec_c_per_x) + 
        1.0/pk->r_pcb_c_per_x + 1.0/(pk->r_pcb2_x + pk->r_pcb);

      wsum = 0.0;
      for(i=0; i < nr; i++)
        wsum += temp->cuboid[pcbidx][i][0];
      wsum /= (l[pcbidx].rx / 2.0 + nr * pk->r_pcb1_x);
      wsum += c->ambient/(pk->r_amb_sec_c_per_x) + 
        v[SOLDER_W]/pk->r_pcb_c_per_x + v[PCB_W]/(pk->r_pcb2_x + pk->r_pcb);
      delta[PCB_C_W] = fabs(v[PCB_C_W] - wsum / csum);
      v[PCB_C_W] = wsum / csum;

      wsum = 0.0;
      for(i=0; i < nr; i++)
        wsum += temp->cuboid[pcbidx][i][nc-1];
      wsum /= (l[pcbidx].rx / 2.0 + nr * pk->r_pcb1_x);
      wsum += c->ambient/(pk->r_amb_sec_c_per_x) + 
        v[SOLDER_E]/pk->r_pcb_c_per_x + v[PCB_E]/(pk->r_pcb2_x + pk->r_pcb);
      delta[PCB_C_E] = fabs(v[PCB_C_E] - wsum / csum);
      v[PCB_C_E] = wsum / csum;

      /* solder north/south	*/
      /* partition r_solder1_y among all the nc grid cells. edge cell has half the ry */
      csum = nc / (l[solderidx].ry / 2.0 + nc * pk->r_solder1_y);
      csum += 1.0/pk->r_solder_per_y + 1.0/pk->r_pcb_c_per_y;

      wsum = 0.0;
      for(j=0; j < nc; j++)
        wsum += temp->cuboid[solderidx][0][j];
      wsum /= (l[solderidx].ry / 2.0 + nc * pk->r_solder1_y);
      wsum += v[PCB_C_N]/pk->r_pcb_c_per_y + v[SUB_N]/pk->r_solder_per_y;
      delta[SOLDER_N] = fabs(v[SOLDER_N] - wsum / csum);
      v[SOLDER_N] = wsum / csum;

      wsum = 0.0;
      for(j=0; j < nc; j++)
        wsum += temp->cuboid[solderidx][nr-1][j];
      wsum /= (l[solderidx].ry / 2.0 + nc * pk->r_solder1_y);
      wsum += v[PCB_C_S]/pk->r_pcb_c_per_y + v[SUB_S]/pk->r_solder_per_y;
      delta[SOLDER_S] = fabs(v[SOLDER_S] - wsum / csum);
      v[SOLDER_S] = wsum / csum;

      /* solder west/east	*/
      /* partition r_solder1_x among all the nr grid cells. edge cell has half the rx */
      csum = nr / (l[solderidx].rx / 2.0 + nr * pk->r_solder1_x);
      csum += 1.0/pk->r_solder_per_x + 1.0/pk->r_pcb_c_per_x;

      wsum = 0.0;
      for(i=0; i < nr; i++)
        wsum += temp->cuboid[solderidx][i][0];
      wsum /= (l[solderidx].rx / 2.0 + nr * pk->r_solder1_x);
      wsum += v[PCB_C_W]/pk->r_pcb_c_per_x + v[SUB_W]/pk->r_solder_per_x;
      delta[SOLDER_W] = fabs(v[SOLDER_W] - wsum / csum);
      v[SOLDER_W] = wsum / csum;

      wsum = 0.0;
      for(i=0; i < nr; i++)
        wsum += temp->cuboid[solderidx][i][nc-1];
      wsum /= (l[solderidx].rx / 2.0 + nr * pk->r_solder1_x);
      wsum += v[PCB_C_E]/pk->r_pcb_c_per_x + v[SUB_E]/pk->r_solder_per_x;
      delta[SOLDER_E] = fabs(v[SOLDER_E] - wsum / csum);
      v[SOLDER_E] = wsum / csum; 

      /* substrate north/south	*/
      /* partition r_sub1_y among all the nc grid cells. edge cell has half the ry */
      csum = nc / (l[subidx].ry / 2.0 + nc * pk->r_sub1_y);
      csum += 1.0/pk->r_solder_per_y;

      wsum = 0.0;
      for(j=0; j < nc; j++)
        wsum += temp->cuboid[subidx][0][j];
      wsum /= (l[subidx].ry / 2.0 + nc * pk->r_sub1_y);
      wsum += v[SOLDER_N]/pk->r_solder_per_y;
      delta[SUB_N] = fabs(v[SUB_N] - wsum / csum);
      v[SUB_N] = wsum / csum;

      wsum = 0.0;
      for(j=0; j < nc; j++)
        wsum += temp->cuboid[subidx][nr-1][j];
      wsum /= (l[subidx].ry / 2.0 + nc * pk->r_sub1_y);
      wsum += v[SOLDER_S]/pk->r_solder_per_y;
      delta[SUB_S] = fabs(v[SUB_S] - wsum / csum);
      v[SUB_S] = wsum / csum;

      /* substrate west/east	*/
      /* partition r_sub1_x among all the nr grid cells. edge cell has half the rx */
      csum = nr / (l[subidx].rx / 2.0 + nr * pk->r_sub1_x);
      csum += 1.0/pk->r_solder_per_x;

      wsum = 0.0;
      for(i=0; i < nr; i++)
        wsum += temp->cuboid[subidx][i][0];
      wsum /= (l[subidx].rx / 2.0 + nr * pk->r_sub1_x);
      wsum += v[SOLDER_W]/pk->r_solder_per_x;
      delta[SUB_W] = fabs(v[SUB_W] - wsum / csum);
      v[SUB_W] = wsum / csum;

      wsum = 0.0;
      for(i=0; i < nr; i++)
        wsum += temp->cuboid[subidx][i][nc-1];
      wsum /= (l[subidx].rx / 2.0 + nr * pk->r_sub1_x);
      wsum += v[SOLDER_E]/pk->r_solder_per_x;
      delta[SUB_E] = fabs(v[SUB_E] - wsum / csum);
      v[SUB_E] = wsum / csum;
  } 

  if (!model->config.model_secondary) {
      for(i=0; i < EXTRA; i++) {
          if (delta[i] > max)
            max = delta[i];
      }
  } else {
      for(i=0; i < EXTRA + EXTRA_SEC; i++) {
          if (delta[i] > max)
            max = delta[i];
      }
  }
  return max;
}

/* macros for calculating conductances	*/
/* conductance to the next cell north. zero if on northern boundary	*/
# define NC(l,n,i,j,nl,nr,nc)		((i > 0) ? (1.0/l[n].ry) : 0.0)
/* conductance to the next cell south. zero if on southern boundary	*/
# define SC(l,n,i,j,nl,nr,nc)		((i < nr-1) ? (1.0/l[n].ry) : 0.0)
/* conductance to the next cell east. zero if on eastern boundary	*/
# define EC(l,n,i,j,nl,nr,nc)		((j < nc-1) ? (1.0/l[n].rx) : 0.0)
/* conductance to the next cell west. zero if on western boundary	*/
# define WC(l,n,i,j,nl,nr,nc)		((j > 0) ? (1.0/l[n].rx) : 0.0)
/* conductance to the next cell below. zero if on bottom face		*/
# define BC(l,n,i,j,nl,nr,nc)		((n < nl-1) ? (1.0/l[n].rz) : 0.0)
/* conductance to the next cell above. zero if on top face			*/
# define AC(l,n,i,j,nl,nr,nc)		((n > 0) ? (1.0/l[n-1].rz) : 0.0)

/* macros for calculating weighted temperatures	*/
/* weighted T of the next cell north. zero if on northern boundary	*/
# define NT(l,v,n,i,j,nl,nr,nc)		((i > 0) ? (v[n][i-1][j]/l[n].ry) : 0.0)
/* weighted T of the next cell south. zero if on southern boundary	*/
# define ST(l,v,n,i,j,nl,nr,nc)		((i < nr-1) ? (v[n][i+1][j]/l[n].ry) : 0.0)
/* weighted T of the next cell east. zero if on eastern boundary	*/
# define ET(l,v,n,i,j,nl,nr,nc)		((j < nc-1) ? (v[n][i][j+1]/l[n].rx) : 0.0)
/* weighted T of the next cell west. zero if on western boundary	*/
# define WT(l,v,n,i,j,nl,nr,nc)		((j > 0) ? (v[n][i][j-1]/l[n].rx) : 0.0)
/* weighted T of the next cell below. zero if on bottom face		*/
# define BT(l,v,n,i,j,nl,nr,nc)		((n < nl-1) ? (v[n+1][i][j]/l[n].rz) : 0.0)
/* weighted T of the next cell above. zero if on top face			*/
# define AT(l,v,n,i,j,nl,nr,nc)		((n > 0) ? (v[n-1][i][j]/l[n-1].rz) : 0.0)

/* single steady state iteration of grid solver - silicon part */
double single_iteration_steady_grid(grid_model_t *model, grid_model_vector_t *power,
                                    grid_model_vector_t *temp)
{
  int n, i, j;
  double prev, delta, max = 0;
  /* sum of the conductances	*/
  double csum;
  /* weighted sum of temperatures	*/
  double wsum;
  double north, south, east, west, above, below;

  /* shortcuts for cell width(cw) and cell height(ch)	*/
  double cw = model->width / model->cols;
  double ch = model->height / model->rows;

  /* shortcuts	*/
  double ***v = temp->cuboid;
  thermal_config_t *c = &model->config;
  layer_t *l = model->layers;
  int nl = model->n_layers;
  int nr = model->rows;
  int nc = model->cols;
  int spidx, hsidx;
  int subidx = LAYER_SUB;
  int solderidx = LAYER_SOLDER;
  int pcbidx = LAYER_PCB;
  spidx = nl - DEFAULT_PACK_LAYERS + LAYER_SP;
  hsidx = nl - DEFAULT_PACK_LAYERS + LAYER_SINK;

  /* for each grid cell	*/
  for(n=0; n < nl; n++) {
      for(i=0; i < nr; i++) {
          for(j=0; j < nc; j++) {
              /* sum the conductances to cells north, south, 
               * east, west, above and below
               */
              // BU_3D: call new macros if detailed_3D model is used 
              // the spreader/heat sink layers will use uniform R
              if(model->config.detailed_3D_used){
#if GATE_GS1_ARM == GS1_ARM_GS_NATIVE
                  /* Preserve the shipped detailed-3D GS coefficient lookup. */
                  north = (i > 0) ?
                    1.0 / find_res_3D(n, i-1, j, model, 2) : 0.0;
                  south = (i < nr-1) ?
                    1.0 / find_res_3D(n, i+1, j, model, 2) : 0.0;
                  east = (j < nc-1) ?
                    1.0 / find_res_3D(n, i, j+1, model, 1) : 0.0;
                  west = (j > 0) ?
                    1.0 / find_res_3D(n, i, j-1, model, 1) : 0.0;
                  above = (n > 0) ?
                    1.0 / find_res_3D(n-1, i, j, model, 3) : 0.0;
                  below = (n < nl-1) ?
                    1.0 / find_res_3D(n, i, j, model, 3) : 0.0;
#else
                  north = STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_NORTH);
                  south = STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_SOUTH);
                  east = STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_EAST);
                  west = STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_WEST);
                  above = STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_ABOVE);
                  below = STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_BELOW);
#endif

                  csum = north + south + east + west + above + below;

                  /* sum of the weighted temperatures of all the neighbours */
                  wsum = ((i > 0) ? v[n][i-1][j] * north : 0.0) +
                    ((i < nr-1) ? v[n][i+1][j] * south : 0.0) +
                    ((j < nc-1) ? v[n][i][j+1] * east : 0.0) +
                    ((j > 0) ? v[n][i][j-1] * west : 0.0) +
                    ((n > 0) ? v[n-1][i][j] * above : 0.0) +
                    ((n < nl-1) ? v[n+1][i][j] * below : 0.0);
              } //end->BU_3D
              else {	
                  csum = NC(l,n,i,j,nl,nr,nc) + SC(l,n,i,j,nl,nr,nc) + 
                    EC(l,n,i,j,nl,nr,nc) + WC(l,n,i,j,nl,nr,nc) + 
                    AC(l,n,i,j,nl,nr,nc) + BC(l,n,i,j,nl,nr,nc);

                  /* sum of the weighted temperatures of all the neighbours	*/
                  wsum = NT(l,v,n,i,j,nl,nr,nc) + ST(l,v,n,i,j,nl,nr,nc) + 
                    ET(l,v,n,i,j,nl,nr,nc) + WT(l,v,n,i,j,nl,nr,nc) + 
                    AT(l,v,n,i,j,nl,nr,nc) + BT(l,v,n,i,j,nl,nr,nc);
              } 

              /* spreader core is connected to its periphery	*/
              if (n == spidx) {
                  /* northern boundary - edge cell has half the ry	*/
                  if (i == 0) {
                      csum += 1.0/(l[n].ry/2.0 + nc*model->pack.r_sp1_y); 
                      wsum += temp->extra[SP_N]/(l[n].ry/2.0 + nc*model->pack.r_sp1_y); 
                  }
                  /* southern boundary - edge cell has half the ry	*/
                  if (i == nr-1) {
                      csum += 1.0/(l[n].ry/2.0 + nc*model->pack.r_sp1_y); 
                      wsum += temp->extra[SP_S]/(l[n].ry/2.0 + nc*model->pack.r_sp1_y); 
                  }
                  /* eastern boundary	 - edge cell has half the rx	*/
                  if (j == nc-1) {
                      csum += 1.0/(l[n].rx/2.0 + nr*model->pack.r_sp1_x); 
                      wsum += temp->extra[SP_E]/(l[n].rx/2.0 + nr*model->pack.r_sp1_x); 
                  }
                  /* western boundary	- edge cell has half the rx		*/
                  if (j == 0) {
                      csum += 1.0/(l[n].rx/2.0 + nr*model->pack.r_sp1_x); 
                      wsum += temp->extra[SP_W]/(l[n].rx/2.0 + nr*model->pack.r_sp1_x); 
                  }
                  /* heatsink core is connected to its inner periphery and ambient	*/
              } else if (n == hsidx) {
                  /* all nodes are connected to the ambient	*/
                  csum += 1.0/l[n].rz;
                  wsum += c->ambient/l[n].rz;
                  /* northern boundary - edge cell has half the ry	*/
                  if (i == 0) {
                      csum += 1.0/(l[n].ry/2.0 + nc*model->pack.r_hs1_y); 
                      wsum += temp->extra[SINK_C_N]/(l[n].ry/2.0 + nc*model->pack.r_hs1_y); 
                  }
                  /* southern boundary - edge cell has half the ry	*/
                  if (i == nr-1) {
                      csum += 1.0/(l[n].ry/2.0 + nc*model->pack.r_hs1_y); 
                      wsum += temp->extra[SINK_C_S]/(l[n].ry/2.0 + nc*model->pack.r_hs1_y); 
                  }
                  /* eastern boundary	 - edge cell has half the rx	*/
                  if (j == nc-1) {
                      csum += 1.0/(l[n].rx/2.0 + nr*model->pack.r_hs1_x); 
                      wsum += temp->extra[SINK_C_E]/(l[n].rx/2.0 + nr*model->pack.r_hs1_x); 
                  }
                  /* western boundary	- edge cell has half the rx		*/
                  if (j == 0) {
                      csum += 1.0/(l[n].rx/2.0 + nr*model->pack.r_hs1_x); 
                      wsum += temp->extra[SINK_C_W]/(l[n].rx/2.0 + nr*model->pack.r_hs1_x); 
                  }
              } else if (model->config.model_secondary && n==subidx) {
                  /* northern boundary - edge cell has half the ry	*/
                  if (i == 0) {
                      csum += 1.0/(l[n].ry/2.0 + nc*model->pack.r_sub1_y); 
                      wsum += temp->extra[SUB_N]/(l[n].ry/2.0 + nc*model->pack.r_sub1_y); 
                  }
                  /* southern boundary - edge cell has half the ry	*/
                  if (i == nr-1) {
                      csum += 1.0/(l[n].ry/2.0 + nc*model->pack.r_sub1_y); 
                      wsum += temp->extra[SUB_S]/(l[n].ry/2.0 + nc*model->pack.r_sub1_y); 
                  }
                  /* eastern boundary	 - edge cell has half the rx	*/
                  if (j == nc-1) {
                      csum += 1.0/(l[n].rx/2.0 + nr*model->pack.r_sub1_x); 
                      wsum += temp->extra[SUB_E]/(l[n].rx/2.0 + nr*model->pack.r_sub1_x); 
                  }
                  /* western boundary	- edge cell has half the rx		*/
                  if (j == 0) {
                      csum += 1.0/(l[n].rx/2.0 + nr*model->pack.r_sub1_x); 
                      wsum += temp->extra[SUB_W]/(l[n].rx/2.0 + nr*model->pack.r_sub1_x); 
                  } 
              } else if (model->config.model_secondary && n==solderidx) {
                  /* northern boundary - edge cell has half the ry	*/
                  if (i == 0) {
                      csum += 1.0/(l[n].ry/2.0 + nc*model->pack.r_solder1_y); 
                      wsum += temp->extra[SOLDER_N]/(l[n].ry/2.0 + nc*model->pack.r_solder1_y); 
                  }
                  /* southern boundary - edge cell has half the ry	*/
                  if (i == nr-1) {
                      csum += 1.0/(l[n].ry/2.0 + nc*model->pack.r_solder1_y); 
                      wsum += temp->extra[SOLDER_S]/(l[n].ry/2.0 + nc*model->pack.r_solder1_y); 
                  }
                  /* eastern boundary	 - edge cell has half the rx	*/
                  if (j == nc-1) {
                      csum += 1.0/(l[n].rx/2.0 + nr*model->pack.r_solder1_x); 
                      wsum += temp->extra[SOLDER_E]/(l[n].rx/2.0 + nr*model->pack.r_solder1_x); 
                  }
                  /* western boundary	- edge cell has half the rx		*/
                  if (j == 0) {
                      csum += 1.0/(l[n].rx/2.0 + nr*model->pack.r_solder1_x); 
                      wsum += temp->extra[SOLDER_W]/(l[n].rx/2.0 + nr*model->pack.r_solder1_x); 
                  } 
              } else if (model->config.model_secondary && n==pcbidx) {
                  /* all nodes are connected to the ambient	*/
                  csum += 1.0/(model->config.r_convec_sec * 
                               (model->config.s_pcb * model->config.s_pcb) / (cw * ch));
                  wsum += c->ambient/(model->config.r_convec_sec * 
                                      (model->config.s_pcb * model->config.s_pcb) / (cw * ch));
                  /* northern boundary - edge cell has half the ry	*/
                  if (i == 0) {
                      csum += 1.0/(l[n].ry/2.0 + nc*model->pack.r_pcb1_y); 
                      wsum += temp->extra[PCB_C_N]/(l[n].ry/2.0 + nc*model->pack.r_pcb1_y); 
                  }
                  /* southern boundary - edge cell has half the ry	*/
                  if (i == nr-1) {
                      csum += 1.0/(l[n].ry/2.0 + nc*model->pack.r_pcb1_y); 
                      wsum += temp->extra[PCB_C_S]/(l[n].ry/2.0 + nc*model->pack.r_pcb1_y); 
                  }
                  /* eastern boundary	 - edge cell has half the rx	*/
                  if (j == nc-1) {
                      csum += 1.0/(l[n].rx/2.0 + nr*model->pack.r_pcb1_x); 
                      wsum += temp->extra[PCB_C_E]/(l[n].rx/2.0 + nr*model->pack.r_pcb1_x); 
                  }
                  /* western boundary	- edge cell has half the rx		*/
                  if (j == 0) {
                      csum += 1.0/(l[n].rx/2.0 + nr*model->pack.r_pcb1_x); 
                      wsum += temp->extra[PCB_C_W]/(l[n].rx/2.0 + nr*model->pack.r_pcb1_x); 
                  }
              }

              /* update the current cell's temperature	*/	   
              prev = v[n][i][j];
              v[n][i][j] = (power->cuboid[n][i][j] + wsum) / csum;

              /* compute maximum delta	*/
              delta =  fabs(prev - v[n][i][j]);
              if (delta > max)
                max = delta;
          }
      }
  }
  /* package part of the iteration	*/
  return (MAX(max, single_iteration_steady_pack(model, power, temp)));
}

/* restriction operator for multigrid solver. given a power vector
 * corresponding to a fine grid, outputs a vector corresponding
 * to one level coarser grid (half the no. of rows and cols)
 * NOTE: model->rows and model->cols denote the size of the
 * coarser grid
 */
void multigrid_restrict_power(grid_model_t *model, grid_model_vector_t *dst, 
                              grid_model_vector_t *src)
{
  /* coarse grid indices	*/
  int n, i, j;

  /* grid cells - add the four nearest neighbours	*/
  for(n=0; n < model->n_layers; n++)
    for(i=0; i < model->rows; i++)
      for(j=0; j < model->cols; j++)
        {
          dst->cuboid[n][i][j] = (src->cuboid[n][2*i][2*j] +
                                  src->cuboid[n][2*i+1][2*j] +
                                  src->cuboid[n][2*i][2*j+1] +
                                  src->cuboid[n][2*i+1][2*j+1]);
        }
  /* package nodes - copy them as it is	*/
  if (!model->config.model_secondary)
    copy_dvector(dst->extra, src->extra, EXTRA);
  else
    copy_dvector(dst->extra, src->extra, EXTRA+EXTRA_SEC);
}

/* prolongation(interpolation) operator for multigrid solver. 
 * given a temperature vector corresponding to a coarse grid,
 * outputs a (bi)linearly interpolated vector corresponding 
 * to one level finer grid (twice the no. of rows and cols)
 * NOTE: model->rows and model->cols denote the size of the
 * coarser grid
 */
void multigrid_prolong_temp(grid_model_t *model, grid_model_vector_t *dst, 
                            grid_model_vector_t *src)
{
  /* coarse grid indices	*/
  int n, i, j;

  /* shortcuts	*/
  int nr = model->rows;
  int nc = model->cols;
  double ***d = dst->cuboid;
  double ***s = src->cuboid;

  /* For the fine grid cells not on the boundary,
   * we want to linearly interpolate in the region 
   * surrounded by the coarse grid cells (i,j),
   * (i+1,j) (i,j+1) and (i+1, j+1). The fine
   * grid cells in that region are (2i+1, 2j+1),
   * (2i+2, 2j+1), (2i+1, 2j+2) and (2i+2, 2j+2).
   * To interpolate bilinearly, we first interpolate
   * along one axis and then along the other. In
   * each axis, due to the proximity of the fine
   * grid cell co-ordinates to one of the coarse
   * grid cells, the weights are 3/4 and 1/4 (not
   * 1/2 and 1/2) So, repeating them along the  
   * other axis also results in weights of 9/16,
   * 3/16, 3/16 and 1/16
   */ 
  for(n=0; n < model->n_layers; n++)
    for(i=0; i < nr-1; i++)
      for(j=0; j < nc-1; j++) {
          d[n][2*i+1][2*j+1] = 9.0/16.0 * s[n][i][j] +
            3.0/16.0 * s[n][i+1][j] +
            1.0/16.0 * s[n][i+1][j+1] +
            3.0/16.0 * s[n][i][j+1];
          d[n][2*i+2][2*j+1] = 3.0/16.0 * s[n][i][j] +
            9.0/16.0 * s[n][i+1][j] +
            3.0/16.0 * s[n][i+1][j+1] +
            1.0/16.0 * s[n][i][j+1];
          d[n][2*i+2][2*j+2] = 1.0/16.0 * s[n][i][j] +
            3.0/16.0 * s[n][i+1][j] +
            9.0/16.0 * s[n][i+1][j+1] +
            3.0/16.0 * s[n][i][j+1];
          d[n][2*i+1][2*j+2] = 3.0/16.0 * s[n][i][j] +
            1.0/16.0 * s[n][i+1][j] +
            3.0/16.0 * s[n][i+1][j+1] +
            9.0/16.0 * s[n][i][j+1];
      }

  /* for the cells on the boundary, we perform
   * a zeroth order interpolation. i.e., copy
   * the nearest coarse grid cell's value
   * as it is
   */
  for(n=0; n < model->n_layers; n++) {
      for(i=0; i < nr; i++) {
          d[n][2*i][0] = d[n][2*i+1][0] = s[n][i][0];
          d[n][2*i][2*nc-1] = d[n][2*i+1][2*nc-1] = s[n][i][nc-1];
      }	
      for(j=0; j < nc; j++) {
          d[n][0][2*j] = d[n][0][2*j+1] = s[n][0][j];
          d[n][2*nr-1][2*j] = d[n][2*nr-1][2*j+1] = s[n][nr-1][j];
      }
  }

  /* package nodes - copy them as it is	*/
  if (!model->config.model_secondary) 
    copy_dvector(dst->extra, src->extra, EXTRA);
  else
    copy_dvector(dst->extra, src->extra, EXTRA+EXTRA_SEC);
}

/* recursive multigrid solver. it uses the Gauss-Seidel (GS) iterative 
 * solver to solve at a particular grid granularity. Although GS removes
 * high frequency errors in a solution estimate pretty quickly, it 
 * takes a lot of time to eliminate the low frequency ones. These 
 * low frequency errors can be eliminated easily by coarsifying the
 * grid. This is the core principle of mulrigrid solvers. The version here
 * is called "nested iteration". It solves the problem (iteratively using GS)
 * first in the coarsest granularity and then utilizes that solution to 
 * estimate the solution in the next finer grid. This is repeated until
 * the solution is found in the finest desired grid. For details, take 
 * a look at Numerical Recipes in C (2nd edition), sections 19.5 and 19.6
 * (http://www.nrbook.com/a/bookcpdf/c19-5.pdf and
 * http://www.nrbook.com/a/bookcpdf/c19-6.pdf). Also refer to Prof.
 * Jayathi Murthy's ME 608 notes from Purdue - Chapter 8, sections 8.7-8.9.
 * (http://meweb.ecn.purdue.edu/~jmurthy/me608/main.pdf - pg. 175-192)
 */ 
void recursive_multigrid(grid_model_t *model, grid_model_vector_t *power,
                         grid_model_vector_t *temp)
{
  double delta;
#if VERBOSE > 1
  unsigned int i = 0;
#endif
  grid_model_vector_t *coarse_power, *coarse_temp;
  int n;
  /* setup heuristic initial temperatures	at the coarsest level*/
  if (model->rows <= 1 || model->cols <= 1) {
      set_heuristic_temp(model, power, temp);

      /* for finer grids. use coarser solutions as estimates	*/
  } else {
      /* make the grid coarser	*/
      model->rows /= 2;
      model->cols /= 2;
      for(n=0; n < model->n_layers; n++) {
          /* only rz's and c's change. rx's and 
           * ry's remain the same	
           */
          model->layers[n].rz /= 4;
          if (model->c_ready)
            model->layers[n].c *= 4;
      }

      /* vectors for the coarse grid	*/
      coarse_power = new_grid_model_vector(model);
      coarse_temp = new_grid_model_vector(model);

      /* coarsen the power vector	*/
      multigrid_restrict_power(model, coarse_power, power);

      /* solve recursively	*/
      recursive_multigrid(model, coarse_power, coarse_temp);

      /* interpolate the solution to the current fine grid	*/
      multigrid_prolong_temp(model, temp, coarse_temp);

      /* cleanup	*/
      free_grid_model_vector(coarse_power);
      free_grid_model_vector(coarse_temp);

      /* restore the grid */
      model->rows *= 2;
      model->cols *= 2;
      for(n=0; n < model->n_layers; n++) {
          model->layers[n].rz *= 4;
          if (model->c_ready)
            model->layers[n].c /= 4;
      }
  }
  /* refine solution iteratively till convergence	*/
  do {
      delta = single_iteration_steady_grid(model, power, temp);
#if VERBOSE > 1
      i++;
#endif
  } while (!eq(delta, 0));
#if VERBOSE > 1
  fprintf(stdout, "no. of iterations for steady state convergence (%d x %d grid): %d\n", 
          model->rows, model->cols, i);
#endif
}

void steady_state_temp_grid(grid_model_t *model, double *power, double *temp)
{
  grid_model_vector_t *p;
  double total;
#if GATE_GS1_ARM > 0
  double gs1_started_s;
#endif

  if (!model->r_ready)
    fatal("R model not ready\n");

  p = new_grid_model_vector(model);

  /* package nodes' power numbers	*/
  set_internal_power_grid(model, power);

  /* total power - package nodes have no power dissipation	*/
  total = sum_dvector(power, model->total_n_blocks);

  /* map the block power numbers to the grid	*/
  xlate_vector_b2g(model, power, p, V_POWER);

#if GATE_GS1_ARM > 0
  if (!(model->gs1_process_start_s > 0.0))
    fatal("G-S1 process-start clock was not initialized\n");
  model->gs1_input_model_s = gs1_now() - model->gs1_process_start_s -
    model->gs1_coefficient_s;
  if (model->gs1_input_model_s < 0.0)
    fatal("invalid G-S1 input/model phase time\n");
  fprintf(stdout, "G-S1 identity: arm=%s gate=%d superlu=%d\n",
          gs1_arm_name(), GATE_GS1_ARM, SUPERLU);
#endif

#if SUPERLU > 0
#if GATE_GS1_ARM == GS1_ARM_CSR_PCG || \
    GATE_GS1_ARM == GS1_ARM_STRUCTURED_PCG
  if(model->config.detailed_3D_used){
      jacobi_pcg_steady_grid(model, p, model->last_steady);
  }
  else{
      fatal("G-S1.5 PCG representations require detailed-3D mode\n");
  }
#else
  /* solve with SuperLU. use grid model's internal 
   * state vector to store the grid temperatures
   */ 
  direct_SLU(model, p, model->last_steady);
#endif
#else
  /* solve recursively. use grid model's internal 
   * state vector to store the grid temperatures
   */ 
  if(model->config.detailed_3D_used){
#if GATE_GS1_ARM == GS1_ARM_GS_NATIVE
      gs1_gauss_seidel_steady_grid(model, p, model->last_steady);
#else
      jacobi_pcg_steady_grid(model, p, model->last_steady);
#endif
  }
  else{
      recursive_multigrid(model, p, model->last_steady);
  }
#endif

  /* map the temperature numbers back	*/
#if GATE_GS1_ARM > 0
  gs1_started_s = gs1_now();
#endif
  xlate_temp_g2b(model, temp, model->last_steady);
#if GATE_GS1_ARM > 0
  gs1_write_state(model, temp, gs1_output_prefix());
  model->gs1_output_s += gs1_now() - gs1_started_s;
  gs1_write_metadata(model, gs1_output_prefix());
#endif

  free_grid_model_vector(p);
}

/* function to access a 1-d array as a 3-d matrix	*/
#define A3D(array,n,i,j,nl,nr,nc)		(array[(n)*(nr)*(nc) + (i)*(nc) + (j)])

/* compute the slope vector for the package nodes	*/
void slope_fn_pack(grid_model_t *model, const double *v,
                   grid_model_vector_t *p, double *dv)
{
  int i, j;
  /* sum of the currents(power values)	*/
  double psum;

  /* shortcuts	*/
  package_RC_t *pk = &model->pack;
  thermal_config_t *c = &model->config;
  layer_t *l = model->layers;
  int nl = model->n_layers;
  int nr = model->rows;
  int nc = model->cols;
  int spidx, hsidx, subidx, solderidx, pcbidx;
  int model_secondary = model->config.model_secondary;

  /* pointer to the starting address of the extra nodes	*/
  const double *x = v + nl*nr*nc;

  spidx = nl - DEFAULT_PACK_LAYERS + LAYER_SP;
  hsidx = nl - DEFAULT_PACK_LAYERS + LAYER_SINK;
  if (model_secondary) {
      subidx = LAYER_SUB;
      solderidx = LAYER_SOLDER;
      pcbidx = LAYER_PCB;	
  }

  /* sink outer north/south	*/
  psum = (c->ambient - x[SINK_N])/(pk->r_hs_per + pk->r_amb_per) + 
    (x[SINK_C_N] - x[SINK_N])/(pk->r_hs2_y + pk->r_hs);
  dv[nl*nr*nc + SINK_N] = psum / (pk->c_hs_per + pk->c_amb_per);
  psum = (c->ambient - x[SINK_S])/(pk->r_hs_per + pk->r_amb_per) + 
    (x[SINK_C_S] - x[SINK_S])/(pk->r_hs2_y + pk->r_hs);
  dv[nl*nr*nc + SINK_S] = psum / (pk->c_hs_per + pk->c_amb_per);

  /* sink outer west/east	*/
  psum = (c->ambient - x[SINK_W])/(pk->r_hs_per + pk->r_amb_per) + 
    (x[SINK_C_W] - x[SINK_W])/(pk->r_hs2_x + pk->r_hs);
  dv[nl*nr*nc + SINK_W] = psum / (pk->c_hs_per + pk->c_amb_per);
  psum = (c->ambient - x[SINK_E])/(pk->r_hs_per + pk->r_amb_per) + 
    (x[SINK_C_E] - x[SINK_E])/(pk->r_hs2_x + pk->r_hs);
  dv[nl*nr*nc + SINK_E] = psum / (pk->c_hs_per + pk->c_amb_per);

  /* sink inner north/south	*/
  /* partition r_hs1_y among all the nc grid cells. edge cell has half the ry	*/
  psum = 0.0;
  for(j=0; j < nc; j++)
    psum += (A3D(v,hsidx,0,j,nl,nr,nc) - x[SINK_C_N]);
  psum /= (l[hsidx].ry / 2.0 + nc * pk->r_hs1_y);
  psum += (c->ambient - x[SINK_C_N])/(pk->r_hs_c_per_y + pk->r_amb_c_per_y) + 
    (x[SP_N] - x[SINK_C_N])/pk->r_sp_per_y +
    (x[SINK_N] - x[SINK_C_N])/(pk->r_hs2_y + pk->r_hs);
  dv[nl*nr*nc + SINK_C_N] = psum / (pk->c_hs_c_per_y + pk->c_amb_c_per_y);

  psum = 0.0;
  for(j=0; j < nc; j++)
    psum += (A3D(v,hsidx,nr-1,j,nl,nr,nc) - x[SINK_C_S]);
  psum /= (l[hsidx].ry / 2.0 + nc * pk->r_hs1_y);
  psum += (c->ambient - x[SINK_C_S])/(pk->r_hs_c_per_y + pk->r_amb_c_per_y) + 
    (x[SP_S] - x[SINK_C_S])/pk->r_sp_per_y +
    (x[SINK_S] - x[SINK_C_S])/(pk->r_hs2_y + pk->r_hs);
  dv[nl*nr*nc + SINK_C_S] = psum / (pk->c_hs_c_per_y + pk->c_amb_c_per_y);

  /* sink inner west/east	*/
  /* partition r_hs1_x among all the nr grid cells. edge cell has half the rx	*/
  psum = 0.0;
  for(i=0; i < nr; i++)
    psum += (A3D(v,hsidx,i,0,nl,nr,nc) - x[SINK_C_W]);
  psum /= (l[hsidx].rx / 2.0 + nr * pk->r_hs1_x);
  psum += (c->ambient - x[SINK_C_W])/(pk->r_hs_c_per_x + pk->r_amb_c_per_x) + 
    (x[SP_W] - x[SINK_C_W])/pk->r_sp_per_x +
    (x[SINK_W] - x[SINK_C_W])/(pk->r_hs2_x + pk->r_hs);
  dv[nl*nr*nc + SINK_C_W] = psum / (pk->c_hs_c_per_x + pk->c_amb_c_per_x);

  psum = 0.0;
  for(i=0; i < nr; i++)
    psum += (A3D(v,hsidx,i,nc-1,nl,nr,nc) - x[SINK_C_E]);
  psum /= (l[hsidx].rx / 2.0 + nr * pk->r_hs1_x);
  psum += (c->ambient - x[SINK_C_E])/(pk->r_hs_c_per_x + pk->r_amb_c_per_x) + 
    (x[SP_E] - x[SINK_C_E])/pk->r_sp_per_x +
    (x[SINK_E] - x[SINK_C_E])/(pk->r_hs2_x + pk->r_hs);
  dv[nl*nr*nc + SINK_C_E] = psum / (pk->c_hs_c_per_x + pk->c_amb_c_per_x);

  /* spreader north/south	*/
  /* partition r_sp1_y among all the nc grid cells. edge cell has half the ry	*/
  psum = 0.0;
  for(j=0; j < nc; j++)
    psum += (A3D(v,spidx,0,j,nl,nr,nc) - x[SP_N]);
  psum /= (l[spidx].ry / 2.0 + nc * pk->r_sp1_y);
  psum += (x[SINK_C_N] - x[SP_N])/pk->r_sp_per_y;
  dv[nl*nr*nc + SP_N] = psum / pk->c_sp_per_y;

  psum = 0.0;
  for(j=0; j < nc; j++)
    psum += (A3D(v,spidx,nr-1,j,nl,nr,nc) - x[SP_S]);
  psum /= (l[spidx].ry / 2.0 + nc * pk->r_sp1_y);
  psum += (x[SINK_C_S] - x[SP_S])/pk->r_sp_per_y;
  dv[nl*nr*nc + SP_S] = psum / pk->c_sp_per_y;

  /* spreader west/east	*/
  /* partition r_sp1_x among all the nr grid cells. edge cell has half the rx	*/
  psum = 0.0;
  for(i=0; i < nr; i++)
    psum += (A3D(v,spidx,i,0,nl,nr,nc) - x[SP_W]);
  psum /= (l[spidx].rx / 2.0 + nr * pk->r_sp1_x);
  psum += (x[SINK_C_W] - x[SP_W])/pk->r_sp_per_x;
  dv[nl*nr*nc + SP_W] = psum / pk->c_sp_per_x;

  psum = 0.0;
  for(i=0; i < nr; i++)
    psum += (A3D(v,spidx,i,nc-1,nl,nr,nc) - x[SP_E]);
  psum /= (l[spidx].rx / 2.0 + nr * pk->r_sp1_x);
  psum += (x[SINK_C_E] - x[SP_E])/pk->r_sp_per_x;
  dv[nl*nr*nc + SP_E] = psum / pk->c_sp_per_x;

  if (model_secondary) {
      /* PCB outer north/south	*/
      psum = (c->ambient - x[PCB_N])/(pk->r_amb_sec_per) + 
        (x[PCB_C_N] - x[PCB_N])/(pk->r_pcb2_y + pk->r_pcb);
      dv[nl*nr*nc + PCB_N] = psum / (pk->c_pcb_per + pk->c_amb_sec_per);
      psum = (c->ambient - x[PCB_S])/(pk->r_amb_sec_per) + 
        (x[PCB_C_S] - x[PCB_S])/(pk->r_pcb2_y + pk->r_pcb);
      dv[nl*nr*nc + PCB_S] = psum / (pk->c_pcb_per + pk->c_amb_sec_per);

      /* PCB outer west/east	*/
      psum = (c->ambient - x[PCB_W])/(pk->r_amb_sec_per) + 
        (x[PCB_C_W] - x[PCB_W])/(pk->r_pcb2_x + pk->r_pcb);
      dv[nl*nr*nc + PCB_W] = psum / (pk->c_pcb_per + pk->c_amb_sec_per);
      psum = (c->ambient - x[PCB_E])/(pk->r_amb_sec_per) + 
        (x[PCB_C_E] - x[PCB_E])/(pk->r_pcb2_x + pk->r_pcb);
      dv[nl*nr*nc + PCB_E] = psum / (pk->c_pcb_per + pk->c_amb_sec_per);

      /* PCB inner north/south	*/
      /* partition r_pcb1_y among all the nc grid cells. edge cell has half the ry	*/
      psum = 0.0;
      for(j=0; j < nc; j++)
        psum += (A3D(v,pcbidx,0,j,nl,nr,nc) - x[PCB_C_N]);
      psum /= (l[pcbidx].ry / 2.0 + nc * pk->r_pcb1_y);
      psum += (c->ambient - x[PCB_C_N])/(pk->r_amb_sec_c_per_y) + 
        (x[SOLDER_N] - x[PCB_C_N])/pk->r_pcb_c_per_y +
        (x[PCB_N] - x[PCB_C_N])/(pk->r_pcb2_y + pk->r_pcb);
      dv[nl*nr*nc + PCB_C_N] = psum / (pk->c_pcb_c_per_y + pk->c_amb_sec_c_per_y);

      psum = 0.0;
      for(j=0; j < nc; j++)
        psum += (A3D(v,pcbidx,nr-1,j,nl,nr,nc) - x[PCB_C_S]);
      psum /= (l[pcbidx].ry / 2.0 + nc * pk->r_pcb1_y);
      psum += (c->ambient - x[PCB_C_S])/(pk->r_amb_sec_c_per_y) + 
        (x[SOLDER_S] - x[PCB_C_S])/pk->r_pcb_c_per_y +
        (x[PCB_S] - x[PCB_C_S])/(pk->r_pcb2_y + pk->r_pcb);
      dv[nl*nr*nc + PCB_C_S] = psum / (pk->c_pcb_c_per_y + pk->c_amb_sec_c_per_y);

      /* PCB inner west/east	*/
      /* partition r_pcb1_x among all the nr grid cells. edge cell has half the rx	*/
      psum = 0.0;
      for(i=0; i < nr; i++)
        psum += (A3D(v,pcbidx,i,0,nl,nr,nc) - x[PCB_C_W]);
      psum /= (l[pcbidx].rx / 2.0 + nr * pk->r_pcb1_x);
      psum += (c->ambient - x[PCB_C_W])/(pk->r_amb_sec_c_per_x) + 
        (x[SOLDER_W] - x[PCB_C_W])/pk->r_pcb_c_per_x +
        (x[PCB_W] - x[PCB_C_W])/(pk->r_pcb2_x + pk->r_pcb);
      dv[nl*nr*nc + PCB_C_W] = psum / (pk->c_pcb_c_per_x + pk->c_amb_sec_c_per_x);

      psum = 0.0;
      for(i=0; i < nr; i++)
        psum += (A3D(v,pcbidx,i,nc-1,nl,nr,nc) - x[PCB_C_E]);
      psum /= (l[pcbidx].rx / 2.0 + nr * pk->r_pcb1_x);
      psum += (c->ambient - x[PCB_C_E])/(pk->r_amb_sec_c_per_x) + 
        (x[SOLDER_E] - x[PCB_C_E])/pk->r_pcb_c_per_x +
        (x[PCB_E] - x[PCB_C_E])/(pk->r_pcb2_x + pk->r_pcb);
      dv[nl*nr*nc + PCB_C_E] = psum / (pk->c_pcb_c_per_x + pk->c_amb_sec_c_per_x);

      /* solder ball north/south	*/
      /* partition r_solder1_y among all the nc grid cells. edge cell has half the ry	*/
      psum = 0.0;
      for(j=0; j < nc; j++)
        psum += (A3D(v,solderidx,0,j,nl,nr,nc) - x[SOLDER_N]);
      psum /= (l[solderidx].ry / 2.0 + nc * pk->r_solder1_y);
      psum += (x[PCB_C_N] - x[SOLDER_N])/pk->r_pcb_c_per_y;
      dv[nl*nr*nc + SOLDER_N] = psum / pk->c_solder_per_y;

      psum = 0.0;
      for(j=0; j < nc; j++)
        psum += (A3D(v,solderidx,nr-1,j,nl,nr,nc) - x[SOLDER_S]);
      psum /= (l[solderidx].ry / 2.0 + nc * pk->r_solder1_y);
      psum += (x[PCB_C_S] - x[SOLDER_S])/pk->r_pcb_c_per_y;
      dv[nl*nr*nc + SOLDER_S] = psum / pk->c_solder_per_y;

      /* solder ball west/east	*/
      /* partition r_solder1_x among all the nr grid cells. edge cell has half the rx	*/
      psum = 0.0;
      for(i=0; i < nr; i++)
        psum += (A3D(v,solderidx,i,0,nl,nr,nc) - x[SOLDER_W]);
      psum /= (l[solderidx].rx / 2.0 + nr * pk->r_solder1_x);
      psum += (x[PCB_C_W] - x[SOLDER_W])/pk->r_pcb_c_per_x;
      dv[nl*nr*nc + SOLDER_W] = psum / pk->c_solder_per_x;

      psum = 0.0;
      for(i=0; i < nr; i++)
        psum += (A3D(v,solderidx,i,nc-1,nl,nr,nc) - x[SOLDER_E]);
      psum /= (l[solderidx].rx / 2.0 + nr * pk->r_solder1_x);
      psum += (x[PCB_C_E] - x[SOLDER_E])/pk->r_pcb_c_per_x;
      dv[nl*nr*nc + SOLDER_E] = psum / pk->c_solder_per_x;

      /* package substrate north/south	*/
      /* partition r_sub1_y among all the nc grid cells. edge cell has half the ry	*/
      psum = 0.0;
      for(j=0; j < nc; j++)
        psum += (A3D(v,subidx,0,j,nl,nr,nc) - x[SUB_N]);
      psum /= (l[subidx].ry / 2.0 + nc * pk->r_sub1_y);
      psum += (x[SOLDER_N] - x[SUB_N])/pk->r_solder_per_y;
      dv[nl*nr*nc + SUB_N] = psum / pk->c_sub_per_y;

      psum = 0.0;
      for(j=0; j < nc; j++)
        psum += (A3D(v,subidx,nr-1,j,nl,nr,nc) - x[SUB_S]);
      psum /= (l[subidx].ry / 2.0 + nc * pk->r_sub1_y);
      psum += (x[SOLDER_S] - x[SUB_S])/pk->r_solder_per_y;
      dv[nl*nr*nc + SUB_S] = psum / pk->c_sub_per_y;

      /* sub ball west/east	*/
      /* partition r_sub1_x among all the nr grid cells. edge cell has half the rx	*/
      psum = 0.0;
      for(i=0; i < nr; i++)
        psum += (A3D(v,subidx,i,0,nl,nr,nc) - x[SUB_W]);
      psum /= (l[subidx].rx / 2.0 + nr * pk->r_sub1_x);
      psum += (x[SOLDER_W] - x[SUB_W])/pk->r_solder_per_x;
      dv[nl*nr*nc + SUB_W] = psum / pk->c_sub_per_x;

      psum = 0.0;
      for(i=0; i < nr; i++)
        psum += (A3D(v,subidx,i,nc-1,nl,nr,nc) - x[SUB_E]);
      psum /= (l[subidx].rx / 2.0 + nr * pk->r_sub1_x);
      psum += (x[SOLDER_E] - x[SUB_E])/pk->r_solder_per_x;
      dv[nl*nr*nc + SUB_E] = psum / pk->c_sub_per_x;
  }
}

/* macros for calculating currents(power values)	*/
/* current(power) from the next cell north. zero if on northern boundary	*/
# define NP(l,v,n,i,j,nl,nr,nc)		((i > 0) ? ((A3D(v,n,i-1,j,nl,nr,nc)-A3D(v,n,i,j,nl,nr,nc))/l[n].ry) : 0.0)
/* current(power) from the next cell south. zero if on southern boundary	*/
# define SP(l,v,n,i,j,nl,nr,nc)		((i < nr-1) ? ((A3D(v,n,i+1,j,nl,nr,nc)-A3D(v,n,i,j,nl,nr,nc))/l[n].ry) : 0.0)
/* current(power) from the next cell east. zero if on eastern boundary	*/
# define EP(l,v,n,i,j,nl,nr,nc)		((j < nc-1) ? ((A3D(v,n,i,j+1,nl,nr,nc)-A3D(v,n,i,j,nl,nr,nc))/l[n].rx) : 0.0)
/* current(power) from the next cell west. zero if on western boundary	*/
# define WP(l,v,n,i,j,nl,nr,nc)		((j > 0) ? ((A3D(v,n,i,j-1,nl,nr,nc)-A3D(v,n,i,j,nl,nr,nc))/l[n].rx) : 0.0)
/* current(power) from the next cell below. zero if on bottom face		*/
# define BP(l,v,n,i,j,nl,nr,nc)		((n < nl-1) ? ((A3D(v,n+1,i,j,nl,nr,nc)-A3D(v,n,i,j,nl,nr,nc))/l[n].rz) : 0.0)
/* current(power) from the next cell above. zero if on top face			*/
# define AP(l,v,n,i,j,nl,nr,nc)		((n > 0) ? ((A3D(v,n-1,i,j,nl,nr,nc)-A3D(v,n,i,j,nl,nr,nc))/l[n-1].rz) : 0.0)

//BU_3D: These are the same macros as above except that they have to check
//for a resistivity value first since the lc model is not uniform across the layer
/* current(power) from the next cell north. zero if on northern boundary	*/
# define NP_det3D(l,v,n,i,j,nl,nr,nc)		((i > 0) ? ((A3D(v,n,i-1,j,nl,nr,nc)-A3D(v,n,i,j,nl,nr,nc))/find_res_3D(n, i-1, j, model,2)) : 0.0)
/* current(power) from the next cell south. zero if on southern boundary	*/
# define SP_det3D(l,v,n,i,j,nl,nr,nc)		((i < nr-1) ? ((A3D(v,n,i+1,j,nl,nr,nc)-A3D(v,n,i,j,nl,nr,nc))/find_res_3D(n, i+1, j, model,2)) : 0.0)
/* current(power) from the next cell east. zero if on eastern boundary	*/
# define EP_det3D(l,v,n,i,j,nl,nr,nc)		((j < nc-1) ? ((A3D(v,n,i,j+1,nl,nr,nc)-A3D(v,n,i,j,nl,nr,nc))/find_res_3D(n, i, j+1, model,1)) : 0.0)
/* current(power) from the next cell west. zero if on western boundary	*/
# define WP_det3D(l,v,n,i,j,nl,nr,nc)		((j > 0) ? ((A3D(v,n,i,j-1,nl,nr,nc)-A3D(v,n,i,j,nl,nr,nc))/find_res_3D(n, i, j-1, model,1)) : 0.0)
/* current(power) from the next cell below. zero if on bottom face		*/
# define BP_det3D(l,v,n,i,j,nl,nr,nc)		((n < nl-1) ? ((A3D(v,n+1,i,j,nl,nr,nc)-A3D(v,n,i,j,nl,nr,nc))/find_res_3D(n, i, j, model,3)) : 0.0)
/* current(power) from the next cell above. zero if on top face			*/
# define AP_det3D(l,v,n,i,j,nl,nr,nc)		((n > 0) ? ((A3D(v,n-1,i,j,nl,nr,nc)-A3D(v,n,i,j,nl,nr,nc))/find_res_3D(n-1, i, j, model,3)) : 0.0)
//end->BU_3D


/* compute the slope vector for the grid cells. the transient
 * equation is CdV + sum{(T - Ti)/Ri} = P 
 * so, slope = dV = [P + sum{(Ti-T)/Ri}]/C
 */
void slope_fn_grid(grid_model_t *model, const double *v,
                   grid_model_vector_t *p, double *dv)
{
  int n, i, j;
  /* sum of the currents(power values)	*/
  double psum;

  /* shortcuts for cell width(cw) and cell height(ch)	*/
  double cw = model->width / model->cols;
  double ch = model->height / model->rows;

  /* shortcuts	*/
  thermal_config_t *c = &model->config;
  layer_t *l = model->layers;
  int nl = model->n_layers;
  int nr = model->rows;
  int nc = model->cols;
  int spidx, hsidx, subidx, solderidx, pcbidx;
  int model_secondary = model->config.model_secondary;

  /* pointer to the starting address of the extra nodes	*/
  const double *x = v + nl*nr*nc;

  spidx = nl - DEFAULT_PACK_LAYERS + LAYER_SP;
  hsidx = nl - DEFAULT_PACK_LAYERS + LAYER_SINK;
  if (model_secondary) {
      subidx = LAYER_SUB;
      solderidx = LAYER_SOLDER;
      pcbidx = LAYER_PCB;	
  }

  /* for each grid cell	*/
  for(n=0; n < nl; n++)
    for(i=0; i < nr; i++)
      for(j=0; j < nc; j++) {
          /* sum the currents(power values) to cells north, south, 
           * east, west, above and below
           */
          // BU_3D: uses grid specific values for all layers 
          // spreader and heat sink layers will use uniform R
          if(model->config.detailed_3D_used == 1){
              psum = NP_det3D(l,v,n,i,j,nl,nr,nc) + SP_det3D(l,v,n,i,j,nl,nr,nc) + 
                EP_det3D(l,v,n,i,j,nl,nr,nc) + WP_det3D(l,v,n,i,j,nl,nr,nc) + 
                AP_det3D(l,v,n,i,j,nl,nr,nc) + BP_det3D(l,v,n,i,j,nl,nr,nc);
          }
          else{
              psum = NP(l,v,n,i,j,nl,nr,nc) + SP(l,v,n,i,j,nl,nr,nc) + 
                EP(l,v,n,i,j,nl,nr,nc) + WP(l,v,n,i,j,nl,nr,nc) + 
                AP(l,v,n,i,j,nl,nr,nc) + BP(l,v,n,i,j,nl,nr,nc);
          }//end->BU_3D

          /* spreader core is connected to its periphery	*/
          if (n == spidx) {
              /* northern boundary - edge cell has half the ry	*/
              if (i == 0)
                psum += (x[SP_N] - A3D(v,n,i,j,nl,nr,nc))/(l[n].ry/2.0 + nc*model->pack.r_sp1_y); 
              /* southern boundary - edge cell has half the ry	*/
              if (i == nr-1)
                psum += (x[SP_S] - A3D(v,n,i,j,nl,nr,nc))/(l[n].ry/2.0 + nc*model->pack.r_sp1_y); 
              /* eastern boundary	 - edge cell has half the rx	*/
              if (j == nc-1)
                psum += (x[SP_E] - A3D(v,n,i,j,nl,nr,nc))/(l[n].rx/2.0 + nr*model->pack.r_sp1_x); 
              /* western boundary	 - edge cell has half the rx	*/
              if (j == 0)
                psum += (x[SP_W] - A3D(v,n,i,j,nl,nr,nc))/(l[n].rx/2.0 + nr*model->pack.r_sp1_x); 
              /* heatsink core is connected to its inner periphery and ambient	*/
          } else if (n == hsidx) {
              /* all nodes are connected to the ambient	*/
              psum += (c->ambient - A3D(v,n,i,j,nl,nr,nc))/l[n].rz;
              /* northern boundary - edge cell has half the ry	*/
              if (i == 0)
                psum += (x[SINK_C_N] - A3D(v,n,i,j,nl,nr,nc))/(l[n].ry/2.0 + nc*model->pack.r_hs1_y); 
              /* southern boundary - edge cell has half the ry	*/
              if (i == nr-1)
                psum += (x[SINK_C_S] - A3D(v,n,i,j,nl,nr,nc))/(l[n].ry/2.0 + nc*model->pack.r_hs1_y); 
              /* eastern boundary	 - edge cell has half the rx	*/
              if (j == nc-1)
                psum += (x[SINK_C_E] - A3D(v,n,i,j,nl,nr,nc))/(l[n].rx/2.0 + nr*model->pack.r_hs1_x); 
              /* western boundary	 - edge cell has half the rx	*/
              if (j == 0)
                psum += (x[SINK_C_W] - A3D(v,n,i,j,nl,nr,nc))/(l[n].rx/2.0 + nr*model->pack.r_hs1_x); 
          }	else if (model_secondary && n == pcbidx) {
              /* all nodes are connected to the ambient	*/
              psum += (c->ambient - A3D(v,n,i,j,nl,nr,nc))/(model->config.r_convec_sec * 
                                                            (model->config.s_pcb * model->config.s_pcb) / (cw * ch));
              /* northern boundary - edge cell has half the ry	*/
              if (i == 0)
                psum += (x[PCB_C_N] - A3D(v,n,i,j,nl,nr,nc))/(l[n].ry/2.0 + nc*model->pack.r_pcb1_y); 
              /* southern boundary - edge cell has half the ry	*/
              if (i == nr-1)
                psum += (x[PCB_C_S] - A3D(v,n,i,j,nl,nr,nc))/(l[n].ry/2.0 + nc*model->pack.r_pcb1_y); 
              /* eastern boundary	 - edge cell has half the rx	*/
              if (j == nc-1)
                psum += (x[PCB_C_E] - A3D(v,n,i,j,nl,nr,nc))/(l[n].rx/2.0 + nr*model->pack.r_pcb1_x); 
              /* western boundary	 - edge cell has half the rx	*/
              if (j == 0)
                psum += (x[PCB_C_W] - A3D(v,n,i,j,nl,nr,nc))/(l[n].rx/2.0 + nr*model->pack.r_pcb1_x); 
          }	else if (model_secondary && n == subidx) {
              /* northern boundary - edge cell has half the ry	*/
              if (i == 0)
                psum += (x[SUB_N] - A3D(v,n,i,j,nl,nr,nc))/(l[n].ry/2.0 + nc*model->pack.r_sub1_y); 
              /* southern boundary - edge cell has half the ry	*/
              if (i == nr-1)
                psum += (x[SUB_S] - A3D(v,n,i,j,nl,nr,nc))/(l[n].ry/2.0 + nc*model->pack.r_sub1_y); 
              /* eastern boundary	 - edge cell has half the rx	*/
              if (j == nc-1)
                psum += (x[SUB_E] - A3D(v,n,i,j,nl,nr,nc))/(l[n].rx/2.0 + nr*model->pack.r_sub1_x); 
              /* western boundary	 - edge cell has half the rx	*/
              if (j == 0)
                psum += (x[SUB_W] - A3D(v,n,i,j,nl,nr,nc))/(l[n].rx/2.0 + nr*model->pack.r_sub1_x); 
          }	else if (model_secondary && n == solderidx) {
              /* northern boundary - edge cell has half the ry	*/
              if (i == 0)
                psum += (x[SOLDER_N] - A3D(v,n,i,j,nl,nr,nc))/(l[n].ry/2.0 + nc*model->pack.r_solder1_y); 
              /* southern boundary - edge cell has half the ry	*/
              if (i == nr-1)
                psum += (x[SOLDER_S] - A3D(v,n,i,j,nl,nr,nc))/(l[n].ry/2.0 + nc*model->pack.r_solder1_y); 
              /* eastern boundary	 - edge cell has half the rx	*/
              if (j == nc-1)
                psum += (x[SOLDER_E] - A3D(v,n,i,j,nl,nr,nc))/(l[n].rx/2.0 + nr*model->pack.r_solder1_x); 
              /* western boundary	 - edge cell has half the rx	*/
              if (j == 0)
                psum += (x[SOLDER_W] - A3D(v,n,i,j,nl,nr,nc))/(l[n].rx/2.0 + nr*model->pack.r_solder1_x); 
          }

          /* update the current cell's temperature	*/	   
          if(model->config.detailed_3D_used == 1)//BU_3D: use find_cap_3D is detailed_3D model is used.
            A3D(dv,n,i,j,nl,nr,nc) = (p->cuboid[n][i][j] + psum) / find_cap_3D(n, i, j, model);
          else
            A3D(dv,n,i,j,nl,nr,nc) = (p->cuboid[n][i][j] + psum) / l[n].c;
      }
  /* for each grid cell	*/
  slope_fn_pack(model, v, p, dv);
}

static int steady_node_count(grid_model_t *model)
{
  int package_nodes = model->config.model_secondary ? EXTRA + EXTRA_SEC : EXTRA;
  size_t plane_nodes;
  size_t maximum_grid_nodes = (size_t)INT_MAX - package_nodes;

  if (model->n_layers <= 0 || model->rows <= 0 || model->cols <= 0)
    fatal("invalid grid dimensions in steady-state solver\n");
  plane_nodes = (size_t)model->rows * (size_t)model->cols;
  if ((size_t)model->n_layers > maximum_grid_nodes / plane_nodes)
    fatal("steady-state node count exceeds supported integer range\n");
  return (int)((size_t)model->n_layers * plane_nodes + package_nodes);
}

#if GATE_GS1_ARM > 0

#define GS1_PATH_SIZE 4096

static const char * const gs1_package_node_names[EXTRA + EXTRA_SEC] = {
  "SP_W", "SP_E", "SP_N", "SP_S",
  "SINK_C_W", "SINK_C_E", "SINK_C_N", "SINK_C_S",
  "SINK_W", "SINK_E", "SINK_N", "SINK_S",
  "SUB_W", "SUB_E", "SUB_N", "SUB_S",
  "SOLDER_W", "SOLDER_E", "SOLDER_N", "SOLDER_S",
  "PCB_C_W", "PCB_C_E", "PCB_C_N", "PCB_C_S",
  "PCB_W", "PCB_E", "PCB_N", "PCB_S"
};

static const char *gs1_arm_name(void)
{
#if GATE_GS1_ARM == GS1_ARM_GS_NATIVE
  return "gs_native";
#elif GATE_GS1_ARM == GS1_ARM_SUPERLU
  return "superlu";
#elif GATE_GS1_ARM == GS1_ARM_MFPCG
  return "mfpcg";
#elif GATE_GS1_ARM == GS1_ARM_CSR_PCG
  return "csr_pcg";
#elif GATE_GS1_ARM == GS1_ARM_STRUCTURED_PCG
  return "structured_pcg";
#else
#error invalid GATE_GS1_ARM
#endif
}

static const char *gs1_output_prefix(void)
{
  const char *prefix = getenv("HOTSPOT_GS1_PREFIX");
  if (!prefix || !prefix[0])
    fatal("HOTSPOT_GS1_PREFIX is required for a G-S1 benchmark build\n");
  return prefix;
}

static FILE *gs1_open_output(const char *prefix, const char *suffix,
                             const char *mode)
{
  char path[GS1_PATH_SIZE];
  char message[GS1_PATH_SIZE + 80];
  int length = snprintf(path, sizeof(path), "%s%s", prefix, suffix);
  FILE *stream;
  if (length < 0 || length >= (int)sizeof(path))
    fatal("HOTSPOT_GS1_PREFIX is too long\n");
  stream = fopen(path, mode);
  if (!stream) {
    snprintf(message, sizeof(message), "unable to open G-S1 output %s\n", path);
    fatal(message);
  }
  return stream;
}

static void gs1_write_bytes(const char *prefix, const char *suffix,
                            const void *values, size_t size, size_t count)
{
  FILE *stream = gs1_open_output(prefix, suffix, "wb");
  if (fwrite(values, size, count, stream) != count) {
    fclose(stream);
    fatal("unable to write complete G-S1 binary output\n");
  }
  if (fclose(stream) != 0)
    fatal("unable to close G-S1 binary output\n");
}

static void gs1_write_completion(const char *prefix, const char *suffix)
{
  static const char marker[] = "complete\n";
  gs1_write_bytes(prefix, suffix, marker, 1, sizeof(marker) - 1);
}

static void gs1_write_node_map(grid_model_t *model, const char *prefix)
{
  int failed;
  int index = 0;
  int layer, row, column, package_index;
  int package_nodes = model->config.model_secondary ? EXTRA + EXTRA_SEC : EXTRA;
  FILE *stream = gs1_open_output(prefix, ".nodes.tsv", "w");
  for (layer = 0; layer < model->n_layers; layer++)
    for (row = 0; row < model->rows; row++)
      for (column = 0; column < model->cols; column++)
        fprintf(stream, "%d\tgrid\t%d\t%d\t%d\n",
                index++, layer, row, column);
  for (package_index = 0; package_index < package_nodes; package_index++)
    fprintf(stream, "%d\tpackage\t%s\n", index++,
            gs1_package_node_names[package_index]);
  failed = ferror(stream);
  if (fclose(stream) != 0)
    failed = 1;
  if (failed)
    fatal("unable to write complete G-S1 node map\n");
}

static void gs1_write_state(grid_model_t *model, const double *block_state,
                            const char *prefix)
{
  const double *state = model->last_steady->cuboid[0][0];
  int package_nodes = model->config.model_secondary ? EXTRA + EXTRA_SEC : EXTRA;
  gs1_write_bytes(prefix, ".state.bin", state, sizeof(double),
                  (size_t)steady_node_count(model));
  gs1_write_bytes(prefix, ".block-state.bin", block_state, sizeof(double),
                  (size_t)(model->total_n_blocks + package_nodes));
  gs1_write_node_map(model, prefix);
  gs1_write_completion(prefix, ".state.complete");
}

static void gs1_write_metadata(grid_model_t *model, const char *prefix)
{
  unsigned short endian_probe = 1;
  int failed;
  int count = steady_node_count(model);
  int package_nodes = model->config.model_secondary ? EXTRA + EXTRA_SEC : EXTRA;
  double managed = model->gs1_coefficient_s + model->gs1_initialization_s +
    model->gs1_assembly_s + model->gs1_factor_s + model->gs1_solve_s;
  const double *state = model->last_steady->cuboid[0][0];
  double peak = state[0];
  int peak_index = 0;
  int index;
  FILE *stream;
  for (index = 1; index < count; index++)
    if (state[index] > peak) {
      peak = state[index];
      peak_index = index;
    }
  stream = gs1_open_output(prefix, ".meta", "w");
  fprintf(stream, "format=hotspot-gs1-v1\n");
  fprintf(stream, "arm=%s\n", gs1_arm_name());
  fprintf(stream, "dispatch=%s\n", gs1_arm_name());
  fprintf(stream, "terminal=%s\n",
          GATE_GS1_ARM == GS1_ARM_GS_NATIVE ? "native_complete" :
          ((GATE_GS1_ARM == GS1_ARM_MFPCG ||
            GATE_GS1_ARM == GS1_ARM_CSR_PCG ||
            GATE_GS1_ARM == GS1_ARM_STRUCTURED_PCG) &&
           (!(model->gs1_relative_residual <= PCG_RELATIVE_RESIDUAL) ||
            !isfinite(model->gs1_relative_residual))) ?
          "iteration_limit" : "converged");
  fprintf(stream, "stopping_contract=%s\n",
          GATE_GS1_ARM == GS1_ARM_GS_NATIVE ? "native_stop" : "candidate_stop");
  fprintf(stream, "scalar=float64-native\n");
  fprintf(stream, "byte_order=%s\n",
          *((unsigned char *)&endian_probe) ? "little" : "big");
  fprintf(stream, "temperature_unit=K\n");
  fprintf(stream, "state_order=layer-row-column-then-package\n");
  fprintf(stream, "node_count=%d\n", count);
  fprintf(stream, "rows=%d\n", model->rows);
  fprintf(stream, "cols=%d\n", model->cols);
  fprintf(stream, "layers=%d\n", model->n_layers);
  fprintf(stream, "package_nodes=%d\n", package_nodes);
  fprintf(stream, "block_state_count=%d\n",
          model->total_n_blocks + package_nodes);
  fprintf(stream, "model_secondary=%d\n", model->config.model_secondary);
  fprintf(stream, "detailed_3D=%d\n", model->config.detailed_3D_used);
  fprintf(stream, "time_input_model_s=%.17g\n", model->gs1_input_model_s);
  fprintf(stream, "time_coefficient_cache_s=%.17g\n", model->gs1_coefficient_s);
  fprintf(stream, "time_initialization_s=%.17g\n", model->gs1_initialization_s);
  fprintf(stream, "time_assembly_s=%.17g\n", model->gs1_assembly_s);
  fprintf(stream, "time_preconditioner_factor_s=%.17g\n", model->gs1_factor_s);
  fprintf(stream, "time_solve_s=%.17g\n", model->gs1_solve_s);
  fprintf(stream, "time_output_s=%.17g\n", model->gs1_output_s);
  fprintf(stream, "time_managed_s=%.17g\n", managed);
  fprintf(stream, "iterations=%lld\n", model->gs1_iterations);
  fprintf(stream, "solution_iterations=%lld\n",
          model->gs1_solution_iterations);
  fprintf(stream, "replay_iterations=%lld\n",
          model->gs1_replay_iterations);
  fprintf(stream, "operator_applications=%lld\n",
          model->gs1_operator_applications);
  fprintf(stream, "residual_evaluations=%lld\n",
          model->gs1_residual_evaluations);
  fprintf(stream, "residual_replacements=%lld\n",
          model->gs1_residual_replacements);
  fprintf(stream, "solver_vector_bytes=%lu\n",
          (unsigned long)model->gs1_solver_vector_bytes);
  fprintf(stream, "operator_bytes=%lu\n",
          (unsigned long)model->gs1_operator_bytes);
  fprintf(stream, "direct_permutation_bytes=%lu\n",
          (unsigned long)model->gs1_direct_permutation_bytes);
  fprintf(stream, "direct_l_nnz=%lld\n", model->gs1_direct_l_nnz);
  fprintf(stream, "direct_u_nnz=%lld\n", model->gs1_direct_u_nnz);
  if (isfinite(model->gs1_relative_residual))
    fprintf(stream, "internal_relative_residual=%.17g\n",
            model->gs1_relative_residual);
  else
    fprintf(stream, "internal_relative_residual=not_available\n");
  if (isfinite(model->gs1_native_delta))
    fprintf(stream, "native_delta=%.17g\n", model->gs1_native_delta);
  else
    fprintf(stream, "native_delta=not_available\n");
  fprintf(stream, "action_audit_probes=%d\n",
          model->gs1_action_audit_probes);
  if (isfinite(model->gs1_action_audit_max_abs_error))
    fprintf(stream, "action_audit_max_abs_error=%.17g\n",
            model->gs1_action_audit_max_abs_error);
  else
    fprintf(stream, "action_audit_max_abs_error=not_run\n");
  if (isfinite(model->gs1_action_audit_max_relative_error))
    fprintf(stream, "action_audit_max_relative_error=%.17g\n",
            model->gs1_action_audit_max_relative_error);
  else
    fprintf(stream, "action_audit_max_relative_error=not_run\n");
  fprintf(stream, "peak_temperature_K=%.17g\n", peak);
  fprintf(stream, "peak_state_index=%d\n", peak_index);
  failed = ferror(stream);
  if (fclose(stream) != 0)
    failed = 1;
  if (failed)
    fatal("unable to write complete G-S1 metadata\n");
  gs1_write_completion(prefix, ".complete");
}

#endif

#if GATE_GS0 > 0

#define GS0_PATH_SIZE 4096
#define GS0_MAX_NODE_COUNT 2048

static const char * const gs0_package_node_names[EXTRA + EXTRA_SEC] = {
  "SP_W", "SP_E", "SP_N", "SP_S",
  "SINK_C_W", "SINK_C_E", "SINK_C_N", "SINK_C_S",
  "SINK_W", "SINK_E", "SINK_N", "SINK_S",
  "SUB_W", "SUB_E", "SUB_N", "SUB_S",
  "SOLDER_W", "SOLDER_E", "SOLDER_N", "SOLDER_S",
  "PCB_C_W", "PCB_C_E", "PCB_C_N", "PCB_C_S",
  "PCB_W", "PCB_E", "PCB_N", "PCB_S"
};

static void gs0_validate_package_node_order(void)
{
  static const int package_node_indices[EXTRA + EXTRA_SEC] = {
    SP_W, SP_E, SP_N, SP_S,
    SINK_C_W, SINK_C_E, SINK_C_N, SINK_C_S,
    SINK_W, SINK_E, SINK_N, SINK_S,
    SUB_W, SUB_E, SUB_N, SUB_S,
    SOLDER_W, SOLDER_E, SOLDER_N, SOLDER_S,
    PCB_C_W, PCB_C_E, PCB_C_N, PCB_C_S,
    PCB_W, PCB_E, PCB_N, PCB_S
  };
  int index;
  for (index = 0; index < EXTRA + EXTRA_SEC; index++)
    if (package_node_indices[index] != index)
      fatal("G-S0 package node names do not match HotSpot indices\n");
}

static const char *gs0_output_prefix(void)
{
  const char *prefix = getenv("HOTSPOT_GS0_PREFIX");
  return prefix && prefix[0] ? prefix : NULL;
}

static FILE *gs0_open_output(const char *prefix, const char *suffix,
                             const char *mode)
{
  char path[GS0_PATH_SIZE];
  char message[GS0_PATH_SIZE + 80];
  int length = snprintf(path, sizeof(path), "%s%s", prefix, suffix);
  FILE *stream;
  if (length < 0 || length >= (int)sizeof(path))
    fatal("HOTSPOT_GS0_PREFIX is too long\n");
  stream = fopen(path, mode);
  if (!stream) {
    snprintf(message, sizeof(message), "unable to open G-S0 output %s\n", path);
    fatal(message);
  }
  return stream;
}

static void gs0_write_doubles(const char *prefix, const char *suffix,
                              const double *values, size_t count)
{
  FILE *stream = gs0_open_output(prefix, suffix, "wb");
  if (fwrite(values, sizeof(double), count, stream) != count) {
    fclose(stream);
    fatal("unable to write complete G-S0 binary output\n");
  }
  if (fclose(stream) != 0)
    fatal("unable to close G-S0 binary output\n");
}

static void gs0_write_completion(const char *prefix, const char *suffix)
{
  static const char marker[] = "complete\n";
  FILE *stream = gs0_open_output(prefix, suffix, "wb");
  if (fwrite(marker, 1, sizeof(marker) - 1, stream) != sizeof(marker) - 1) {
    fclose(stream);
    fatal("unable to write G-S0 completion marker\n");
  }
  if (fclose(stream) != 0)
    fatal("unable to close G-S0 completion marker\n");
}

static void gs0_write_node_map(grid_model_t *model, const char *prefix)
{
  int failed;
  int index = 0;
  int layer, row, column, package_index;
  int package_nodes = model->config.model_secondary ? EXTRA + EXTRA_SEC : EXTRA;
  FILE *stream = gs0_open_output(prefix, ".nodes.tsv", "w");
  for (layer = 0; layer < model->n_layers; layer++)
    for (row = 0; row < model->rows; row++)
      for (column = 0; column < model->cols; column++)
        fprintf(stream, "%d\tgrid\t%d\t%d\t%d\n",
                index++, layer, row, column);
  for (package_index = 0; package_index < package_nodes; package_index++)
    fprintf(stream, "%d\tpackage\t%s\n", index++,
            gs0_package_node_names[package_index]);
  failed = ferror(stream);
  if (fclose(stream) != 0)
    failed = 1;
  if (failed)
    fatal("unable to write complete G-S0 node map\n");
}

static void gs0_write_metadata(grid_model_t *model, const char *prefix,
                               const char *arm, int probe_count,
                               int basis_action_count, int input_check_count,
                               int input_mutation_count,
                               double coupling_asymmetry)
{
  unsigned short endian_probe = 1;
  int failed;
  int package_nodes = model->config.model_secondary ? EXTRA + EXTRA_SEC : EXTRA;
  FILE *stream = gs0_open_output(prefix, ".meta", "w");
  fprintf(stream, "format=hotspot-gs0-v1\n");
  fprintf(stream, "arm=%s\n", arm);
  fprintf(stream, "dispatch=%s\n", SUPERLU > 0 ? "superlu" : "mfpcg");
  fprintf(stream, "scalar=float64-native\n");
  fprintf(stream, "byte_order=%s\n",
          *((unsigned char *)&endian_probe) ? "little" : "big");
  fprintf(stream, "matrix_order=row-major\n");
  fprintf(stream, "node_count=%d\n", steady_node_count(model));
  fprintf(stream, "rows=%d\n", model->rows);
  fprintf(stream, "cols=%d\n", model->cols);
  fprintf(stream, "layers=%d\n", model->n_layers);
  fprintf(stream, "package_nodes=%d\n", package_nodes);
  fprintf(stream, "probe_count=%d\n", probe_count);
  fprintf(stream, "basis_action_count=%d\n", basis_action_count);
  fprintf(stream, "input_check_count=%d\n", input_check_count);
  fprintf(stream, "input_mutation_count=%d\n", input_mutation_count);
  fprintf(stream, "input_unchanged=%d\n", input_mutation_count == 0);
  fprintf(stream, "action_entry=%s\n",
          probe_count > 0 ? "steady_matvec" : "superlu_builders");
  fprintf(stream, "pcg_started_at_operator_export=0\n");
  fprintf(stream, "package_node_order_checked=1\n");
#if SUPERLU > 0
  fprintf(stream, "superlu_index_bytes=%lu\n", (unsigned long)sizeof(int_t));
#endif
  if (isfinite(coupling_asymmetry))
    fprintf(stream, "grid_coupling_asymmetry=%.17g\n", coupling_asymmetry);
  else
    fprintf(stream, "grid_coupling_asymmetry=not_applicable\n");
  if (probe_count > 0) {
    fprintf(stream, "linearity_0=4,0,0.25,1,-1.5\n");
    fprintf(stream, "linearity_1=5,2,-2,3,0.75\n");
  }
  failed = ferror(stream);
  if (fclose(stream) != 0)
    failed = 1;
  if (failed)
    fatal("unable to write complete G-S0 metadata\n");
  gs0_validate_package_node_order();
  gs0_write_node_map(model, prefix);
}

static void gs0_write_state(grid_model_t *model, const double *state)
{
  const char *prefix = gs0_output_prefix();
  if (prefix) {
    gs0_write_doubles(prefix, ".state.bin", state,
                      (size_t)steady_node_count(model));
    gs0_write_completion(prefix, ".state.complete");
  }
}

#endif

#if SUPERLU < 1 || GATE_GS1_ARM == GS1_ARM_CSR_PCG || \
    GATE_GS1_ARM == GS1_ARM_STRUCTURED_PCG

#define PCG_MAX_ITERATIONS 5000
#define PCG_SYMMETRY_TOLERANCE 1.0e-10

static double dot_product(const double *a, const double *b, int count)
{
  double sum = 0.0;
  int i;
  for (i = 0; i < count; i++)
    sum += a[i] * b[i];
  return sum;
}

/* slope_fn_grid returns C^-1(b-Ax).  Multiplying by the exact node
 * capacitance recovers the steady-state residual without duplicating the
 * detailed-grid or package conductance equations. */
static void build_capacitance_vector(grid_model_t *model, double *cap)
{
  int n, i, j;
  int nl = model->n_layers;
  int nr = model->rows;
  int nc = model->cols;
  int offset = nl * nr * nc;
  package_RC_t *pk = &model->pack;

  for (n = 0; n < nl; n++)
    for (i = 0; i < nr; i++)
      for (j = 0; j < nc; j++)
        cap[n * nr * nc + i * nc + j] =
          find_cap_3D(n, i, j, model);

  cap[offset + SP_W] = cap[offset + SP_E] = pk->c_sp_per_x;
  cap[offset + SP_N] = cap[offset + SP_S] = pk->c_sp_per_y;
  cap[offset + SINK_C_W] = cap[offset + SINK_C_E] =
    pk->c_hs_c_per_x + pk->c_amb_c_per_x;
  cap[offset + SINK_C_N] = cap[offset + SINK_C_S] =
    pk->c_hs_c_per_y + pk->c_amb_c_per_y;
  cap[offset + SINK_W] = cap[offset + SINK_E] =
    pk->c_hs_per + pk->c_amb_per;
  cap[offset + SINK_N] = cap[offset + SINK_S] =
    pk->c_hs_per + pk->c_amb_per;

  if (model->config.model_secondary) {
    cap[offset + SUB_W] = cap[offset + SUB_E] = pk->c_sub_per_x;
    cap[offset + SUB_N] = cap[offset + SUB_S] = pk->c_sub_per_y;
    cap[offset + SOLDER_W] = cap[offset + SOLDER_E] = pk->c_solder_per_x;
    cap[offset + SOLDER_N] = cap[offset + SOLDER_S] = pk->c_solder_per_y;
    cap[offset + PCB_C_W] = cap[offset + PCB_C_E] =
      pk->c_pcb_c_per_x + pk->c_amb_sec_c_per_x;
    cap[offset + PCB_C_N] = cap[offset + PCB_C_S] =
      pk->c_pcb_c_per_y + pk->c_amb_sec_c_per_y;
    cap[offset + PCB_W] = cap[offset + PCB_E] =
      pk->c_pcb_per + pk->c_amb_sec_per;
    cap[offset + PCB_N] = cap[offset + PCB_S] =
      pk->c_pcb_per + pk->c_amb_sec_per;
  }
}

static void steady_residual(grid_model_t *model, grid_model_vector_t *power,
                            const double *cap, const double *x, double *residual)
{
  int count = steady_node_count(model);
  int i;
#if GATE_GS1_ARM == GS1_ARM_MFPCG || \
    GATE_GS1_ARM == GS1_ARM_CSR_PCG || \
    GATE_GS1_ARM == GS1_ARM_STRUCTURED_PCG
  model->gs1_residual_evaluations++;
#endif
  slope_fn_grid(model, x, power, residual);
  for (i = 0; i < count; i++)
    residual[i] *= cap[i];
}

static void steady_matvec(grid_model_t *model, grid_model_vector_t *power,
                          const double *cap, const double *rhs,
                          const double *x, double *product, double *scaled)
{
  int count = steady_node_count(model);
  int i;
  double magnitude = 0.0;
  double scale;
#if GATE_GS1_ARM == GS1_ARM_MFPCG
  model->gs1_operator_applications++;
#endif
  for (i = 0; i < count; i++)
    magnitude = MAX(magnitude, fabs(x[i]));
  if (magnitude == 0.0) {
    zero_dvector(product, count);
    return;
  }
  scale = 1.0 / magnitude;
  for (i = 0; i < count; i++)
    scaled[i] = x[i] * scale;
  steady_residual(model, power, cap, scaled, product);
  for (i = 0; i < count; i++)
    product[i] = (rhs[i] - product[i]) / scale;
}

#if SUPERLU > 0
/* Generic sparse storage for the G-S1.5 CSR arm.  This is deliberately not
 * SuperLU's CSC: the candidate is a row-oriented SpMV PCG solve after the
 * exact HotSpot matrix builder has completed. */
typedef struct {
  int count;
  int_t nnz;
  int_t *row_offsets;
  int_t *columns;
  double *values;
} steady_csr_t;

static void free_steady_csr(steady_csr_t *csr)
{
  free(csr->row_offsets);
  free(csr->columns);
  free(csr->values);
  memset(csr, 0, sizeof(*csr));
}

static void csr_from_superlu(const SuperMatrix *matrix, steady_csr_t *csr,
                             double *diagonal)
{
  NCformat *source = (NCformat *)matrix->Store;
  int_t *cursor;
  int column;
  int_t index;
  if (!source || matrix->nrow != matrix->ncol || matrix->nrow <= 0)
    fatal("invalid explicit matrix for CSR-PCG\n");
  csr->count = matrix->nrow;
  csr->nnz = source->nnz;
  csr->row_offsets = (int_t *)calloc((size_t)csr->count + 1, sizeof(int_t));
  csr->columns = (int_t *)calloc((size_t)csr->nnz, sizeof(int_t));
  csr->values = (double *)calloc((size_t)csr->nnz, sizeof(double));
  if (!csr->row_offsets || !csr->columns || !csr->values)
    fatal("memory allocation failed building CSR-PCG representation\n");
  for (index = 0; index < csr->nnz; index++) {
    int_t row = source->rowind[index];
    if (row < 0 || row >= csr->count)
      fatal("invalid row index in explicit HotSpot matrix\n");
    csr->row_offsets[row + 1]++;
  }
  for (column = 0; column < csr->count; column++)
    csr->row_offsets[column + 1] += csr->row_offsets[column];
  cursor = (int_t *)calloc((size_t)csr->count, sizeof(int_t));
  if (!cursor)
    fatal("memory allocation failed building CSR row cursors\n");
  memcpy(cursor, csr->row_offsets, (size_t)csr->count * sizeof(int_t));
  for (column = 0; column < csr->count; column++)
    for (index = source->colptr[column]; index < source->colptr[column + 1];
         index++) {
      int_t row = source->rowind[index];
      int_t destination = cursor[row]++;
      csr->columns[destination] = column;
      csr->values[destination] = ((double *)source->nzval)[index];
    }
  free(cursor);
  for (column = 0; column < csr->count; column++) {
    int_t entry;
    diagonal[column] = 0.0;
    for (entry = csr->row_offsets[column];
         entry < csr->row_offsets[column + 1]; entry++)
      if (csr->columns[entry] == column) {
        diagonal[column] = csr->values[entry];
        break;
      }
    if (!(diagonal[column] > 0.0) || !isfinite(diagonal[column]))
      fatal("invalid CSR Jacobi diagonal\n");
  }
}

static void csr_matvec(const steady_csr_t *csr, const double *x,
                       double *product)
{
  int row;
  for (row = 0; row < csr->count; row++) {
    int_t entry;
    double sum = 0.0;
    for (entry = csr->row_offsets[row]; entry < csr->row_offsets[row + 1];
         entry++)
      sum += csr->values[entry] * x[csr->columns[entry]];
    product[row] = sum;
  }
}

/* Exact hybrid structured representation.  Every regular grid-to-grid
 * coupling is one of the positive-coordinate directions; all connections
 * involving package nodes remain explicit compact edges.  Thus the format
 * preserves package semantics rather than treating HotSpot as a seven-point
 * stencil. */
typedef struct {
  int count;
  int grid_count;
  double *east;
  double *south;
  double *below;
  int_t package_edge_count;
  int_t *edge_left;
  int_t *edge_right;
  double *edge_value;
} steady_structured_t;

static void free_steady_structured(steady_structured_t *structured)
{
  free(structured->east);
  free(structured->south);
  free(structured->below);
  free(structured->edge_left);
  free(structured->edge_right);
  free(structured->edge_value);
  memset(structured, 0, sizeof(*structured));
}

static int structured_grid_direction(int column, int row, int grid_count,
                                     int rows, int cols)
{
  int plane = rows * cols;
  int column_layer = column / plane;
  int column_in_plane = column % plane;
  if (column < 0 || row < 0 || column >= grid_count || row >= grid_count)
    return -1;
  if (row == column + 1 && column_in_plane % cols + 1 < cols)
    return 0;
  if (row == column + cols && column_in_plane / cols + 1 < rows)
    return 1;
  if (row == column + plane && column_layer + 1 < grid_count / plane)
    return 2;
  return -1;
}

static void structured_from_superlu(const SuperMatrix *matrix,
                                    steady_structured_t *structured,
                                    double *diagonal, int rows, int cols,
                                    int layers)
{
  NCformat *source = (NCformat *)matrix->Store;
  int grid_count = rows * cols * layers;
  int column;
  int_t index, edge_cursor = 0;
  if (!source || matrix->nrow != matrix->ncol || matrix->nrow <= grid_count)
    fatal("invalid explicit matrix for structured-PCG\n");
  structured->count = matrix->nrow;
  structured->grid_count = grid_count;
  structured->east = (double *)calloc((size_t)grid_count, sizeof(double));
  structured->south = (double *)calloc((size_t)grid_count, sizeof(double));
  structured->below = (double *)calloc((size_t)grid_count, sizeof(double));
  if (!structured->east || !structured->south || !structured->below)
    fatal("memory allocation failed building structured grid edges\n");
  for (column = 0; column < structured->count; column++)
    for (index = source->colptr[column]; index < source->colptr[column + 1];
         index++) {
      int row = source->rowind[index];
      if (row < column)
        continue;
      if (row == column)
        continue;
      if (column >= grid_count || row >= grid_count) {
        structured->package_edge_count++;
      } else if (structured_grid_direction(column, row, grid_count, rows, cols) < 0) {
        fatal("explicit HotSpot matrix has unsupported grid coupling\n");
      }
    }
  structured->edge_left = (int_t *)calloc((size_t)structured->package_edge_count,
                                           sizeof(int_t));
  structured->edge_right = (int_t *)calloc((size_t)structured->package_edge_count,
                                            sizeof(int_t));
  structured->edge_value = (double *)calloc((size_t)structured->package_edge_count,
                                             sizeof(double));
  if (structured->package_edge_count && (!structured->edge_left ||
      !structured->edge_right || !structured->edge_value))
    fatal("memory allocation failed building structured package edges\n");
  for (column = 0; column < structured->count; column++) {
    diagonal[column] = 0.0;
    for (index = source->colptr[column]; index < source->colptr[column + 1];
         index++) {
      int row = source->rowind[index];
      double value = ((double *)source->nzval)[index];
      if (row == column) {
        diagonal[column] = value;
      } else if (row > column) {
        if (column >= grid_count || row >= grid_count) {
          structured->edge_left[edge_cursor] = column;
          structured->edge_right[edge_cursor] = row;
          structured->edge_value[edge_cursor++] = value;
        } else {
          int direction = structured_grid_direction(column, row, grid_count,
                                                    rows, cols);
          double *destination = direction == 0 ? structured->east :
            (direction == 1 ? structured->south : structured->below);
          if (direction < 0 || destination[column] != 0.0)
            fatal("duplicate or invalid structured grid edge\n");
          destination[column] = value;
        }
      }
    }
    if (!(diagonal[column] > 0.0) || !isfinite(diagonal[column]))
      fatal("invalid structured Jacobi diagonal\n");
  }
  if (edge_cursor != structured->package_edge_count)
    fatal("incomplete structured package-edge conversion\n");
}

static void structured_matvec(const steady_structured_t *structured,
                              const double *diagonal, const double *x,
                              double *product, int rows, int cols)
{
  int index;
  int plane = rows * cols;
  for (index = 0; index < structured->count; index++)
    product[index] = diagonal[index] * x[index];
  for (index = 0; index < structured->grid_count; index++) {
    if (structured->east[index] != 0.0) {
      int neighbor = index + 1;
      product[index] += structured->east[index] * x[neighbor];
      product[neighbor] += structured->east[index] * x[index];
    }
    if (structured->south[index] != 0.0) {
      int neighbor = index + cols;
      product[index] += structured->south[index] * x[neighbor];
      product[neighbor] += structured->south[index] * x[index];
    }
    if (structured->below[index] != 0.0) {
      int neighbor = index + plane;
      product[index] += structured->below[index] * x[neighbor];
      product[neighbor] += structured->below[index] * x[index];
    }
  }
  for (index = 0; index < structured->package_edge_count; index++) {
    int left = structured->edge_left[index];
    int right = structured->edge_right[index];
    double value = structured->edge_value[index];
    product[left] += value * x[right];
    product[right] += value * x[left];
  }
}
#endif

static double relative_difference(double a, double b)
{
  double scale = MAX(MAX(fabs(a), fabs(b)), 1.0);
  return fabs(a - b) / scale;
}

#if SUPERLU > 0 && (GATE_GS1_ARM == GS1_ARM_CSR_PCG || \
                     GATE_GS1_ARM == GS1_ARM_STRUCTURED_PCG)

#define GS1_ACTION_AUDIT_PROBE_COUNT 3
#define GS1_ACTION_AUDIT_TOLERANCE 1.0e-12

/* This admission-only check compares each explicit PCG representation with
 * the original detailed-3D operator on deterministic vectors.  It is opt-in
 * so measured runs do not charge the audit to solver time. */
static void gs1_audit_explicit_representation(
    grid_model_t *model, grid_model_vector_t *power, const double *cap,
    const double *rhs, const double *diagonal,
#if GATE_GS1_ARM == GS1_ARM_CSR_PCG
    const steady_csr_t *csr,
#else
    const steady_structured_t *structured,
#endif
    double *input, double *reference, double *candidate, double *scaled)
{
  const char *enabled = getenv("HOTSPOT_GS1_ACTION_AUDIT");
  int count = steady_node_count(model);
  int probe, index;
  double max_abs_error = 0.0;
  double max_relative_error = 0.0;

  if (!enabled || strcmp(enabled, "1") != 0)
    return;
  for (probe = 0; probe < GS1_ACTION_AUDIT_PROBE_COUNT; probe++) {
    for (index = 0; index < count; index++) {
      if (probe == 0)
        input[index] = (double)((index % 17) - 8);
      else if (probe == 1)
        input[index] = (double)(((index * 37) % 101) - 50) / 8.0;
      else
        input[index] = 0.25 * (double)((index % 17) - 8) -
          1.5 * (double)(((index * 37) % 101) - 50) / 8.0;
    }
    steady_matvec(model, power, cap, rhs, input, reference, scaled);
#if GATE_GS1_ARM == GS1_ARM_CSR_PCG
    csr_matvec(csr, input, candidate);
#else
    structured_matvec(structured, diagonal, input, candidate,
                      model->rows, model->cols);
#endif
    for (index = 0; index < count; index++) {
      double absolute_error = fabs(candidate[index] - reference[index]);
      max_abs_error = MAX(max_abs_error, absolute_error);
      max_relative_error = MAX(max_relative_error,
                               relative_difference(candidate[index],
                                                   reference[index]));
    }
  }
  model->gs1_action_audit_probes = GS1_ACTION_AUDIT_PROBE_COUNT;
  model->gs1_action_audit_max_abs_error = max_abs_error;
  model->gs1_action_audit_max_relative_error = max_relative_error;
  if (max_relative_error > GS1_ACTION_AUDIT_TOLERANCE)
    fatal("explicit PCG representation differs from HotSpot operator\n");
}

#endif

/* CG requires a symmetric operator.  Detailed-3D HotSpot permits spatially
 * varying resistance values whose historical one-sided interface rule can be
 * nonsymmetric, so reject such instances before entering PCG. */
static double maximum_grid_coupling_asymmetry(grid_model_t *model)
{
  int n, i, j;
  int nl = model->n_layers;
  int nr = model->rows;
  int nc = model->cols;
  double mismatch = 0.0;

  for (n = 0; n < nl; n++)
    for (i = 0; i < nr; i++)
      for (j = 0; j < nc; j++) {
        if (j + 1 < nc)
          mismatch = MAX(mismatch, relative_difference(
                STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_EAST),
                STEADY_CONDUCTANCE(model, n, i, j + 1, CONDUCTANCE_WEST)));
        if (i + 1 < nr)
          mismatch = MAX(mismatch, relative_difference(
                STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_SOUTH),
                STEADY_CONDUCTANCE(model, n, i + 1, j, CONDUCTANCE_NORTH)));
        if (n + 1 < nl)
          mismatch = MAX(mismatch, relative_difference(
                STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_BELOW),
                STEADY_CONDUCTANCE(model, n + 1, i, j, CONDUCTANCE_ABOVE)));
      }
  return mismatch;
}

static void build_jacobi_diagonal(grid_model_t *model,
                                  grid_model_vector_t *power,
                                  const double *cap, const double *rhs,
                                  double *diagonal, double *basis,
                                  double *work, double *scaled)
{
  int n, i, j, k;
  int nl = model->n_layers;
  int nr = model->rows;
  int nc = model->cols;
  int grid_nodes = nl * nr * nc;
  int count = steady_node_count(model);
  int spidx = nl - DEFAULT_PACK_LAYERS + LAYER_SP;
  int hsidx = nl - DEFAULT_PACK_LAYERS + LAYER_SINK;
  double cw = model->width / nc;
  double ch = model->height / nr;
  layer_t *layer = model->layers;

  for (n = 0; n < nl; n++)
    for (i = 0; i < nr; i++)
      for (j = 0; j < nc; j++) {
        k = n * nr * nc + i * nc + j;
        diagonal[k] =
          STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_NORTH) +
          STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_SOUTH) +
          STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_EAST) +
          STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_WEST) +
          STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_ABOVE) +
          STEADY_CONDUCTANCE(model, n, i, j, CONDUCTANCE_BELOW);
        if (n == spidx) {
          if (i == 0 || i == nr - 1)
            diagonal[k] += 1.0 / (layer[n].ry / 2.0 + nc * model->pack.r_sp1_y);
          if (j == 0 || j == nc - 1)
            diagonal[k] += 1.0 / (layer[n].rx / 2.0 + nr * model->pack.r_sp1_x);
        } else if (n == hsidx) {
          diagonal[k] += 1.0 / layer[n].rz;
          if (i == 0 || i == nr - 1)
            diagonal[k] += 1.0 / (layer[n].ry / 2.0 + nc * model->pack.r_hs1_y);
          if (j == 0 || j == nc - 1)
            diagonal[k] += 1.0 / (layer[n].rx / 2.0 + nr * model->pack.r_hs1_x);
        } else if (model->config.model_secondary && n == LAYER_SUB) {
          if (i == 0 || i == nr - 1)
            diagonal[k] += 1.0 / (layer[n].ry / 2.0 + nc * model->pack.r_sub1_y);
          if (j == 0 || j == nc - 1)
            diagonal[k] += 1.0 / (layer[n].rx / 2.0 + nr * model->pack.r_sub1_x);
        } else if (model->config.model_secondary && n == LAYER_SOLDER) {
          if (i == 0 || i == nr - 1)
            diagonal[k] += 1.0 / (layer[n].ry / 2.0 + nc * model->pack.r_solder1_y);
          if (j == 0 || j == nc - 1)
            diagonal[k] += 1.0 / (layer[n].rx / 2.0 + nr * model->pack.r_solder1_x);
        } else if (model->config.model_secondary && n == LAYER_PCB) {
          diagonal[k] += 1.0 / (model->config.r_convec_sec *
              model->config.s_pcb * model->config.s_pcb / (cw * ch));
          if (i == 0 || i == nr - 1)
            diagonal[k] += 1.0 / (layer[n].ry / 2.0 + nc * model->pack.r_pcb1_y);
          if (j == 0 || j == nc - 1)
            diagonal[k] += 1.0 / (layer[n].rx / 2.0 + nr * model->pack.r_pcb1_x);
        }
      }

  /* There are only 12 (or 28) package nodes.  Extracting their exact
   * diagonal entries with basis matvecs keeps a single package equation path
   * and is negligible beside hundreds of Krylov iterations. */
  for (k = grid_nodes; k < count; k++) {
    zero_dvector(basis, count);
    basis[k] = 1.0;
    steady_matvec(model, power, cap, rhs, basis, work, scaled);
    diagonal[k] = work[k];
  }

  for (k = 0; k < count; k++)
    if (!(diagonal[k] > 0.0) || !isfinite(diagonal[k]))
      fatal("invalid Jacobi diagonal in detailed-3D PCG solver\n");
}

#if GATE_GS1_ARM == GS1_ARM_GS_NATIVE
static void gs1_gauss_seidel_steady_grid(grid_model_t *model,
                                         grid_model_vector_t *power,
                                         grid_model_vector_t *temp)
{
  double started_s;
  double delta;
  long long first_pass_sweeps = 0;

  started_s = gs1_now();
  set_heuristic_temp(model, power, temp);
  model->gs1_initialization_s += gs1_now() - started_s;
  started_s = gs1_now();
  do {
    delta = single_iteration_steady_grid(model, power, temp);
    first_pass_sweeps++;
  } while (!eq(delta, 0));
  model->gs1_solve_s += gs1_now() - started_s;
  model->gs1_iterations = first_pass_sweeps;
  model->gs1_solution_iterations = first_pass_sweeps;
  model->gs1_native_delta = delta;
}

#endif

#if GATE_GS0 > 0

#define GS0_PROBE_COUNT 8

static void gs0_build_probes(double *probes, int count)
{
  int i;
  for (i = 0; i < count; i++) {
    double sign = (i & 1) ? -1.0 : 1.0;
    probes[0 * count + i] = sign * (double)(1 + i % 7) * 1.0e-12;
    probes[1 * count + i] = -sign * (1.0 + (double)(i % 5) / 8.0) * 1.0e12;
    probes[2 * count + i] = (double)(i % 17 - 8);
    probes[3 * count + i] = 280.0 + (double)(i % 19) / 2.0;
    probes[6 * count + i] = sign * (1.0 + (double)(i % 11) / 4.0) *
      (i % 3 == 0 ? 1.0e-12 : (i % 3 == 1 ? 1.0 : 1.0e12));
    probes[7 * count + i] = (double)((i * 37) % 101 - 50) / 8.0;
  }
  for (i = 0; i < count; i++) {
    probes[4 * count + i] =
      0.25 * probes[0 * count + i] - 1.5 * probes[1 * count + i];
    probes[5 * count + i] =
      -2.0 * probes[2 * count + i] + 0.75 * probes[3 * count + i];
  }
}

static void gs0_export_matrix_free(grid_model_t *model,
                                   grid_model_vector_t *power,
                                   const double *cap, const double *rhs,
                                   const double *diagonal,
                                   double coupling_asymmetry,
                                   int residual_input_unchanged,
                                   double *basis, double *work,
                                   double *scaled)
{
  const char *prefix = gs0_output_prefix();
  int count = steady_node_count(model);
  int input_mutation_count = residual_input_unchanged ? 0 : 1;
  int i, j;
  double *matrix;
  double *probes;
  double *actions;
  double *input_copy;

  if (!prefix)
    return;
  if (!model->config.detailed_3D_used)
    fatal("G-S0 export requires detailed-3D mode\n");
  if (count > GS0_MAX_NODE_COUNT)
    fatal("G-S0 dense export is restricted to at most 2048 nodes\n");

  matrix = (double *)calloc((size_t)count * count, sizeof(double));
  probes = (double *)calloc((size_t)GS0_PROBE_COUNT * count, sizeof(double));
  actions = (double *)calloc((size_t)GS0_PROBE_COUNT * count, sizeof(double));
  input_copy = (double *)calloc(count, sizeof(double));
  if (!matrix || !probes || !actions || !input_copy)
    fatal("memory allocation failed in G-S0 matrix-free export\n");

  for (j = 0; j < count; j++) {
    zero_dvector(basis, count);
    basis[j] = 1.0;
    memcpy(input_copy, basis, (size_t)count * sizeof(double));
    steady_matvec(model, power, cap, rhs, basis, work, scaled);
    if (memcmp(input_copy, basis, (size_t)count * sizeof(double)) != 0)
      input_mutation_count++;
    for (i = 0; i < count; i++)
      matrix[(size_t)i * count + j] = work[i];
  }

  gs0_build_probes(probes, count);
  for (j = 0; j < GS0_PROBE_COUNT; j++) {
    double *input = &probes[(size_t)j * count];
    double *action = &actions[(size_t)j * count];
    memcpy(input_copy, input, (size_t)count * sizeof(double));
    steady_matvec(model, power, cap, rhs, input, action, scaled);
    if (memcmp(input_copy, input, (size_t)count * sizeof(double)) != 0)
      input_mutation_count++;
  }

  gs0_write_doubles(prefix, ".matrix.bin", matrix,
                    (size_t)count * count);
  gs0_write_doubles(prefix, ".rhs.bin", rhs, (size_t)count);
  gs0_write_doubles(prefix, ".diag.bin", diagonal, (size_t)count);
  gs0_write_doubles(prefix, ".probes.bin", probes,
                    (size_t)GS0_PROBE_COUNT * count);
  gs0_write_doubles(prefix, ".probe_actions.bin", actions,
                    (size_t)GS0_PROBE_COUNT * count);
  gs0_write_metadata(model, prefix, "mfpcg", GS0_PROBE_COUNT, count,
                     count + GS0_PROBE_COUNT + 1, input_mutation_count,
                     coupling_asymmetry);
  gs0_write_completion(prefix, ".operator.complete");

  free(matrix);
  free(probes);
  free(actions);
  free(input_copy);
  if (input_mutation_count != 0)
    fatal("steady_matvec modified a G-S0 input vector\n");
}

#endif

static void jacobi_pcg_steady_grid(grid_model_t *model,
                                   grid_model_vector_t *power,
                                   grid_model_vector_t *temp)
{
  int count = steady_node_count(model);
  int i, iteration;
  int converged = 0;
#if GATE_GS0 > 0
  int residual_input_unchanged = 1;
#endif
  const char *termination_status;
  double coupling_asymmetry;
  double rhs_norm, residual_norm;
  double rho, next_rho, alpha, beta, pap;
#if GATE_GS1_ARM == GS1_ARM_MFPCG || \
    GATE_GS1_ARM == GS1_ARM_CSR_PCG || \
    GATE_GS1_ARM == GS1_ARM_STRUCTURED_PCG
  double gs1_started_s = gs1_now();
#endif
  double *cap = (double *)calloc(count, sizeof(double));
  double *rhs = (double *)calloc(count, sizeof(double));
  double *r = (double *)calloc(count, sizeof(double));
  double *z = (double *)calloc(count, sizeof(double));
  double *direction = (double *)calloc(count, sizeof(double));
  double *product = (double *)calloc(count, sizeof(double));
  double *diagonal = (double *)calloc(count, sizeof(double));
  double *probe = (double *)calloc(count, sizeof(double));
  double *scaled = (double *)calloc(count, sizeof(double));
  double *x = temp->cuboid[0][0];
#if GATE_GS1_ARM == GS1_ARM_CSR_PCG
  steady_csr_t csr = {0};
  SuperMatrix explicit_matrix, explicit_rhs_matrix;
  double *explicit_rhs = NULL;
#elif GATE_GS1_ARM == GS1_ARM_STRUCTURED_PCG
  steady_structured_t structured = {0};
  SuperMatrix explicit_matrix, explicit_rhs_matrix;
  double *explicit_rhs = NULL;
#endif
#if GATE_GS1_ARM == GS1_ARM_MFPCG
  model->gs1_solver_vector_bytes = (size_t)count * sizeof(double) * 9;
  model->gs1_operator_bytes = (size_t)model->n_layers * model->rows *
    model->cols * CONDUCTANCE_COUNT * sizeof(double);
#elif GATE_GS1_ARM == GS1_ARM_CSR_PCG
  model->gs1_solver_vector_bytes = (size_t)count * sizeof(double) * 9;
#elif GATE_GS1_ARM == GS1_ARM_STRUCTURED_PCG
  model->gs1_solver_vector_bytes = (size_t)count * sizeof(double) * 9;
#endif
#if GATE_GS0 > 0
  double bitwise_zero = 0.0;
#endif

  /* Steady-only HotSpot normally leaves C uninitialized.  PCG does not use
   * thermal capacitance physically; populate it only as an algebraic scaling
   * that lets slope_fn_grid expose the already implemented residual. */
  if (!model->c_ready)
    populate_C_model_grid(model, NULL);
  if (!cap || !rhs || !r || !z || !direction || !product || !diagonal ||
      !probe || !scaled)
    fatal("memory allocation failed in detailed-3D PCG solver\n");

  build_capacitance_vector(model, cap);
  for (i = 0; i < count; i++)
    if (!(cap[i] > 0.0) || !isfinite(cap[i]))
      fatal("invalid capacitance scaling in detailed-3D PCG solver\n");
  coupling_asymmetry = maximum_grid_coupling_asymmetry(model);
  zero_dvector(probe, count);
  steady_residual(model, power, cap, probe, rhs);
#if GATE_GS0 > 0
  for (i = 0; i < count; i++)
    if (memcmp(&probe[i], &bitwise_zero, sizeof(double)) != 0)
      residual_input_unchanged = 0;
#endif
#if GATE_GS1_ARM == GS1_ARM_CSR_PCG || \
    GATE_GS1_ARM == GS1_ARM_STRUCTURED_PCG
  gs1_started_s = gs1_now();
  explicit_matrix = build_steady_grid_matrix(model);
  explicit_rhs_matrix = build_steady_rhs_vector(model, power, &explicit_rhs);
  if (!explicit_rhs)
    fatal("missing explicit CSR-PCG right-hand side\n");
  for (i = 0; i < count; i++)
    if (relative_difference(rhs[i], explicit_rhs[i]) > 1.0e-12)
      fatal("native and explicit HotSpot right-hand sides differ\n");
#if GATE_GS1_ARM == GS1_ARM_CSR_PCG
  csr_from_superlu(&explicit_matrix, &csr, diagonal);
  model->gs1_operator_bytes =
    ((size_t)csr.count + 1 + (size_t)csr.nnz * 2) * sizeof(int_t) +
    (size_t)csr.nnz * sizeof(double);
#else
  structured_from_superlu(&explicit_matrix, &structured, diagonal,
                          model->rows, model->cols, model->n_layers);
  model->gs1_operator_bytes =
    (size_t)structured.grid_count * sizeof(double) * 3 +
    (size_t)structured.package_edge_count *
      (sizeof(int_t) * 2 + sizeof(double));
#endif
  gs1_audit_explicit_representation(model, power, cap, rhs, diagonal,
#if GATE_GS1_ARM == GS1_ARM_CSR_PCG
                                    &csr,
#else
                                    &structured,
#endif
                                    direction, r, product, scaled);
  SUPERLU_FREE(explicit_rhs);
  Destroy_CompCol_Matrix(&explicit_matrix);
  Destroy_SuperMatrix_Store(&explicit_rhs_matrix);
  model->gs1_assembly_s += gs1_now() - gs1_started_s;
#else
  build_jacobi_diagonal(model, power, cap, rhs,
                        diagonal, direction, product, scaled);
#endif
#if GATE_GS1_ARM == GS1_ARM_MFPCG || \
    GATE_GS1_ARM == GS1_ARM_CSR_PCG || \
    GATE_GS1_ARM == GS1_ARM_STRUCTURED_PCG
  model->gs1_factor_s += gs1_now() - gs1_started_s;
  {
    const char *diagonal_prefix = getenv("HOTSPOT_GS1_DIAGONAL_PREFIX");
    if (diagonal_prefix && diagonal_prefix[0]) {
      gs1_write_bytes(diagonal_prefix, ".jacobi.bin", diagonal,
                      sizeof(double), (size_t)count);
      gs1_write_completion(diagonal_prefix, ".jacobi.complete");
    }
  }
#endif
#if GATE_GS0 < 1
  /* Export the G-S1 Jacobi receipt first, but still reject the operator before
   * the first Krylov step.  This lets the independent validator distinguish a
   * scientific supported-domain STOP from a missing-artifact wiring failure. */
  if (coupling_asymmetry > PCG_SYMMETRY_TOLERANCE)
    fatal("detailed-3D steady operator has asymmetric grid couplings\n");
#endif
#if GATE_GS0 > 0
  gs0_export_matrix_free(model, power, cap, rhs, diagonal,
                         coupling_asymmetry, residual_input_unchanged,
                         direction, product, scaled);
  if (!residual_input_unchanged)
    fatal("steady residual modified its G-S0 zero input vector\n");
  if (coupling_asymmetry > PCG_SYMMETRY_TOLERANCE) {
    const char *prefix = gs0_output_prefix();
    if (prefix)
      gs0_write_completion(prefix, ".rejected_before_pcg");
    fatal("detailed-3D steady operator has asymmetric grid couplings\n");
  }
#endif

  #if GATE_GS1_ARM == GS1_ARM_MFPCG || \
      GATE_GS1_ARM == GS1_ARM_CSR_PCG || \
      GATE_GS1_ARM == GS1_ARM_STRUCTURED_PCG
  gs1_started_s = gs1_now();
  #endif
  set_heuristic_temp(model, power, temp);
  #if GATE_GS1_ARM == GS1_ARM_MFPCG || \
      GATE_GS1_ARM == GS1_ARM_CSR_PCG || \
      GATE_GS1_ARM == GS1_ARM_STRUCTURED_PCG
  model->gs1_initialization_s += gs1_now() - gs1_started_s;
  gs1_started_s = gs1_now();
  #endif
  steady_residual(model, power, cap, x, r);
  rhs_norm = sqrt(dot_product(rhs, rhs, count));
  if (!(rhs_norm > 0.0) || !isfinite(rhs_norm))
    fatal("invalid right-hand side in detailed-3D PCG solver\n");
  residual_norm = sqrt(dot_product(r, r, count));
  if (!isfinite(residual_norm))
    fatal("invalid initial residual in detailed-3D PCG solver\n");
  iteration = 0;
  if (residual_norm <= PCG_RELATIVE_RESIDUAL * rhs_norm) {
    converged = 1;
  } else {
    for (i = 0; i < count; i++) {
      z[i] = r[i] / diagonal[i];
      direction[i] = z[i];
    }
    rho = dot_product(r, z, count);
    if (!(rho > 0.0) || !isfinite(rho))
      fatal("invalid initial residual in detailed-3D PCG solver\n");

    for (iteration = 0; iteration < PCG_MAX_ITERATIONS; iteration++) {
#if GATE_GS1_ARM == GS1_ARM_CSR_PCG
      model->gs1_operator_applications++;
      csr_matvec(&csr, direction, product);
#elif GATE_GS1_ARM == GS1_ARM_STRUCTURED_PCG
      model->gs1_operator_applications++;
      structured_matvec(&structured, diagonal, direction, product,
                        model->rows, model->cols);
#else
      steady_matvec(model, power, cap, rhs, direction, product, scaled);
#endif
      pap = dot_product(direction, product, count);
      if (!(pap > 0.0) || !isfinite(pap))
        fatal("detailed-3D PCG encountered a non-positive direction\n");
      alpha = rho / pap;
      for (i = 0; i < count; i++) {
        x[i] += alpha * direction[i];
        r[i] -= alpha * product[i];
      }
      residual_norm = sqrt(dot_product(r, r, count));
      if (!isfinite(residual_norm))
        fatal("invalid residual in detailed-3D PCG solver\n");
      if (residual_norm <= PCG_RELATIVE_RESIDUAL * rhs_norm) {
        /* Validate apparent convergence with the physical residual.  If the
         * recursively updated residual has drifted, replace it and restart
         * the Krylov direction instead of accepting or rejecting the same
         * otherwise valid iterate. */
        steady_residual(model, power, cap, x, r);
        residual_norm = sqrt(dot_product(r, r, count));
        if (!isfinite(residual_norm))
          fatal("invalid physical residual in detailed-3D PCG solver\n");
        if (residual_norm <= PCG_RELATIVE_RESIDUAL * rhs_norm) {
          converged = 1;
          iteration++;
          break;
        }
#if GATE_GS1_ARM == GS1_ARM_MFPCG || \
    GATE_GS1_ARM == GS1_ARM_CSR_PCG || \
    GATE_GS1_ARM == GS1_ARM_STRUCTURED_PCG
        model->gs1_residual_replacements++;
#endif
        for (i = 0; i < count; i++) {
          z[i] = r[i] / diagonal[i];
          direction[i] = z[i];
        }
        rho = dot_product(r, z, count);
        if (!(rho > 0.0) || !isfinite(rho))
          fatal("invalid physical residual in detailed-3D PCG solver\n");
        continue;
      }
      for (i = 0; i < count; i++)
        z[i] = r[i] / diagonal[i];
      next_rho = dot_product(r, z, count);
      if (!(next_rho > 0.0) || !isfinite(next_rho))
        fatal("invalid residual in detailed-3D PCG solver\n");
      beta = next_rho / rho;
      for (i = 0; i < count; i++)
        direction[i] = z[i] + beta * direction[i];
      rho = next_rho;
    }
  }

  if (!converged) {
    steady_residual(model, power, cap, x, r);
    residual_norm = sqrt(dot_product(r, r, count));
    if (!isfinite(residual_norm))
      fatal("invalid final physical residual in detailed-3D PCG solver\n");
    converged = residual_norm <= PCG_RELATIVE_RESIDUAL * rhs_norm;
  }
  termination_status = converged ? "converged" : "iteration_limit";
#if GATE_GS1_ARM == GS1_ARM_MFPCG || \
    GATE_GS1_ARM == GS1_ARM_CSR_PCG || \
    GATE_GS1_ARM == GS1_ARM_STRUCTURED_PCG
  model->gs1_iterations = iteration;
  model->gs1_solution_iterations = iteration;
  model->gs1_relative_residual = residual_norm / rhs_norm;
#endif
#if VERBOSE > 0
  fprintf(stdout,
          "detailed-3D PCG termination: status=%s iterations=%d limit=%d relative_residual=%.9e\n",
          termination_status, iteration, PCG_MAX_ITERATIONS,
          residual_norm / rhs_norm);
#endif
  if (!converged)
    fatal("detailed-3D PCG reached its iteration limit\n");

#if GATE_GS0 > 0
  gs0_write_state(model, temp->cuboid[0][0]);
#endif

  free(cap);
  free(rhs);
  free(r);
  free(z);
  free(direction);
  free(product);
  free(diagonal);
  free(probe);
  free(scaled);
#if GATE_GS1_ARM == GS1_ARM_CSR_PCG
  free_steady_csr(&csr);
#elif GATE_GS1_ARM == GS1_ARM_STRUCTURED_PCG
  free_steady_structured(&structured);
#endif
#if GATE_GS1_ARM == GS1_ARM_MFPCG || \
    GATE_GS1_ARM == GS1_ARM_CSR_PCG || \
    GATE_GS1_ARM == GS1_ARM_STRUCTURED_PCG
  model->gs1_solve_s += gs1_now() - gs1_started_s;
#endif
}

#endif

void compute_temp_grid(grid_model_t *model, double *power, double *temp, double time_elapsed)
{
  double t, h, new_h;
  int extra_nodes;
  grid_model_vector_t *p;
#if VERBOSE > 1
  unsigned int i = 0;
#endif

  if (model->config.model_secondary)
    extra_nodes = EXTRA + EXTRA_SEC;
  else
    extra_nodes = EXTRA;

  if (!model->r_ready || !model->c_ready)
    fatal("grid model not ready\n");

  p = new_grid_model_vector(model);

  /* package nodes' power numbers	*/
  set_internal_power_grid(model, power);

  /* map the block power/temp numbers to the grid	*/
  xlate_vector_b2g(model, power, p, V_POWER);

  /* if temp is NULL, re-use the temperature from the
   * last call. otherwise, translate afresh and remember 
   * the grid and block temperature arrays for future use
   */
  if (temp != NULL) {
      xlate_vector_b2g(model, temp, model->last_trans, V_TEMP);
      model->last_temp = temp;
  }

  /* Obtain temp at time (t+time_elapsed). 
   * Instead of getting the temperature at t+time_elapsed directly, we
   * do it in multiple steps with the correct step size at each time 
   * provided by rk4. 
   */
  for (t = 0, new_h = MIN_STEP; t < time_elapsed && new_h >= MIN_STEP*DELTA; t+=h) {
      h = new_h;
      /* pass the entire grid and the tail of package nodes 
       * as a 1-d array
       */
      new_h = rk4(model, model->last_trans->cuboid[0][0],  p, 
                  /* array size = grid size + EXTRA	*/
                  model->rows * model->cols * model->n_layers + extra_nodes, &h,
                  model->last_trans->cuboid[0][0], 
                  /* the slope function callback is typecast accordingly */
                  (slope_fn_ptr) slope_fn_grid);
      new_h = MIN(new_h, time_elapsed-t-h);
#if VERBOSE > 1
      i++;
#endif	
  }

#if VERBOSE > 1
  fprintf(stdout, "no. of rk4 calls during compute_temp: %d\n", i+1);
#endif

  /* map the temperature numbers back	*/
  xlate_temp_g2b(model, model->last_temp, model->last_trans);

  free_grid_model_vector(p);
}

/* debug print	*/
void debug_print_blist(blist_t *head, flp_t *flp)
{
  blist_t *ptr;
  fprintf(stdout, "printing blist information...\n");
  for(ptr = head; ptr; ptr = ptr->next) {
      fprintf(stdout, "unit: %s\n", flp->units[ptr->idx].name);
      fprintf(stdout, "occupancy: %f\n", ptr->occupancy);
  }
}

void debug_print_glist(glist_t *array, flp_t *flp)
{
  int i;
  fprintf(stdout, "printing glist information...\n");
  for(i=0; i < flp->n_units; i++)
    fprintf(stdout, "unit: %s\tstartx: %d\tendx: %d\tstarty: %d\tendy: %d\n",
            flp->units[i].name, array[i].j1, array[i].j2, array[i].i1, array[i].i2);
}
/*BU_3D: Debugging Functions */
void debug_print_layer_det3D(grid_model_t *model, layer_t *layer)
{
  int i, j;
  for(i=0; i < model->rows; i++)
    for(j=0; j < model->cols; j++) {
        fprintf(stdout, "%d,%d,%lf,%lf,%lf\n", i, j,layer->b2gmap[i][j]->rx,layer->b2gmap[i][j]->ry,layer->b2gmap[i][j]->rz);//row,col,rx,ry,rz
    }

}
/*end->BU_3D*/

void debug_print_layer(grid_model_t *model, layer_t *layer)
{
  int i, j;
  fprintf(stdout, "printing layer information...\n");
  fprintf(stdout, "no: %d\n", layer->no);
  fprintf(stdout, "has_lateral: %d\n", layer->has_lateral);
  fprintf(stdout, "has_power: %d\n", layer->has_power);
  fprintf(stdout, "k: %f\n", layer->k);
  fprintf(stdout, "thickness: %f\n", layer->thickness);
  fprintf(stdout, "sp: %f\n", layer->sp);
  fprintf(stdout, "rx: %f\try: %f\trz: %f\tc: %f\n", 
          layer->rx, layer->ry, layer->rz, layer->c);

  fprintf(stdout, "printing b2gmap information...\n");
  for(i=0; i < model->rows; i++)
    for(j=0; j < model->cols; j++) {
        fprintf(stdout, "row: %d, col: %d\n", i, j);
        debug_print_blist(layer->b2gmap[i][j], layer->flp);
    }

  fprintf(stdout, "printing g2bmap information...\n");
  debug_print_glist(layer->g2bmap, layer->flp);
}

//BU_3D added test_b2gmap function
/* test the block-grid map data structure	*/
void test_b2gmap(grid_model_t *model, layer_t *layer)
{
  int i, j;
  blist_t *ptr;
  double sum;

  /* a correctly formed b2gmap should have the 
   * sum of occupancies in each linked list
   * to be equal to 1.0
   */
  for (i=0; i < model->rows; i++)
    for(j=0; j < model->cols; j++) {
        sum = 0.0;
        for(ptr = layer->b2gmap[i][j]; ptr; ptr = ptr->next)
          sum += ptr->occupancy;
        if (!eq(floor(sum*1e5 + 0.5)/1e5, 1.0)) {
            fprintf(stdout, "i: %d\tj: %d\n", i, j);
            debug_print_blist(layer->b2gmap[i][j], layer->flp);
            fatal("erroneous b2gmap data structure. invalid floorplan?\n");
        }	
    }
}//end->BU_3D

void debug_print_grid_model_vector(grid_model_t *model, grid_model_vector_t *v, int nl, int nr, int nc)
{
  int n;
  int extra_nodes;

  if (model->config.model_secondary)
    extra_nodes = EXTRA + EXTRA_SEC;
  else
    extra_nodes = EXTRA;

  fprintf(stdout, "printing cuboid information...\n");
  for(n=0; n < nl; n++)
    dump_dmatrix(v->cuboid[n], nr, nc);
  fprintf(stdout, "printing extra information...\n");
  dump_dvector(v->extra, extra_nodes);
}

void debug_print_grid(grid_model_t *model)
{
  int i;
  int extra_nodes;

  if (model->config.model_secondary)
    extra_nodes = EXTRA + EXTRA_SEC;
  else
    extra_nodes = EXTRA;

  fprintf(stdout, "printing grid model information...\n");
  fprintf(stdout, "rows: %d\n", model->rows);
  fprintf(stdout, "cols: %d\n", model->cols);
  fprintf(stdout, "width: %f\n", model->width);
  fprintf(stdout, "height: %f\n", model->height);

  debug_print_package_RC(&model->pack);

  fprintf(stdout, "total_n_blocks: %d\n", model->total_n_blocks);
  fprintf(stdout, "map_mode: %d\n", model->map_mode);
  fprintf(stdout, "r_ready: %d\n", model->r_ready);
  fprintf(stdout, "c_ready: %d\n", model->c_ready);
  fprintf(stdout, "has_lcf: %d\n", model->has_lcf);

  fprintf(stdout, "printing last_steady information...\n");
  debug_print_grid_model_vector(model, model->last_steady, model->n_layers, 
                                model->rows, model->cols);

  fprintf(stdout, "printing last_trans information...\n");
  debug_print_grid_model_vector(model, model->last_trans, model->n_layers,
                                model->rows, model->cols);

  fprintf(stdout, "printing last_temp information...\n");
  if (model->last_temp)
    dump_dvector(model->last_temp, model->total_n_blocks + extra_nodes);
  else
    fprintf(stdout, "(null)\n");

  for(i=0; i < model->n_layers; i++)
    debug_print_layer(model, &model->layers[i]);

  fprintf(stdout, "base_n_units: %d\n", model->base_n_units);
}

#if SUPERLU > 0
SuperMatrix build_steady_grid_matrix(grid_model_t *model)
{
  SuperMatrix A;
  double   *a;
  int      *asub, *xa;
  double   *cooV;
  int      *cooX, *cooY;

  int      i, j, l, m, n, nnz;
  double   dia_val;
  int      curidx, grididx;
  int      xoffset, yoffset;
  double   Rn, Rs, Rw, Re, Ra, Rb;
  double   R_temp;

  /* shortcuts	*/
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int spidx, hsidx;
  int subidx = LAYER_SUB;
  int solderidx = LAYER_SOLDER;
  int pcbidx = LAYER_PCB;
  int model_secondary = model->config.model_secondary;
  double cw = model->width / model->cols;
  double ch = model->height / model->rows;
  layer_t *lyr = model->layers;
  package_RC_t *pk = &model->pack;

  spidx = nl - DEFAULT_PACK_LAYERS + LAYER_SP;
  hsidx = nl - DEFAULT_PACK_LAYERS + LAYER_SINK;
  /* Initialize matrix A. */
  if(model_secondary){
      m = n = (nl*nc*nr + EXTRA + EXTRA_SEC);
      /* Num of non-zeros
       * each layer Five diagonal  : c*r, c*(r-1), (c-1)*r, c*(r-1), (c-1)*r
       * vertical connections      : 2*(nl-1)*nr*nc
       * layer peripheral          : 5 * 2 * (2*nr+2*nc)
       * diaganal for pkg nodes    : EXTRA + EXTRA_SEC
       * between prim pkg nodes    : 2 * 2 * 4 
       * between sec  pkg nodes    : 3 * 2 * 4 
       */
      nnz = 5*nr*nc-2*(nr+nc);
      nnz *= nl;
      nnz += 2*(nl-1)*nr*nc;
      nnz += 5*2*2*(nr+nc);
      nnz += EXTRA + EXTRA_SEC;
      nnz += 16 + 24;
  }
  else{
      m = n = (nl*nc*nr + EXTRA);
      /* Num of non-zeros
       * each layer Five diagonal  : c*r, c*(r-1), (c-1)*r, c*(r-1), (c-1)*r
       * vertical connections      : 2*(nl-1)*nr*nc
       * sp and hs layer peripheral: 2 * 2 * (2*nr+2*nc)
       * diaganal for pkg nodes    : EXTRA
       * between pkg nodes         : 2 * 2 * 4 
       */
      nnz = 5*nr*nc-2*(nr+nc);
      nnz *= nl;
      nnz += 2*(nl-1)*nr*nc;
      nnz += 2*2*2*(nr+nc);
      nnz += EXTRA + 16;
  }

  if ( !(cooV = doubleMalloc(nnz)) ) fatal("Malloc fails for cooV[].\n");
  if ( !(cooX = intMalloc(nnz)) ) fatal("Malloc fails for cooX[].\n");
  if ( !(cooY = intMalloc(nnz)) ) fatal("Malloc fails for cooY[].\n");

  if ( !(a = doubleMalloc(nnz)) ) fatal("Malloc fails for a[].\n");
  if ( !(asub = intMalloc(nnz)) ) fatal("Malloc fails for asub[].\n");
  if ( !(xa = intMalloc(n+1)) ) fatal("Malloc fails for xa[].\n");

  curidx = 0;
  for(l=0; l<nl; l++){
      xoffset = l*nr*nc;
      yoffset = l*nr*nc;
      for(i=0; i<nr; i++)
        for(j=0; j<nc; j++){
            grididx = i*nc + j;

            if(model->config.detailed_3D_used == 1){
                if(j > 0)    Rw = find_res_3D(l,i,j-1,model,1); else Rw = LARGENUM;
                if(j < nc-1) Re = find_res_3D(l,i,j+1,model,1); else Re = LARGENUM;
                if(i > 0)    Rn = find_res_3D(l,i-1,j,model,2); else Rn = LARGENUM;
                if(i < nr-1) Rs = find_res_3D(l,i+1,j,model,2); else Rs = LARGENUM;
                if(l > 0)    Ra = find_res_3D(l-1,i,j,model,3); else Ra = LARGENUM;
                if(l < nl-1) Rb = find_res_3D(l,i,j,model,3);   else Rb = LARGENUM;
            }
            else{
                if(j > 0)    Rw = lyr[l].rx;   else Rw = LARGENUM;
                if(j < nc-1) Re = lyr[l].rx;   else Re = LARGENUM;
                if(i > 0)    Rn = lyr[l].ry;   else Rn = LARGENUM;
                if(i < nr-1) Rs = lyr[l].ry;   else Rs = LARGENUM;
                if(l > 0)    Ra = lyr[l-1].rz; else Ra = LARGENUM;
                if(l < nl-1) Rb = lyr[l].rz;   else Rb = LARGENUM;
            }

            dia_val = 0;
            if(0 == grididx){//top left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = -1.0/Re; curidx++; dia_val += 1.0/Re;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = -1.0/Rs; curidx++; dia_val += 1.0/Rs;

                if(l<nl-1){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nr*nc+yoffset;
                    cooV[curidx] = -1.0/Rb; curidx++; dia_val += 1.0/Rb;
                }
                if(l>0){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nr*nc+yoffset;
                    cooV[curidx] = -1.0/Ra; curidx++; dia_val += 1.0/Ra;
                }

                if(l == spidx){
                    /* northern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_sp1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SP_N;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                    /* western boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_sp1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SP_W;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if(l == hsidx){
                    /* northern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_hs1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SINK_C_N;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                    /* western boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_hs1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SINK_C_W;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

                    /* all nodes are connected to the ambient	*/
                    R_temp = lyr[l].rz;
                    dia_val += 1.0/R_temp;
                }
                else if ((l==subidx) && model_secondary) {
                    /* northern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_sub1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SUB_N;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                    /* western boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_sub1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SUB_W;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if ((l==solderidx) && model_secondary) {
                    /* northern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_solder1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SOLDER_N;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                    /* western boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_solder1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SOLDER_W;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if ((l==pcbidx) && model_secondary) {
                    /* northern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_pcb1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + PCB_C_N;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                    /* western boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_pcb1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + PCB_C_W;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

                    /* all nodes are connected to the ambient	*/
                    R_temp = model->config.r_convec_sec * \
                             (model->config.s_pcb * model->config.s_pcb) / (cw * ch);
                    dia_val += 1.0/R_temp;
                }

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = dia_val; curidx++;
            }
            else if((nc-1) == grididx ){//top right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = -1.0/Rw; curidx++; dia_val += 1.0/Rw;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = -1.0/Rs; curidx++; dia_val += 1.0/Rs;

                if(l<nl-1){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nr*nc+yoffset;
                    cooV[curidx] = -1.0/Rb; curidx++; dia_val += 1.0/Rb;
                }
                if(l>0){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nr*nc+yoffset;
                    cooV[curidx] = -1.0/Ra; curidx++; dia_val += 1.0/Ra;
                }

                if(l == spidx){
                    /* northern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_sp1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SP_N;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                    /* eastern boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_sp1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SP_E;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if(l == hsidx){
                    /* northern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_hs1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SINK_C_N;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                    /* eastern boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_hs1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SINK_C_E;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

                    /* all nodes are connected to the ambient	*/
                    R_temp = lyr[l].rz;
                    dia_val += 1.0/R_temp;
                }
                else if ((l==subidx) && model_secondary) {
                    /* northern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_sub1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SUB_N;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                    /* eastern boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_sub1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SUB_E;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if ((l==solderidx) && model_secondary) {
                    /* northern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_solder1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SOLDER_N;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                    /* eastern boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_solder1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SOLDER_E;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if ((l==pcbidx) && model_secondary) {
                    /* northern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_pcb1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + PCB_C_N;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                    /* eastern boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_pcb1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + PCB_C_E;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

                    /* all nodes are connected to the ambient	*/
                    R_temp = model->config.r_convec_sec * \
                             (model->config.s_pcb * model->config.s_pcb) / (cw * ch);
                    dia_val += 1.0/R_temp;
                }

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = dia_val; curidx++;
            }
            else if((nr*nc-nc) == grididx){//bottom left corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = -1.0/Re; curidx++; dia_val += 1.0/Re;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = -1.0/Rn; curidx++; dia_val += 1.0/Rn;

                if(l<nl-1){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nr*nc+yoffset;
                    cooV[curidx] = -1.0/Rb; curidx++; dia_val += 1.0/Rb;
                }
                if(l>0){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nr*nc+yoffset;
                    cooV[curidx] = -1.0/Ra; curidx++; dia_val += 1.0/Ra;
                }

                if(l == spidx){
                    /* southern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_sp1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SP_S;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                    /* western boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_sp1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SP_W;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if(l == hsidx){
                    /* southern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_hs1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SINK_C_S;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                    /* western boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_hs1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SINK_C_W;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

                    /* all nodes are connected to the ambient	*/
                    R_temp = lyr[l].rz;
                    dia_val += 1.0/R_temp;
                }
                else if ((l==subidx) && model_secondary) {
                    /* southern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_sub1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SUB_S;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                    /* western boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_sub1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SUB_W;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if ((l==solderidx) && model_secondary) {
                    /* southern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_solder1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SOLDER_S;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                    /* western boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_solder1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SOLDER_W;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if ((l==pcbidx) && model_secondary) {
                    /* southern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_pcb1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + PCB_C_S;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                    /* western boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_pcb1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + PCB_C_W;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

                    /* all nodes are connected to the ambient	*/
                    R_temp = model->config.r_convec_sec * \
                             (model->config.s_pcb * model->config.s_pcb) / (cw * ch);
                    dia_val += 1.0/R_temp;
                }

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = dia_val; curidx++;
            }
            else if((nr*nc-1) == grididx){//bottom right corner
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = -1.0/Rw; curidx++; dia_val += 1.0/Rw;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = -1.0/Rn; curidx++; dia_val += 1.0/Rn;

                if(l<nl-1){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nr*nc+yoffset;
                    cooV[curidx] = -1.0/Rb; curidx++; dia_val += 1.0/Rb;
                }
                if(l>0){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nr*nc+yoffset;
                    cooV[curidx] = -1.0/Ra; curidx++; dia_val += 1.0/Ra;
                }

                if(l == spidx){
                    /* southern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_sp1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SP_S;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                    /* eastern boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_sp1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SP_E;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if(l == hsidx){
                    /* southern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_hs1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SINK_C_S;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                    /* eastern boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_hs1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SINK_C_E;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

                    /* all nodes are connected to the ambient	*/
                    R_temp = lyr[l].rz;
                    dia_val += 1.0/R_temp;
                }
                else if ((l==subidx) && model_secondary) {
                    /* southern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_sub1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SUB_S;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                    /* eastern boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_sub1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SUB_E;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if ((l==solderidx) && model_secondary) {
                    /* southern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_solder1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SOLDER_S;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                    /* eastern boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_solder1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SOLDER_E;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if ((l==pcbidx) && model_secondary) {
                    /* southern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_pcb1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + PCB_C_S;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                    /* eastern boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_pcb1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + PCB_C_E;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

                    /* all nodes are connected to the ambient	*/
                    R_temp = model->config.r_convec_sec * \
                             (model->config.s_pcb * model->config.s_pcb) / (cw * ch);
                    dia_val += 1.0/R_temp;
                }

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = dia_val; curidx++;
            }
            else if((grididx > 0) && (grididx < nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = -1.0/Rw; curidx++; dia_val += 1.0/Rw;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = -1.0/Re; curidx++; dia_val += 1.0/Re;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = -1.0/Rs; curidx++; dia_val += 1.0/Rs;

                if(l<nl-1){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nr*nc+yoffset;
                    cooV[curidx] = -1.0/Rb; curidx++; dia_val += 1.0/Rb;
                }
                if(l>0){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nr*nc+yoffset;
                    cooV[curidx] = -1.0/Ra; curidx++; dia_val += 1.0/Ra;
                }

                if(l == spidx){
                    /* northern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_sp1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SP_N;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if(l == hsidx){
                    /* northern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_hs1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SINK_C_N;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

                    /* all nodes are connected to the ambient	*/
                    R_temp = lyr[l].rz;
                    dia_val += 1.0/R_temp;
                }
                else if ((l==subidx) && model_secondary) {
                    /* northern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_sub1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SUB_N;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if ((l==solderidx) && model_secondary) {
                    /* northern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_solder1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SOLDER_N;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if ((l==pcbidx) && model_secondary) {
                    /* northern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_pcb1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + PCB_C_N;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

                    /* all nodes are connected to the ambient	*/
                    R_temp = model->config.r_convec_sec * \
                             (model->config.s_pcb * model->config.s_pcb) / (cw * ch);
                    dia_val += 1.0/R_temp;
                }

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = dia_val; curidx++;
            }
            else if((grididx > nr*nc-nc) && (grididx < nr*nc-1)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = -1.0/Rw; curidx++; dia_val += 1.0/Rw;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = -1.0/Re; curidx++; dia_val += 1.0/Re;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = -1.0/Rn; curidx++; dia_val += 1.0/Rn;

                if(l<nl-1){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nr*nc+yoffset;
                    cooV[curidx] = -1.0/Rb; curidx++; dia_val += 1.0/Rb;
                }
                if(l>0){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nr*nc+yoffset;
                    cooV[curidx] = -1.0/Ra; curidx++; dia_val += 1.0/Ra;
                }

                if(l == spidx){
                    /* southern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_sp1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SP_S;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if(l == hsidx){
                    /* southern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_hs1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SINK_C_S;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

                    /* all nodes are connected to the ambient	*/
                    R_temp = lyr[l].rz;
                    dia_val += 1.0/R_temp;
                }
                else if ((l==subidx) && model_secondary) {
                    /* southern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_sub1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SUB_S;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if ((l==solderidx) && model_secondary) {
                    /* southern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_solder1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SOLDER_S;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if ((l==pcbidx) && model_secondary) {
                    /* southern boundary - edge cell has half the ry	*/
                    R_temp = lyr[l].ry/2.0 + nc*model->pack.r_pcb1_y;
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + PCB_C_S;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

                    /* all nodes are connected to the ambient	*/
                    R_temp = model->config.r_convec_sec * \
                             (model->config.s_pcb * model->config.s_pcb) / (cw * ch);
                    dia_val += 1.0/R_temp;
                }

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = dia_val; curidx++;
            }
            else if(0 == (grididx%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = -1.0/Rn; curidx++; dia_val += 1.0/Rn;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = -1.0/Rs; curidx++; dia_val += 1.0/Rs;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = -1.0/Re; curidx++; dia_val += 1.0/Re;

                if(l<nl-1){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nr*nc+yoffset;
                    cooV[curidx] = -1.0/Rb; curidx++; dia_val += 1.0/Rb;
                }
                if(l>0){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nr*nc+yoffset;
                    cooV[curidx] = -1.0/Ra; curidx++; dia_val += 1.0/Ra;
                }

                if(l == spidx){
                    /* western boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_sp1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SP_W;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if(l == hsidx){
                    /* western boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_hs1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SINK_C_W;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

                    /* all nodes are connected to the ambient	*/
                    R_temp = lyr[l].rz;
                    dia_val += 1.0/R_temp;
                }
                else if ((l==subidx) && model_secondary) {
                    /* western boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_sub1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SUB_W;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if ((l==solderidx) && model_secondary) {
                    /* western boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_solder1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SOLDER_W;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if ((l==pcbidx) && model_secondary) {
                    /* western boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_pcb1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + PCB_C_W;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

                    /* all nodes are connected to the ambient	*/
                    R_temp = model->config.r_convec_sec * \
                             (model->config.s_pcb * model->config.s_pcb) / (cw * ch);
                    dia_val += 1.0/R_temp;
                }

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = dia_val; curidx++;
            }
            else if(0 == ((grididx+1)%nc)){
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = -1.0/Rn; curidx++; dia_val += 1.0/Rn;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = -1.0/Rs; curidx++; dia_val += 1.0/Rs;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = -1.0/Rw; curidx++; dia_val += 1.0/Rw;

                if(l<nl-1){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nr*nc+yoffset;
                    cooV[curidx] = -1.0/Rb; curidx++; dia_val += 1.0/Rb;
                }
                if(l>0){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nr*nc+yoffset;
                    cooV[curidx] = -1.0/Ra; curidx++; dia_val += 1.0/Ra;
                }

                if(l == spidx){
                    /* eastern boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_sp1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SP_E;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if(l == hsidx){
                    /* eastern boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_hs1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SINK_C_E;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

                    /* all nodes are connected to the ambient	*/
                    R_temp = lyr[l].rz;
                    dia_val += 1.0/R_temp;
                }
                else if ((l==subidx) && model_secondary) {
                    /* eastern boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_sub1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SUB_E;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if ((l==solderidx) && model_secondary) {
                    /* eastern boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_solder1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + SOLDER_E;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
                }
                else if ((l==pcbidx) && model_secondary) {
                    /* eastern boundary	- edge cell has half the rx		*/
                    R_temp = lyr[l].rx/2.0 + nr*model->pack.r_pcb1_x; 
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = nl*nr*nc + PCB_C_E;
                    cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

                    /* all nodes are connected to the ambient	*/
                    R_temp = model->config.r_convec_sec * \
                             (model->config.s_pcb * model->config.s_pcb) / (cw * ch);
                    dia_val += 1.0/R_temp;
                }

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = dia_val; curidx++;
            }
            else{
                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nc+yoffset;
                cooV[curidx] = -1.0/Rn; curidx++; dia_val += 1.0/Rn;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nc+yoffset;
                cooV[curidx] = -1.0/Rs; curidx++; dia_val += 1.0/Rs;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-1+yoffset;
                cooV[curidx] = -1.0/Rw; curidx++; dia_val += 1.0/Rw;

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+1+yoffset;
                cooV[curidx] = -1.0/Re; curidx++; dia_val += 1.0/Re;

                if(l<nl-1){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+nr*nc+yoffset;
                    cooV[curidx] = -1.0/Rb; curidx++; dia_val += 1.0/Rb;
                }
                if(l>0){
                    cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx-nr*nc+yoffset;
                    cooV[curidx] = -1.0/Ra; curidx++; dia_val += 1.0/Ra;
                }

                if(l == hsidx){
                    /* all nodes are connected to the ambient	*/
                    R_temp = lyr[l].rz;
                    dia_val += 1.0/R_temp;
                }
                else if ((l==pcbidx) && model_secondary) {
                    /* all nodes are connected to the ambient	*/
                    R_temp = model->config.r_convec_sec * \
                             (model->config.s_pcb * model->config.s_pcb) / (cw * ch);
                    dia_val += 1.0/R_temp;
                }

                cooX[curidx] = grididx+xoffset; cooY[curidx] = grididx+yoffset;
                cooV[curidx] = dia_val; curidx++;
            }
        }
  }

  /* Package nodes */
  /* sink outer north/south	*/
  xoffset = nl*nr*nc;
  yoffset = nl*nr*nc;
  dia_val = 0;
  R_temp = pk->r_hs2_y + pk->r_hs; 
  cooX[curidx] = SINK_N+xoffset; cooY[curidx] = SINK_C_N+yoffset;
  cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
  R_temp = pk->r_hs_per + pk->r_amb_per; dia_val += 1.0/R_temp;
  cooX[curidx] = SINK_N+xoffset; cooY[curidx] = SINK_N+yoffset;
  cooV[curidx] = dia_val; curidx++;

  dia_val = 0;
  R_temp = pk->r_hs2_y + pk->r_hs; 
  cooX[curidx] = SINK_S+xoffset; cooY[curidx] = SINK_C_S+yoffset;
  cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
  R_temp = pk->r_hs_per + pk->r_amb_per; dia_val += 1.0/R_temp;
  cooX[curidx] = SINK_S+xoffset; cooY[curidx] = SINK_S+yoffset;
  cooV[curidx] = dia_val; curidx++;

  /* sink outer west/east	*/
  dia_val = 0;
  R_temp = pk->r_hs2_x + pk->r_hs; 
  cooX[curidx] = SINK_W+xoffset; cooY[curidx] = SINK_C_W+yoffset;
  cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
  R_temp = pk->r_hs_per + pk->r_amb_per; dia_val += 1.0/R_temp;
  cooX[curidx] = SINK_W+xoffset; cooY[curidx] = SINK_W+yoffset;
  cooV[curidx] = dia_val; curidx++;

  dia_val = 0;
  R_temp = pk->r_hs2_x + pk->r_hs; 
  cooX[curidx] = SINK_E+xoffset; cooY[curidx] = SINK_C_E+yoffset;
  cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
  R_temp = pk->r_hs_per + pk->r_amb_per; dia_val += 1.0/R_temp;
  cooX[curidx] = SINK_E+xoffset; cooY[curidx] = SINK_E+yoffset;
  cooV[curidx] = dia_val; curidx++;

  /* sink inner north/south	*/
  xoffset = nl*nr*nc;
  yoffset = hsidx*nr*nc;
  dia_val = 0;
  for(i=0, j=0; j < nc; j++){
      grididx = i*nc + j;
      R_temp = lyr[hsidx].ry / 2.0 + nc * pk->r_hs1_y; 
      cooX[curidx] = SINK_C_N+xoffset; cooY[curidx] = grididx+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
  }
  xoffset = nl*nr*nc;
  yoffset = nl*nr*nc;
  R_temp = pk->r_sp_per_y; 
  cooX[curidx] = SINK_C_N+xoffset; cooY[curidx] = SP_N+yoffset;
  cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
  R_temp = pk->r_hs2_y + pk->r_hs; 
  cooX[curidx] = SINK_C_N+xoffset; cooY[curidx] = SINK_N+yoffset;
  cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
  R_temp = pk->r_hs_c_per_y + pk->r_amb_c_per_y; dia_val += 1.0/R_temp;

  cooX[curidx] = SINK_C_N+xoffset; cooY[curidx] = SINK_C_N+yoffset;
  cooV[curidx] = dia_val; curidx++;

  xoffset = nl*nr*nc;
  yoffset = hsidx*nr*nc;
  dia_val = 0;
  for(i=nr-1, j=0; j < nc; j++){
      grididx = i*nc + j;
      R_temp = lyr[hsidx].ry / 2.0 + nc * pk->r_hs1_y; 
      cooX[curidx] = SINK_C_S+xoffset; cooY[curidx] = grididx+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
  }
  xoffset = nl*nr*nc;
  yoffset = nl*nr*nc;
  R_temp = pk->r_sp_per_y; 
  cooX[curidx] = SINK_C_S+xoffset; cooY[curidx] = SP_S+yoffset;
  cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
  R_temp = pk->r_hs2_y + pk->r_hs; 
  cooX[curidx] = SINK_C_S+xoffset; cooY[curidx] = SINK_S+yoffset;
  cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
  R_temp = pk->r_hs_c_per_y + pk->r_amb_c_per_y; dia_val += 1.0/R_temp;

  cooX[curidx] = SINK_C_S+xoffset; cooY[curidx] = SINK_C_S+yoffset;
  cooV[curidx] = dia_val; curidx++;

  /* sink inner west/east	*/
  xoffset = nl*nr*nc;
  yoffset = hsidx*nr*nc;
  dia_val = 0;
  for(i=0, j=0; i < nr; i++){
      grididx = i*nc + j;
      R_temp = lyr[hsidx].rx / 2.0 + nr * pk->r_hs1_x; 
      cooX[curidx] = SINK_C_W+xoffset; cooY[curidx] = grididx+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
  }
  xoffset = nl*nr*nc;
  yoffset = nl*nr*nc;
  R_temp = pk->r_sp_per_x; 
  cooX[curidx] = SINK_C_W+xoffset; cooY[curidx] = SP_W+yoffset;
  cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
  R_temp = pk->r_hs2_x + pk->r_hs; 
  cooX[curidx] = SINK_C_W+xoffset; cooY[curidx] = SINK_W+yoffset;
  cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
  R_temp = pk->r_hs_c_per_x + pk->r_amb_c_per_x; dia_val += 1.0/R_temp;

  cooX[curidx] = SINK_C_W+xoffset; cooY[curidx] = SINK_C_W+yoffset;
  cooV[curidx] = dia_val; curidx++;

  xoffset = nl*nr*nc;
  yoffset = hsidx*nr*nc;
  dia_val = 0;
  for(i=0, j=nc-1; i < nr; i++){
      grididx = i*nc + j;
      R_temp = lyr[hsidx].rx / 2.0 + nr * pk->r_hs1_x; 
      cooX[curidx] = SINK_C_E+xoffset; cooY[curidx] = grididx+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
  }
  xoffset = nl*nr*nc;
  yoffset = nl*nr*nc;
  R_temp = pk->r_sp_per_x; 
  cooX[curidx] = SINK_C_E+xoffset; cooY[curidx] = SP_E+yoffset;
  cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
  R_temp = pk->r_hs2_x + pk->r_hs; 
  cooX[curidx] = SINK_C_E+xoffset; cooY[curidx] = SINK_E+yoffset;
  cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
  R_temp = pk->r_hs_c_per_x + pk->r_amb_c_per_x; dia_val += 1.0/R_temp;

  cooX[curidx] = SINK_C_E+xoffset; cooY[curidx] = SINK_C_E+yoffset;
  cooV[curidx] = dia_val; curidx++;

  /* spreader north/south	*/
  xoffset = nl*nr*nc;
  yoffset = spidx*nr*nc;
  dia_val = 0;
  for(i=0, j=0; j < nc; j++){
      grididx = i*nc + j;
      R_temp = lyr[spidx].ry / 2.0 + nc * pk->r_sp1_y; 
      cooX[curidx] = SP_N+xoffset; cooY[curidx] = grididx+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
  }
  xoffset = nl*nr*nc;
  yoffset = nl*nr*nc;
  R_temp = pk->r_sp_per_y; 
  cooX[curidx] = SP_N+xoffset; cooY[curidx] = SINK_C_N+yoffset;
  cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

  cooX[curidx] = SP_N+xoffset; cooY[curidx] = SP_N+yoffset;
  cooV[curidx] = dia_val; curidx++;

  xoffset = nl*nr*nc;
  yoffset = spidx*nr*nc;
  dia_val = 0;
  for(i=nr-1, j=0; j < nc; j++){
      grididx = i*nc + j;
      R_temp = lyr[spidx].ry / 2.0 + nc * pk->r_sp1_y; 
      cooX[curidx] = SP_S+xoffset; cooY[curidx] = grididx+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
  }
  xoffset = nl*nr*nc;
  yoffset = nl*nr*nc;
  R_temp = pk->r_sp_per_y; 
  cooX[curidx] = SP_S+xoffset; cooY[curidx] = SINK_C_S+yoffset;
  cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

  cooX[curidx] = SP_S+xoffset; cooY[curidx] = SP_S+yoffset;
  cooV[curidx] = dia_val; curidx++;

  /* spreader west/east	*/
  xoffset = nl*nr*nc;
  yoffset = spidx*nr*nc;
  dia_val = 0;
  for(i=0, j=0; i < nr; i++){
      grididx = i*nc + j;
      R_temp = lyr[spidx].rx / 2.0 + nr * pk->r_sp1_x; 
      cooX[curidx] = SP_W+xoffset; cooY[curidx] = grididx+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
  }
  xoffset = nl*nr*nc;
  yoffset = nl*nr*nc;
  R_temp = pk->r_sp_per_x; 
  cooX[curidx] = SP_W+xoffset; cooY[curidx] = SINK_C_W+yoffset;
  cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

  cooX[curidx] = SP_W+xoffset; cooY[curidx] = SP_W+yoffset;
  cooV[curidx] = dia_val; curidx++;

  xoffset = nl*nr*nc;
  yoffset = spidx*nr*nc;
  dia_val = 0;
  for(i=0, j=nc-1; i < nr; i++){
      grididx = i*nc + j;
      R_temp = lyr[spidx].rx / 2.0 + nr * pk->r_sp1_x; 
      cooX[curidx] = SP_E+xoffset; cooY[curidx] = grididx+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
  }
  xoffset = nl*nr*nc;
  yoffset = nl*nr*nc;
  R_temp = pk->r_sp_per_x; 
  cooX[curidx] = SP_E+xoffset; cooY[curidx] = SINK_C_E+yoffset;
  cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

  cooX[curidx] = SP_E+xoffset; cooY[curidx] = SP_E+yoffset;
  cooV[curidx] = dia_val; curidx++;

  if(model_secondary){
      /* secondary path package nodes */
      /* PCB outer north/south	*/
      xoffset = nl*nr*nc;
      yoffset = nl*nr*nc;
      dia_val = 0;
      R_temp = pk->r_pcb2_y + pk->r_pcb; 
      cooX[curidx] = PCB_N+xoffset; cooY[curidx] = PCB_C_N+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      R_temp = pk->r_amb_sec_per; dia_val += 1.0/R_temp;
      cooX[curidx] = PCB_N+xoffset; cooY[curidx] = PCB_N+yoffset;
      cooV[curidx] = dia_val; curidx++;

      dia_val = 0;
      R_temp = pk->r_pcb2_y + pk->r_pcb; 
      cooX[curidx] = PCB_S+xoffset; cooY[curidx] = PCB_C_S+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      R_temp = pk->r_amb_sec_per; dia_val += 1.0/R_temp;
      cooX[curidx] = PCB_S+xoffset; cooY[curidx] = PCB_S+yoffset;
      cooV[curidx] = dia_val; curidx++;

      /* PCB outer west/east	*/
      dia_val = 0;
      R_temp = pk->r_pcb2_x + pk->r_pcb; 
      cooX[curidx] = PCB_W+xoffset; cooY[curidx] = PCB_C_W+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      R_temp = pk->r_amb_sec_per; dia_val += 1.0/R_temp;
      cooX[curidx] = PCB_W+xoffset; cooY[curidx] = PCB_W+yoffset;
      cooV[curidx] = dia_val; curidx++;

      dia_val = 0;
      R_temp = pk->r_pcb2_x + pk->r_pcb; 
      cooX[curidx] = PCB_E+xoffset; cooY[curidx] = PCB_C_E+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      R_temp = pk->r_amb_sec_per; dia_val += 1.0/R_temp;
      cooX[curidx] = PCB_E+xoffset; cooY[curidx] = PCB_E+yoffset;
      cooV[curidx] = dia_val; curidx++;

      /* PCB inner north/south	*/
      xoffset = nl*nr*nc;
      yoffset = pcbidx*nr*nc;
      dia_val = 0;
      for(i=0, j=0; j < nc; j++){
          grididx = i*nc + j;
          R_temp = lyr[pcbidx].ry / 2.0 + nc * pk->r_pcb1_y; 
          cooX[curidx] = PCB_C_N+xoffset; cooY[curidx] = grididx+yoffset;
          cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      }
      xoffset = nl*nr*nc;
      yoffset = nl*nr*nc;
      R_temp = pk->r_pcb_c_per_y; 
      cooX[curidx] = PCB_C_N+xoffset; cooY[curidx] = SOLDER_N+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      R_temp = pk->r_pcb2_y + pk->r_pcb; 
      cooX[curidx] = PCB_C_N+xoffset; cooY[curidx] = PCB_N+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      R_temp = pk->r_amb_sec_c_per_y; dia_val += 1.0/R_temp;

      cooX[curidx] = PCB_C_N+xoffset; cooY[curidx] = PCB_C_N+yoffset;
      cooV[curidx] = dia_val; curidx++;

      xoffset = nl*nr*nc;
      yoffset = pcbidx*nr*nc;
      dia_val = 0;
      for(i=nr-1, j=0; j < nc; j++){
          grididx = i*nc + j;
          R_temp = lyr[pcbidx].ry / 2.0 + nc * pk->r_pcb1_y; 
          cooX[curidx] = PCB_C_S+xoffset; cooY[curidx] = grididx+yoffset;
          cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      }
      xoffset = nl*nr*nc;
      yoffset = nl*nr*nc;
      R_temp = pk->r_pcb_c_per_y; 
      cooX[curidx] = PCB_C_S+xoffset; cooY[curidx] = SOLDER_S+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      R_temp = pk->r_pcb2_y + pk->r_pcb; 
      cooX[curidx] = PCB_C_S+xoffset; cooY[curidx] = PCB_S+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      R_temp = pk->r_amb_sec_c_per_y; dia_val += 1.0/R_temp;

      cooX[curidx] = PCB_C_S+xoffset; cooY[curidx] = PCB_C_S+yoffset;
      cooV[curidx] = dia_val; curidx++;

      /* PCB inner west/east	*/
      xoffset = nl*nr*nc;
      yoffset = pcbidx*nr*nc;
      dia_val = 0;
      for(i=0, j=0; i < nr; i++){
          grididx = i*nc + j;
          R_temp = lyr[pcbidx].rx / 2.0 + nr * pk->r_pcb1_x; 
          cooX[curidx] = PCB_C_W+xoffset; cooY[curidx] = grididx+yoffset;
          cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      }
      xoffset = nl*nr*nc;
      yoffset = nl*nr*nc;
      R_temp = pk->r_pcb_c_per_x; 
      cooX[curidx] = PCB_C_W+xoffset; cooY[curidx] = SOLDER_W+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      R_temp = pk->r_pcb2_x + pk->r_pcb; 
      cooX[curidx] = PCB_C_W+xoffset; cooY[curidx] = PCB_W+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      R_temp = pk->r_amb_sec_c_per_x; dia_val += 1.0/R_temp;

      cooX[curidx] = PCB_C_W+xoffset; cooY[curidx] = PCB_C_W+yoffset;
      cooV[curidx] = dia_val; curidx++;

      xoffset = nl*nr*nc;
      yoffset = pcbidx*nr*nc;
      dia_val = 0;
      for(i=0, j=nc-1; i < nr; i++){
          grididx = i*nc + j;
          R_temp = lyr[pcbidx].rx / 2.0 + nr * pk->r_pcb1_x; 
          cooX[curidx] = PCB_C_E+xoffset; cooY[curidx] = grididx+yoffset;
          cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      }
      xoffset = nl*nr*nc;
      yoffset = nl*nr*nc;
      R_temp = pk->r_pcb_c_per_x; 
      cooX[curidx] = PCB_C_E+xoffset; cooY[curidx] = SOLDER_E+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      R_temp = pk->r_pcb2_x + pk->r_pcb; 
      cooX[curidx] = PCB_C_E+xoffset; cooY[curidx] = PCB_E+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      R_temp = pk->r_amb_sec_c_per_x; dia_val += 1.0/R_temp;

      cooX[curidx] = PCB_C_E+xoffset; cooY[curidx] = PCB_C_E+yoffset;
      cooV[curidx] = dia_val; curidx++;

      /* solder north/south	*/
      xoffset = nl*nr*nc;
      yoffset = solderidx*nr*nc;
      dia_val = 0;
      for(i=0, j=0; j < nc; j++){
          grididx = i*nc + j;
          R_temp = lyr[solderidx].ry / 2.0 + nc * pk->r_solder1_y; 
          cooX[curidx] = SOLDER_N+xoffset; cooY[curidx] = grididx+yoffset;
          cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      }
      xoffset = nl*nr*nc;
      yoffset = nl*nr*nc;
      R_temp = pk->r_solder_per_y; 
      cooX[curidx] = SOLDER_N+xoffset; cooY[curidx] = SUB_N+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      R_temp = pk->r_pcb_c_per_y; 
      cooX[curidx] = SOLDER_N+xoffset; cooY[curidx] = PCB_C_N+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

      cooX[curidx] = SOLDER_N+xoffset; cooY[curidx] = SOLDER_N+yoffset;
      cooV[curidx] = dia_val; curidx++;

      xoffset = nl*nr*nc;
      yoffset = solderidx*nr*nc;
      dia_val = 0;
      for(i=nr-1, j=0; j < nc; j++){
          grididx = i*nc + j;
          R_temp = lyr[solderidx].ry / 2.0 + nc * pk->r_solder1_y; 
          cooX[curidx] = SOLDER_S+xoffset; cooY[curidx] = grididx+yoffset;
          cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      }
      xoffset = nl*nr*nc;
      yoffset = nl*nr*nc;
      R_temp = pk->r_solder_per_y; 
      cooX[curidx] = SOLDER_S+xoffset; cooY[curidx] = SUB_S+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      R_temp = pk->r_pcb_c_per_y; 
      cooX[curidx] = SOLDER_S+xoffset; cooY[curidx] = PCB_C_S+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

      cooX[curidx] = SOLDER_S+xoffset; cooY[curidx] = SOLDER_S+yoffset;
      cooV[curidx] = dia_val; curidx++;

      /* solder west/east	*/
      xoffset = nl*nr*nc;
      yoffset = solderidx*nr*nc;
      dia_val = 0;
      for(i=0, j=0; i < nr; i++){
          grididx = i*nc + j;
          R_temp = lyr[solderidx].rx / 2.0 + nr * pk->r_solder1_x; 
          cooX[curidx] = SOLDER_W+xoffset; cooY[curidx] = grididx+yoffset;
          cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      }
      xoffset = nl*nr*nc;
      yoffset = nl*nr*nc;
      R_temp = pk->r_solder_per_x; 
      cooX[curidx] = SOLDER_W+xoffset; cooY[curidx] = SUB_W+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      R_temp = pk->r_pcb_c_per_x; 
      cooX[curidx] = SOLDER_W+xoffset; cooY[curidx] = PCB_C_W+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

      cooX[curidx] = SOLDER_W+xoffset; cooY[curidx] = SOLDER_W+yoffset;
      cooV[curidx] = dia_val; curidx++;

      xoffset = nl*nr*nc;
      yoffset = solderidx*nr*nc;
      dia_val = 0;
      for(i=0, j=nc-1; i < nr; i++){
          grididx = i*nc + j;
          R_temp = lyr[solderidx].rx / 2.0 + nr * pk->r_solder1_x; 
          cooX[curidx] = SOLDER_E+xoffset; cooY[curidx] = grididx+yoffset;
          cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      }
      xoffset = nl*nr*nc;
      yoffset = nl*nr*nc;
      R_temp = pk->r_solder_per_x; 
      cooX[curidx] = SOLDER_E+xoffset; cooY[curidx] = SUB_E+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      R_temp = pk->r_pcb_c_per_x; 
      cooX[curidx] = SOLDER_E+xoffset; cooY[curidx] = PCB_C_E+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

      cooX[curidx] = SOLDER_E+xoffset; cooY[curidx] = SOLDER_E+yoffset;
      cooV[curidx] = dia_val; curidx++;

      /* substrate north/south	*/
      xoffset = nl*nr*nc;
      yoffset = subidx*nr*nc;
      dia_val = 0;
      for(i=0, j=0; j < nc; j++){
          grididx = i*nc + j;
          R_temp = lyr[subidx].ry / 2.0 + nc * pk->r_sub1_y; 
          cooX[curidx] = SUB_N+xoffset; cooY[curidx] = grididx+yoffset;
          cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      }
      xoffset = nl*nr*nc;
      yoffset = nl*nr*nc;
      R_temp = pk->r_solder_per_y; 
      cooX[curidx] = SUB_N+xoffset; cooY[curidx] = SOLDER_N+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

      cooX[curidx] = SUB_N+xoffset; cooY[curidx] = SUB_N+yoffset;
      cooV[curidx] = dia_val; curidx++;

      xoffset = nl*nr*nc;
      yoffset = subidx*nr*nc;
      dia_val = 0;
      for(i=nr-1, j=0; j < nc; j++){
          grididx = i*nc + j;
          R_temp = lyr[subidx].ry / 2.0 + nc * pk->r_sub1_y; 
          cooX[curidx] = SUB_S+xoffset; cooY[curidx] = grididx+yoffset;
          cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      }
      xoffset = nl*nr*nc;
      yoffset = nl*nr*nc;
      R_temp = pk->r_solder_per_y; 
      cooX[curidx] = SUB_S+xoffset; cooY[curidx] = SOLDER_S+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

      cooX[curidx] = SUB_S+xoffset; cooY[curidx] = SUB_S+yoffset;
      cooV[curidx] = dia_val; curidx++;

      /* substrate west/east	*/
      xoffset = nl*nr*nc;
      yoffset = subidx*nr*nc;
      dia_val = 0;
      for(i=0, j=0; i < nr; i++){
          grididx = i*nc + j;
          R_temp = lyr[subidx].rx / 2.0 + nr * pk->r_sub1_x; 
          cooX[curidx] = SUB_W+xoffset; cooY[curidx] = grididx+yoffset;
          cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      }
      xoffset = nl*nr*nc;
      yoffset = nl*nr*nc;
      R_temp = pk->r_solder_per_x; 
      cooX[curidx] = SUB_W+xoffset; cooY[curidx] = SOLDER_W+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

      cooX[curidx] = SUB_W+xoffset; cooY[curidx] = SUB_W+yoffset;
      cooV[curidx] = dia_val; curidx++;

      xoffset = nl*nr*nc;
      yoffset = subidx*nr*nc;
      dia_val = 0;
      for(i=0, j=nc-1; i < nr; i++){
          grididx = i*nc + j;
          R_temp = lyr[subidx].rx / 2.0 + nr * pk->r_sub1_x; 
          cooX[curidx] = SUB_E+xoffset; cooY[curidx] = grididx+yoffset;
          cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;
      }
      xoffset = nl*nr*nc;
      yoffset = nl*nr*nc;
      R_temp = pk->r_solder_per_x; 
      cooX[curidx] = SUB_E+xoffset; cooY[curidx] = SOLDER_E+yoffset;
      cooV[curidx] = -1.0/R_temp; curidx++; dia_val += 1.0/R_temp;

      cooX[curidx] = SUB_E+xoffset; cooY[curidx] = SUB_E+yoffset;
      cooV[curidx] = dia_val; curidx++;
  }

  if(curidx != nnz)
    fatal("Steady-state Matrix build error: less elements than nnz\n");

  coo2csc(n, nnz, cooX, cooY, cooV, asub, xa, a);

  /* Create matrix A in the format expected by SuperLU. */
  dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

  free(cooV);
  free(cooX);
  free(cooY);

  return A;
}

SuperMatrix build_steady_rhs_vector(grid_model_t *model, grid_model_vector_t *power, double **rhs)
{
  SuperMatrix B;
  int      idx, nrhs;
  int      i, j, l, m;
  int      base_idx;

  /* shortcuts	*/
  int nr = model->rows;
  int nc = model->cols;
  int nl = model->n_layers;
  int hsidx = nl - DEFAULT_PACK_LAYERS + LAYER_SINK;
  int pcbidx = LAYER_PCB; 
  int model_secondary = model->config.model_secondary;
  double cw = model->width / model->cols;
  double ch = model->height / model->rows;
  thermal_config_t *c = &model->config;
  layer_t *lyr = model->layers;
  package_RC_t *pk = &model->pack;

  if(model_secondary)
    m = (nl*nc*nr + EXTRA + EXTRA_SEC);
  else
    m = (nl*nc*nr + EXTRA);

  nrhs = 1;
  if ( !(*rhs = doubleMalloc(m*nrhs)) ) fatal("Malloc fails for rhs[].\n");

  for(l=0; l<nl; l++)
    for(i=0; i<nr; i++)
      for(j=0; j<nc; j++){
          idx = l*nr*nc + i*nc + j;
          if(l == hsidx){
              (*rhs)[idx] = c->ambient/lyr[l].rz + power->cuboid[l][i][j];
          }
          else if((l == pcbidx) && model_secondary){
              (*rhs)[idx] = power->cuboid[l][i][j] + c->ambient/(model->config.r_convec_sec *
                                                                 (model->config.s_pcb * model->config.s_pcb) / (cw * ch));
          }
          else{
              (*rhs)[idx] = power->cuboid[l][i][j];
          }
      }
  /* Package nodes */
  base_idx = nl*nc*nr;
  (*rhs)[base_idx + SP_W] = 0;
  (*rhs)[base_idx + SP_E] = 0;
  (*rhs)[base_idx + SP_N] = 0;
  (*rhs)[base_idx + SP_S] = 0;
  (*rhs)[base_idx + SINK_C_W] = c->ambient/(pk->r_hs_c_per_x + pk->r_amb_c_per_x);
  (*rhs)[base_idx + SINK_C_E] = c->ambient/(pk->r_hs_c_per_x + pk->r_amb_c_per_x);
  (*rhs)[base_idx + SINK_C_N] = c->ambient/(pk->r_hs_c_per_y + pk->r_amb_c_per_y);
  (*rhs)[base_idx + SINK_C_S] = c->ambient/(pk->r_hs_c_per_y + pk->r_amb_c_per_y);
  (*rhs)[base_idx + SINK_W] = c->ambient/(pk->r_hs_per + pk->r_amb_per);
  (*rhs)[base_idx + SINK_E] = c->ambient/(pk->r_hs_per + pk->r_amb_per);
  (*rhs)[base_idx + SINK_N] = c->ambient/(pk->r_hs_per + pk->r_amb_per);
  (*rhs)[base_idx + SINK_S] = c->ambient/(pk->r_hs_per + pk->r_amb_per);
  if(model_secondary){
      (*rhs)[base_idx + SUB_W] = 0;
      (*rhs)[base_idx + SUB_E] = 0;
      (*rhs)[base_idx + SUB_N] = 0;
      (*rhs)[base_idx + SUB_S] = 0;
      (*rhs)[base_idx + SOLDER_W] = 0;
      (*rhs)[base_idx + SOLDER_E] = 0;
      (*rhs)[base_idx + SOLDER_N] = 0;
      (*rhs)[base_idx + SOLDER_S] = 0;
      (*rhs)[base_idx + PCB_C_W] = c->ambient/(pk->r_amb_sec_c_per_x);
      (*rhs)[base_idx + PCB_C_E] = c->ambient/(pk->r_amb_sec_c_per_x);
      (*rhs)[base_idx + PCB_C_N] = c->ambient/(pk->r_amb_sec_c_per_y);
      (*rhs)[base_idx + PCB_C_S] = c->ambient/(pk->r_amb_sec_c_per_y);
      (*rhs)[base_idx + PCB_W] = c->ambient/(pk->r_amb_sec_per);
      (*rhs)[base_idx + PCB_E] = c->ambient/(pk->r_amb_sec_per);
      (*rhs)[base_idx + PCB_N] = c->ambient/(pk->r_amb_sec_per);
      (*rhs)[base_idx + PCB_S] = c->ambient/(pk->r_amb_sec_per);
  }

  dCreate_Dense_Matrix(&B, m, nrhs, *rhs, m, SLU_DN, SLU_D, SLU_GE);

  return B;
}

#if GATE_GS0 > 0

static void gs0_export_superlu(grid_model_t *model, SuperMatrix *A,
                               const double *rhs)
{
  const char *prefix = gs0_output_prefix();
  NCformat *store;
  double *values;
  int_t *row_indices;
  int_t *column_pointers;
  double *matrix;
  int count = steady_node_count(model);
  int column;

  if (!prefix)
    return;
  if (!model->config.detailed_3D_used)
    fatal("G-S0 export requires detailed-3D mode\n");
  if (count > GS0_MAX_NODE_COUNT)
    fatal("G-S0 dense export is restricted to at most 2048 nodes\n");
  if (A->Stype != SLU_NC || A->Dtype != SLU_D ||
      A->nrow != count || A->ncol != count)
    fatal("unexpected SuperLU matrix format in G-S0 export\n");

  store = (NCformat *)A->Store;
  if (!store || !store->nzval || !store->rowind || !store->colptr)
    fatal("incomplete SuperLU storage in G-S0 export\n");
  values = (double *)store->nzval;
  row_indices = (int_t *)store->rowind;
  column_pointers = (int_t *)store->colptr;
  if (store->nnz < 0 || column_pointers[0] != 0 ||
      column_pointers[count] != store->nnz)
    fatal("invalid SuperLU column pointers in G-S0 export\n");
  matrix = (double *)calloc((size_t)count * count, sizeof(double));
  if (!matrix)
    fatal("memory allocation failed in G-S0 SuperLU export\n");

  for (column = 0; column < count; column++) {
    int_t position;
    if (column_pointers[column] < 0 ||
        column_pointers[column] > column_pointers[column + 1] ||
        column_pointers[column + 1] > store->nnz)
      fatal("invalid SuperLU column range in G-S0 export\n");
    for (position = column_pointers[column];
         position < column_pointers[column + 1]; position++) {
      int_t row = row_indices[position];
      if (row < 0 || row >= count)
        fatal("invalid row index in G-S0 SuperLU export\n");
      matrix[(size_t)row * count + (size_t)column] += values[position];
    }
  }

  gs0_write_doubles(prefix, ".matrix.bin", matrix,
                    (size_t)count * count);
  gs0_write_doubles(prefix, ".rhs.bin", rhs, (size_t)count);
  gs0_write_metadata(model, prefix, "superlu", 0, 0, 0, 0, NAN);
  gs0_write_completion(prefix, ".operator.complete");
  free(matrix);
}

#endif

#if GATE_GS1_ARM == GS1_ARM_SUPERLU

static void gs1_export_superlu_operator(grid_model_t *model, SuperMatrix *A,
                                        const double *rhs)
{
  const char *prefix = getenv("HOTSPOT_GS1_OPERATOR_PREFIX");
  unsigned short endian_probe = 1;
  NCformat *store;
  double *values;
  int_t *row_indices;
  int_t *column_pointers;
  int64_t *rows64;
  int64_t *columns64;
  int count = steady_node_count(model);
  int failed;
  int_t position;
  FILE *stream;
  if (!prefix || !prefix[0])
    return;
  if (A->Stype != SLU_NC || A->Dtype != SLU_D ||
      A->nrow != count || A->ncol != count)
    fatal("unexpected SuperLU matrix format in G-S1 export\n");
  store = (NCformat *)A->Store;
  if (!store || !store->nzval || !store->rowind || !store->colptr ||
      store->nnz < 0)
    fatal("incomplete SuperLU storage in G-S1 export\n");
  values = (double *)store->nzval;
  row_indices = (int_t *)store->rowind;
  column_pointers = (int_t *)store->colptr;
  if (column_pointers[0] != 0 || column_pointers[count] != store->nnz)
    fatal("invalid SuperLU column pointers in G-S1 export\n");
  rows64 = (int64_t *)calloc((size_t)store->nnz, sizeof(int64_t));
  columns64 = (int64_t *)calloc((size_t)count + 1, sizeof(int64_t));
  if (!rows64 || !columns64)
    fatal("memory allocation failed in G-S1 operator export\n");
  for (position = 0; position < store->nnz; position++)
    rows64[position] = (int64_t)row_indices[position];
  for (position = 0; position <= count; position++)
    columns64[position] = (int64_t)column_pointers[position];
  gs1_write_bytes(prefix, ".values.bin", values, sizeof(double),
                  (size_t)store->nnz);
  gs1_write_bytes(prefix, ".rowind.i64.bin", rows64, sizeof(int64_t),
                  (size_t)store->nnz);
  gs1_write_bytes(prefix, ".colptr.i64.bin", columns64, sizeof(int64_t),
                  (size_t)count + 1);
  gs1_write_bytes(prefix, ".rhs.bin", rhs, sizeof(double), (size_t)count);
  gs1_write_node_map(model, prefix);
  stream = gs1_open_output(prefix, ".meta", "w");
  fprintf(stream, "format=hotspot-gs1-operator-v1\n");
  fprintf(stream, "storage=superlu-nc-csc\n");
  fprintf(stream, "scalar=float64-native\n");
  fprintf(stream, "index_storage=int64-native\n");
  fprintf(stream, "byte_order=%s\n",
          *((unsigned char *)&endian_probe) ? "little" : "big");
  fprintf(stream, "source_index_bytes=%lu\n", (unsigned long)sizeof(int_t));
  fprintf(stream, "node_count=%d\n", count);
  fprintf(stream, "nnz=%lld\n", (long long)store->nnz);
  fprintf(stream, "rows=%d\n", model->rows);
  fprintf(stream, "cols=%d\n", model->cols);
  fprintf(stream, "layers=%d\n", model->n_layers);
  failed = ferror(stream);
  if (fclose(stream) != 0)
    failed = 1;
  if (failed)
    fatal("unable to write complete G-S1 operator metadata\n");
  gs1_write_completion(prefix, ".operator.complete");
  free(rows64);
  free(columns64);
}

#endif

void direct_SLU(grid_model_t *model, grid_model_vector_t *power, grid_model_vector_t *temp)
{
  SuperMatrix A, L, U, B;
  double   *rhs;
  int      *perm_r; /* row permutations from partial pivoting */
  int      *perm_c; /* column permutation vector */
  int      info;
  superlu_options_t options;
  SuperLUStat_t stat;

  int          i, dim;
  DNformat     *Astore;
  double       *dp;
#if GATE_GS1_ARM == GS1_ARM_SUPERLU
  double gs1_started_s;
#endif

  dim = steady_node_count(model);
#if GATE_GS0 > 0
  if (dim > GS0_MAX_NODE_COUNT)
    fatal("G-S0 dense export is restricted to at most 2048 nodes\n");
#endif

#if GATE_GS1_ARM == GS1_ARM_SUPERLU
  gs1_started_s = gs1_now();
#endif
  A = build_steady_grid_matrix(model);
  B = build_steady_rhs_vector(model, power, &rhs);
#if GATE_GS1_ARM == GS1_ARM_SUPERLU
  model->gs1_assembly_s += gs1_now() - gs1_started_s;
  {
    NCformat *store = (NCformat *)A.Store;
    model->gs1_operator_bytes =
      (size_t)store->nnz * (sizeof(double) + sizeof(int_t)) +
      ((size_t)dim + 1) * sizeof(int_t) + (size_t)dim * sizeof(double);
  }
  gs1_export_superlu_operator(model, &A, rhs);
#endif

#if GATE_GS0 > 0
  gs0_export_superlu(model, &A, rhs);
#endif

#if GATE_GS1_ARM == GS1_ARM_SUPERLU
  /* dgssv does not expose consistently comparable wall-clock subphases.
   * Charge every post-assembly setup and solve operation to one combined,
   * monotonic-clock direct-factor/solve phase without inventing a split. */
  gs1_started_s = gs1_now();
  model->gs1_direct_permutation_bytes = (size_t)dim * sizeof(int_t) * 2;
#endif
  if ( !(perm_r = intMalloc(dim)) ) fatal("Malloc fails for perm_r[].\n");
  if ( !(perm_c = intMalloc(dim)) ) fatal("Malloc fails for perm_c[].\n");

  /* Set the default input options. */
  set_default_options(&options);
  options.ColPerm = MMD_AT_PLUS_A;
  options.DiagPivotThresh = 0.01;
  options.SymmetricMode = YES;
  options.Equil = YES;

  /* Initialize the statistics variables. */
  StatInit(&stat);

  /* Solve the linear system. */
  dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
  if (info != 0)
    fatal("SuperLU failed to solve the detailed-3D steady system\n");
#if GATE_GS1_ARM == GS1_ARM_SUPERLU
  {
    SCformat *l_store = (SCformat *)L.Store;
    NCformat *u_store = (NCformat *)U.Store;
    if (!l_store || !u_store)
      fatal("SuperLU returned incomplete factors\n");
    model->gs1_direct_l_nnz = (long long)l_store->nnz;
    model->gs1_direct_u_nnz = (long long)u_store->nnz;
  }
#endif

  Astore = (DNformat *) B.Store;
  dp = (double *) Astore->nzval;
  /* Copy the solution into the caller-provided destination vector. */
  for(i=0; i<dim; ++i){
      temp->cuboid[0][0][i] = dp[i];
  }
#if GATE_GS0 > 0
  gs0_write_state(model, temp->cuboid[0][0]);
#endif

  SUPERLU_FREE (rhs);
  SUPERLU_FREE (perm_r);
  SUPERLU_FREE (perm_c);
  Destroy_CompCol_Matrix(&A);
  Destroy_SuperMatrix_Store(&B);
  Destroy_SuperNode_Matrix(&L);
  Destroy_CompCol_Matrix(&U);
  StatFree(&stat);
#if GATE_GS1_ARM == GS1_ARM_SUPERLU
  model->gs1_factor_s += gs1_now() - gs1_started_s;
#endif
}
#endif
