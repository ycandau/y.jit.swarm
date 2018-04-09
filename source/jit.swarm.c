/**
*  jit.swarm - A Max MSP object for swarm sonification
*  Yves Candau - ycandau@gmail.com
*/

/****************************************************************
*  @todo
*
*  Optimization:
*    - vectorialize position and velocity instead of per agent
*    - implement a nearest neighbor search instead of brute force
*    - check cost of randomization function, rand() vs jit_rand() vs by hand
*/

/****************************************************************
*  Header files
*/
#include "jit.swarm.h"

/****************************************************************
*  Jitter class pointer.
*  Declared as extern in header file.
*/
void *_jit_swarm_class = NULL;

/****************************************************************
*  Jitter class initialization.
*
*  Create the Jitter class and the MOP wrappers.
*  Create and setup the class attributes.
*  Create and setup the class methods.
*
*  @return A Jitter error type: JIT_ERR_NONE
*/
t_jit_err jit_swarm_init(void)
{
  TRACE("jit_swarm_init");

	void *mop;

  // Create the Jitter class
	_jit_swarm_class = jit_class_new(
    "jit_swarm",              // Jitter class name
    (method)jit_swarm_new,    // Jitter class constructor
    (method)jit_swarm_free,   // Jitter class destructor
		sizeof(t_swarm),          // Jitter class size
    0L);

	// Create a MOP wrapper
	mop = jit_object_new(_jit_sym_jit_mop, 0, 4);

  // Unlink everything
  jit_mop_output_nolink(mop, 1);
  jit_mop_output_nolink(mop, 2);
  jit_mop_output_nolink(mop, 3);
  jit_mop_output_nolink(mop, 4);

  // Setup individual jit_mop_io objects for each matrix output
  // Position matrix: 3 float N
  void *output1 = jit_object_method(mop, _jit_sym_getoutput, 1);
  jit_attr_setlong(output1, _jit_sym_minplanecount, 3);
  jit_attr_setlong(output1, _jit_sym_maxplanecount, 3);
  jit_attr_setlong(output1, _jit_sym_mindim, 0);
  jit_attr_setlong(output1, _jit_sym_maxdim, 2*AGENT_MAX_DEF);
  jit_attr_setlong(output1, _jit_sym_mindimcount, 1);
  jit_attr_setlong(output1, _jit_sym_maxdimcount, 1);
  jit_attr_setsym(output1, _jit_sym_types, SWR_FLOAT_SYM);

  // Signed normalized position matrix [-1, 1] : 3 float N
  void *output2 = jit_object_method(mop, _jit_sym_getoutput, 2);
  jit_attr_setlong(output2, _jit_sym_minplanecount, 3);
  jit_attr_setlong(output2, _jit_sym_maxplanecount, 3);
  jit_attr_setlong(output2, _jit_sym_mindim, 0);
  jit_attr_setlong(output2, _jit_sym_maxdim, 2*AGENT_MAX_DEF);
  jit_attr_setlong(output2, _jit_sym_mindimcount, 1);
  jit_attr_setlong(output2, _jit_sym_maxdimcount, 1);
  jit_attr_setsym(output2, _jit_sym_types, SWR_FLOAT_SYM);

  // Scale matrix: 3 float N
  void *output3 = jit_object_method(mop, _jit_sym_getoutput, 3);
  jit_attr_setlong(output3, _jit_sym_minplanecount, 3);
  jit_attr_setlong(output3, _jit_sym_maxplanecount, 3);
  jit_attr_setlong(output3, _jit_sym_mindim, 0);
  jit_attr_setlong(output3, _jit_sym_maxdim, 2*AGENT_MAX_DEF);
  jit_attr_setlong(output3, _jit_sym_mindimcount, 1);
  jit_attr_setlong(output3, _jit_sym_maxdimcount, 1);
  jit_attr_setsym(output3, _jit_sym_types, SWR_FLOAT_SYM);

  // Colour matrix: 4 float N (RGBA)
  void *output4 = jit_object_method(mop, _jit_sym_getoutput, 4);
  jit_attr_setlong(output4, _jit_sym_minplanecount, 4);
  jit_attr_setlong(output4, _jit_sym_maxplanecount, 4);
  jit_attr_setlong(output4, _jit_sym_mindim, 0);
  jit_attr_setlong(output4, _jit_sym_maxdim, 2*AGENT_MAX_DEF);
  jit_attr_setlong(output4, _jit_sym_mindimcount, 1);
  jit_attr_setlong(output4, _jit_sym_maxdimcount, 1);
  jit_attr_setsym(output4, _jit_sym_types, SWR_FLOAT_SYM);

  // Add the MOP to the Jitter class
	jit_class_addadornment(_jit_swarm_class, (t_jit_object *)mop);

  // Initialize the object methods
  swarm_init_method(_jit_swarm_class);

  // Initialize the object attributes
  swarm_init_attr(_jit_swarm_class);

  // Register the Jitter class
	jit_class_register(_jit_swarm_class);

	return JIT_ERR_NONE;
}

/****************************************************************
*  Jitter class constructor.
*  
*  @param x A pointer to the t_swarm structure.
*/
t_swarm *jit_swarm_new(void)
{
  TRACE("jit_swarm_new");

  t_swarm *x = NULL;
  x = (t_swarm *)jit_object_alloc(_jit_swarm_class);

  // Test the allocation and return NULL on error
  if (!x) {
    jit_object_error(NULL, "jit.swarm:  Allocation error.");
    return NULL;
  }

  // Initialize the random function
  srand((unsigned int)time(NULL));

  // Initialize the evolutionary algorithm
  // @NB: to do before initializing the agents
  ea_init_gene_repr(x);
  ea_init_gene_rec(x);

  // Initialize the parameters
  swarm_init_param(x);   

  // Allocate and initialize the agents
  x->agent_max = 0;
  x->agent_arr = NULL;
  x->agent_act = NULL;
  swarm_new_agents(x, 0, AGENT_MAX_DEF);
  
  return x;
}

/****************************************************************
*  Jitter class destructor.
*  
*  @param x A pointer to the t_swarm structure.
*/
void jit_swarm_free(t_swarm *x)
{
  TRACE("jit_swarm_free");

  if (x->agent_arr) { jit_disposeptr((char *)x->agent_arr); }
  if (x->agent_act) { jit_disposeptr((char *)x->agent_act); }

  // The Jitter object is freed from the Max wrapper

  return;
}

/****************************************************************
*  Jitter matrices calculation.
*
*  Ensure that the matrices are valid, compatible, and if so, process.
*  Called by matrix_calc and bang messages to the object.
*
*  @param x A pointer to the t_swarm structure.
*  @param inputs A pointer to the list of input matrices.
*  @param outputs A pointer to the list of output matrices.
*
*  @return A Jitter error type:  JIT_ERR_NONE or
*    JIT_ERR_INVALID_PTR: Invalid pointers for the object or the matrices
*    JIT_ERR_INVALID_OUTPUT: Invalid pointers for the data arrays
*/
t_jit_err jit_swarm_matrix_calc(t_swarm *x, void *inputs, void *outputs)
{
  TRACE_F("jit_swarm_matrix_calc");

  t_jit_err err = JIT_ERR_NONE;
  void *out_matr_1, *out_matr_2, *out_matr_3, *out_matr_4;
  t_ptr_int out_lock_1, out_lock_2, out_lock_3, out_lock_4;   // using t_ptr_int instead of long to avoid warning
  t_jit_matrix_info out_minfo_1, out_minfo_2, out_minfo_3, out_minfo_4;
  t_swr_float *out_data_1, *out_data_2, *out_data_3, *out_data_4;

  // Get the matrices
  out_matr_1 = jit_object_method(outputs, _jit_sym_getindex, 0);
  out_matr_2 = jit_object_method(outputs, _jit_sym_getindex, 1);
  out_matr_3 = jit_object_method(outputs, _jit_sym_getindex, 2);
  out_matr_4 = jit_object_method(outputs, _jit_sym_getindex, 3);

  // If the object and all the matrices are valid, then process
  if (x && out_matr_1 && out_matr_2 && out_matr_3 && out_matr_4) {

    // Lock the matrices
    out_lock_1 = (t_ptr_int)jit_object_method(out_matr_1, _jit_sym_lock, 1);
    out_lock_2 = (t_ptr_int)jit_object_method(out_matr_2, _jit_sym_lock, 1);
    out_lock_3 = (t_ptr_int)jit_object_method(out_matr_3, _jit_sym_lock, 1);
    out_lock_4 = (t_ptr_int)jit_object_method(out_matr_4, _jit_sym_lock, 1);

    // Get the matrices info and
    jit_object_method(out_matr_1, _jit_sym_getinfo, &out_minfo_1);
    jit_object_method(out_matr_2, _jit_sym_getinfo, &out_minfo_2);
    jit_object_method(out_matr_3, _jit_sym_getinfo, &out_minfo_3);
    jit_object_method(out_matr_4, _jit_sym_getinfo, &out_minfo_4);

    out_minfo_1.dim[0] = x->agent_cnt;
    out_minfo_2.dim[0] = x->agent_cnt;
    out_minfo_3.dim[0] = x->agent_cnt;
    out_minfo_4.dim[0] = x->agent_cnt;

    jit_object_method(out_matr_1, _jit_sym_setinfo, &out_minfo_1);
    jit_object_method(out_matr_2, _jit_sym_setinfo, &out_minfo_2);
    jit_object_method(out_matr_3, _jit_sym_setinfo, &out_minfo_3);
    jit_object_method(out_matr_4, _jit_sym_setinfo, &out_minfo_4);

    // Get the matrices data pointers and test
    jit_object_method(out_matr_1, _jit_sym_getdata, &out_data_1);
    jit_object_method(out_matr_2, _jit_sym_getdata, &out_data_2);
    jit_object_method(out_matr_3, _jit_sym_getdata, &out_data_3);
    jit_object_method(out_matr_4, _jit_sym_getdata, &out_data_4);

    if (!out_data_1 || !out_data_2 || !out_data_3 || !out_data_4) {
      err = JIT_ERR_INVALID_OUTPUT;
      goto out;
    }

    // Copy the data to the Jitter matrices
    t_agent *agent;
    t_swr_float *pos, *col, size;
    t_swr_float *out_data_1_ref = out_data_1;

    for (t_swr_ind ind = 0; ind < x->agent_cnt; ind++) {
      agent = x->agent_arr + (x->agent_act[ind]);
      pos = agent->pos;
      col = agent->col;
      size = agent->size * SIZE_MUL_BASE;

      // Copy the position matrix
      *out_data_1++ = *pos++;
      *out_data_1++ = *pos++;
      *out_data_1++ = *pos++;

      // Copy the size matrix
      *out_data_3++ = size;
      *out_data_3++ = size;
      *out_data_3++ = size;

      // Copy the colour matrix
      *out_data_4++ = *col++;
      *out_data_4++ = *col++;
      *out_data_4++ = *col++;
      *out_data_4++ = *col++;
    }

    // Copy the signed normalized position matrix
    out_data_1 = out_data_1_ref;
    for (t_swr_ind ind = 0; ind < x->agent_cnt; ind++) {
      *out_data_2++ = (2.0f * *out_data_1++ - x->length_x) / x->length_x;
      *out_data_2++ = (2.0f * *out_data_1++ - x->length_y) / x->length_x;
      *out_data_2++ = (2.0f * *out_data_1++ - x->length_z) / x->length_x;
    }

    // ... else return an error
  } else {
    return JIT_ERR_INVALID_PTR;
  }

out:
  // Unlock the matrices and return
  jit_object_method(out_matr_1, _jit_sym_lock, out_lock_1);
  jit_object_method(out_matr_2, _jit_sym_lock, out_lock_2);
  jit_object_method(out_matr_3, _jit_sym_lock, out_lock_3);
  jit_object_method(out_matr_4, _jit_sym_lock, out_lock_4);

  return err;
}

/****************************************************************
*  Recursive matrix calculation.
*
*  A recursive function to handle higher dimension matrices, by processing 2D sections.
*  Not used here. Included to reference the general Jitter object template.
*  Could be called directly or through a parallelization function.
*  See jit_parallel_ndim_simplecalc1() for instance.
*
*  @param x A pointer to the t_swarm structure.
*  @param dimcount The number of dimensions, recursively decreased.
*  @param dim A pointer to the list of sizes per dimension list.
*  @param planecount The planecount of the matrix.
*  @param in_minfo A pointer to the input matrix information structure.
*  @param bip A base pointer to the input matrix data array.
*  @param out_minfo A pointer to the output matrix information structure.
*  @param bop A base pointer to the output matrix data array.
*/
void jit_swarm_calculate_ndim(t_swarm *x, long dimcount, long *dim, long planecount,
  t_jit_matrix_info *in_minfo, char *bip, t_jit_matrix_info *out_minfo, char *bop)
{
  return;
}

/****************************************************************
*   Post information on a Jitter matrix information structure.
*
*   @param minfo A pointer to a matrix information structure.
*
*   @todo Currently only posts information for up to 3 planes. Generalize to any number.
*/
void jit_minfo_post(t_jit_matrix_info *minfo)
{
  TRACE("jit_minfo_post");

  if (minfo->dimcount < 1) {
    jit_object_error(NULL, "matrix info:  Invalid dimcount:  %i", minfo->dimcount);
    return;
  }

  switch (minfo->dimcount) {
  case 1:
    post("matrix info:  planecount: %i - type: %s - dimcount: %i - dim: %i - dimstride: %i - size : %i",
      minfo->planecount, minfo->type->s_name, minfo->dimcount, minfo->dim[0], minfo->dimstride[0], minfo->size);
    break;

  case 2:
    post("matrix info:  planecount: %i - type: %s - dimcount: %i - dim: %i %i - dimstride: %i %i - size : %i",
      minfo->planecount, minfo->type->s_name, minfo->dimcount, minfo->dim[0], minfo->dim[1], minfo->dimstride[0],
      minfo->dimstride[1], minfo->size);
    break;

  default:
    post("matrix info:  planecount: %i - type: %s - dimcount: %i - dim: %i %i %i - dimstride: %i %i %i - size : %i",
      minfo->planecount, minfo->type->s_name, minfo->dimcount, minfo->dim[0], minfo->dim[1], minfo->dim[2],
      minfo->dimstride[0], minfo->dimstride[1], minfo->dimstride[2], minfo->size);
    break;
  }

  return;
}

/****************************************************************
*  Initialize the swarm system parameters to defaults.
*  
*  @param x A pointer to the t_swarm structure.
*/
void swarm_init_param(t_swarm *x)
{
  TRACE("swarm_init_param");

  // Attributes:  Basic - Swarm parameters

  jit_attr_setfloat(x, gensym("dt"),       1.0f);
  jit_attr_setfloat(x, gensym("vel_max"), 12.0f);
  jit_attr_setfloat(x, gensym("size"),     1.0f);

  // Attributes:  Evolutionary parameters
  
  jit_attr_setfloat(x, gensym("fit_xy_mul"), 1.0f);
  jit_attr_setfloat(x, gensym("fit_vel_mul"), 1.0f);
  jit_attr_setfloat(x, gensym("mut_d"), 0.1f);
  jit_attr_setfloat(x, gensym("mut_p"), 0.3f);
  jit_attr_setfloat(x, gensym("select_fi"), 1);
  jit_attr_setfloat(x, gensym("cross_fi"), 1);
  jit_attr_setfloat(x, gensym("rank_cnt"), 20);
  jit_attr_setfloat(x, gensym("init_d"), 0.5f);

  // Attributes:  Basic - Timing

  jit_attr_setfloat(x, gensym("dorm_t"),   1.0f);
  jit_attr_setfloat(x, gensym("dorm_dt"),  0.1f);
  jit_attr_setfloat(x, gensym("hatch_t"),  1.0f);
  jit_attr_setfloat(x, gensym("hatch_dt"), 0.1f);
  jit_attr_setfloat(x, gensym("living_t"),  10.0f);
  jit_attr_setfloat(x, gensym("living_dt"),  0.1f);
  jit_attr_setfloat(x, gensym("dying_t"),  1.0f);
  jit_attr_setfloat(x, gensym("dying_dt"), 0.1f);

  // Attributes:  Basic - Agent modifiers

  jit_attr_setfloat(x, gensym("inert_mul"),  1.0f);
  jit_attr_setfloat(x, gensym("cohes_mul"),  1.0f);
  jit_attr_setfloat(x, gensym("align_mul"),  1.0f);
  jit_attr_setfloat(x, gensym("attra_mul"),  1.0f);
  jit_attr_setfloat(x, gensym("attra_thr"),  1.0f);
  jit_attr_setfloat(x, gensym("separ_mul"),  1.0f);
  jit_attr_setfloat(x, gensym("separ_thr"),  1.0f);
  jit_attr_setfloat(x, gensym("bound_mul"),  1.0f);
  jit_attr_setfloat(x, gensym("bound_thr"),  1.0f);
  jit_attr_setfloat(x, gensym("brown_mul"),  1.0f);
  jit_attr_setfloat(x, gensym("brown_prob"), 1.0f);

  // Attributes:  Basic - Display

  jit_attr_setlong(x, gensym("length_x"), 1280);
  jit_attr_setlong(x, gensym("length_y"), 720);
  jit_attr_setlong(x, gensym("length_z"), 500);

  // Dumpout all attribute values
  jit_object_method(x, gensym("getstate"));   // not working currently

  return;
}

/****************************************************************
*  Allocate and initialize all the agents in the array.
*
*  @param x A pointer to the t_swarm structure.
*  @param max The length of the array to allocate.
*
*  @return A Jitter error type:  JIT_ERR_NONE or
*    JIT_ERR_OUT_OF_MEM: if there is an allocation error
*/
t_jit_err swarm_new_agents(t_swarm *x, t_swr_ind cur, t_swr_ind max)
{
  TRACE("swarm_new_agents");

  // If the number is unchanged and the arrays already allocated: initialize
  if ((max == x->agent_max) && (x->agent_arr) && (x->agent_act)) { return JIT_ERR_NONE; }

  // If the arrays are already allocated: free and reset to NULL
  if (x->agent_arr || x->agent_act) {
    jit_disposeptr((char *)x->agent_arr);
    x->agent_arr = NULL;
    jit_disposeptr((char *)x->agent_act);
    x->agent_act = NULL;
  }

  // Allocate and test the allocation
  x->agent_max = max;
  x->agent_arr = (t_agent *)jit_newptr(sizeof(t_agent) * x->agent_max);
  x->agent_act = (t_swr_ind *)jit_newptr(sizeof(t_swr_ind) * x->agent_max);

  if ((!x->agent_arr) || (!x->agent_act)) {
    x->agent_cnt = 0;
    x->agent_max = 0;
    jit_object_error((t_object *)x, "swarm_new_agents:  Allocation error for the agent arrays.");
    return JIT_ERR_OUT_OF_MEM;
  }

  // Set all agents inactive first
  x->agent_cnt = 0;
  for (t_swr_ind ind = 0; ind < x->agent_max; ind++) {
    agent_init(x, x->agent_arr + ind);
    x->agent_act[ind] = 0;
  }

  // Add the required number of active agents
  for (t_swr_ind cnt = 0; cnt < cur; cnt++) { swarm_add(x, NULL, 0, NULL); }
  
  return JIT_ERR_NONE;
}

/****************************************************************
*  Iterate the swarm system.
*
*  Calculate one iteration of the swarm system, and update the Jitter matrices.
*
*  @param x A pointer to the t_swarm structure.
*/
void swarm_iter(t_swarm *x)
{
  TRACE_F("jit_swarm_iter");

  t_agent *agent = NULL;

  // Calculate the average position and velocity
  swarm_average(x);

  // Loop through all the active agents
  for (t_swr_ind ind = 0; ind < x->agent_cnt; ind++) {
    agent = x->agent_arr + (x->agent_act[ind]);

    // Change the state if necessary
    if (!agent->cntd) { agent_state_change(x, agent); }
    if (IS_INACTIVE(agent)) { continue; }

    switch (agent->state & ~ST_ACTIVE) {

    case ST_DORMANT:  break;
    case ST_HATCHING: agent_iter_hatching(x, agent); break;
    case ST_FLYING:   agent_iter_flying(x, agent); break;
    case ST_DYING:    agent_iter_dying(x, agent); break;
    default: break;
    }
    agent->age++;
    agent->cntd--;
  }

  // Calculate the new positions and track statistics
  swarm_move(x);
  swarm_track_stat(x);

  return;
}

/****************************************************************
*  Calculate the average position and velocity of the swarm.
*/
void swarm_average(t_swarm *x)
{
  t_agent     *agent = NULL;
  t_swr_float *X = NULL;
  t_swr_float *dX = NULL;
  t_swr_float *Xa = x->pos_avg;
  t_swr_float *dXa = x->vel_avg;
  t_swr_float  mult;
  t_swr_ind    cnt = 0;

  // Initialize the vectors to 0
  V2D_V_e_S(Xa, 0.0f);
  V2D_V_e_S(dXa, 0.0f);

  // Loop through the active agents and add the vectors
  for (t_swr_ind ind = 0; ind < x->agent_cnt; ind++) {
    agent = x->agent_arr + (x->agent_act[ind]);
    if (agent->state & ST_DORMANT) { continue; }
    X = agent->pos;
    dX = agent->vel;
    V2D_V_pe_V(Xa, X);
    V2D_V_pe_V(dXa, dX);
    cnt++;
  }

  // Divide by the number of active agents
  mult = 1.0f / cnt;
  V2D_V_te_S(Xa, mult);
  V2D_V_te_S(dXa, mult);

  return;
}

/****************************************************************
*  Calculate the new position vectors of the swarm's agents.
*/
void swarm_move(t_swarm *x)
{
  t_agent     *agent = NULL;
  t_swr_float *X = NULL;
  t_swr_float *dX = NULL;

  // Calculate new positions
  for (t_swr_ind ind = 0; ind < x->agent_cnt; ind++) {
    agent = x->agent_arr + (x->agent_act[ind]);
    X = agent->pos;
    dX = agent->vel;
    V2D_V_pe_S_V(X, x->dt, dX);
    X[2] = x->length_z / 2.0f; dX[2] = 0.0f;
  }

  return;
}

/****************************************************************
*  Track statistics on the agents.
*
*  Then used to evaluate the fitness of the agents.
*
*  @param x A pointer to the t_swarm structure.
*/
void swarm_track_stat(t_swarm *x)
{
  TRACE_F("swarm_track_stat");

  // Loop through the active agents
  t_agent *agent = NULL;
  t_swr_float i_DX, i_DY, i_DV;
  int X, Y, V, bin_ind;

  i_DX = (t_swr_float)BIN_CNT_X / x->length_x;
  i_DY = (t_swr_float)BIN_CNT_Y / x->length_y;
  i_DV = (t_swr_float)BIN_CNT_VEL / x->vel_max;

  // Loop through the active and flying agents
  for (t_swr_ind ind = 0; ind < x->agent_cnt; ind++) {
    agent = x->agent_arr + (x->agent_act[ind]);
    if (!(agent->state & ST_FLYING)) { continue; }

    // Bin the agent's position
    X = (int)(agent->pos[0] * i_DX);
    Y = (int)(agent->pos[1] * i_DY);

    // Clamp the index, use the last index for out of bounds
    bin_ind = ((X < 0) || (Y < 0) || (X >= BIN_CNT_X) || (Y >= BIN_CNT_Y)) ?
      BIN_CNT_X * BIN_CNT_Y :
      X + Y * BIN_CNT_X;

    // Increment the corresponding bin counter
    (agent->bin_xy[bin_ind])++;

    // Bin the agent's velocity
    V = (int)(agent->nvel * i_DV);
    
    // Clamp the index
    bin_ind = (V >= BIN_CNT_VEL) ?
      BIN_CNT_VEL - 1 :
      (V < 0) ? 0 : V;
    
    // Increment the corresponding bin counter
    (agent->bin_vel[bin_ind])++;
  }

  return;
}

/****************************************************************
*  Refresh the array of active agents.
*
*  Called every time there is a change in active status.
*
*  @param x A pointer to the t_swarm structure.
*/
void swarm_refresh_act(t_swarm *x)
{
  t_agent *agent;
  t_swr_ind *agent_act = x->agent_act;

  for (t_swr_ind ind = 0; ind < x->agent_max; ind++) {
    agent = x->agent_arr + ind;
    if (IS_ACTIVE(agent)) { *agent_act++ = ind; }
  }

  if (agent_act - x->agent_act != x->agent_cnt) {
    jit_object_error((t_object *)x, "swarm_refresh_act:  Invalid active agents array.");
  }

  return;
}

/****************************************************************
*  Add an agent to the active list.
*/
t_agent *agent_add(t_swarm *x)
{
  TRACE("agent_add");

  if (x->agent_cnt >= x->agent_max) { return NULL; }

  // Loop to find the first inactive agent
  t_agent *agent = x->agent_arr;
  while ((agent != x->agent_arr + x->agent_max) && (IS_ACTIVE(agent))) { agent++; }

  if (agent == x->agent_arr + x->agent_max) {
    jit_object_error((t_object *)x, "agent_add:  No inactive agent found.");
    return NULL;
  }

  // Set the agent to active
  agent_init(x, agent);
  agent->state = ST_ACTIVE | ST_DORMANT;
  x->agent_cnt++;

  // Refresh the array
  swarm_refresh_act(x);

  return agent;
}

/****************************************************************
*  Remove an agent from the active list.
*
*   @param agent Either NULL to remove the first active agent from the back
*     Or the pointer of an active agent
*/
void agent_remove(t_swarm *x, t_agent *agent)
{
  TRACE("agent_remove");

  // If a pointer to an agent is provided
  if (agent) {
    if (IS_INACTIVE(agent)) {
      jit_object_error((t_object *)x, "agent_remove:  Invalid action. Agent is inactive.");
      return;
    }
  }

  // Otherwise look for an active agent to remove
  else {
    if (x->agent_cnt < 1) { return; }

    // Loop from the end to find the first active agent
    agent = x->agent_arr + x->agent_max;
    do { agent--; } while (IS_INACTIVE(agent) && (agent != x->agent_arr));

    if (IS_INACTIVE(agent)) {
      jit_object_error((t_object *)x, "swarm_remove:  No active agent found.");
      return;
    }
  }

  // Set the agent to inactive
  agent->state = ST_NULL;
  x->agent_cnt--;

  // Refresh the array
  swarm_refresh_act(x);

  return;
}

/****************************************************************
*  Initialize an agent.
*
*  Reset an agent by setting its parameters to default values.
*  Does not allocate the structure or change its active status.
*
*  @param x A pointer to the t_swarm structure.
*  @param agent A pointer to the t_agent structure.
*/
void agent_init(t_swarm *x, t_agent *agent)
{
  TRACE_F("agent_init");

  agent->state = ST_NULL;
  agent->age = 0;
  agent->cntd = 0;

  agent->pos[0] = x->length_x * (t_swr_float)rand() / RAND_MAX;
  agent->pos[1] = x->length_y * (t_swr_float)rand() / RAND_MAX;
  agent->pos[2] = 0.0f;
  agent->vel[0] = 0.0f;
  agent->vel[1] = 0.0f;
  agent->vel[2] = 0.0f;
  agent->nvel = V2D_N(agent->vel);
  agent->size_ref = x->size;
  agent->size = agent->size_ref;

  // Use three out of phase sines for color
  t_swr_ind ind = (t_swr_ind)(agent - x->agent_arr);
  agent->col_ref[0] = 0.5f * (sinf(0.4f * ind + 0 * PI / 3) + 1);
  agent->col_ref[1] = 0.5f * (sinf(0.4f * ind + 2 * PI / 3) + 1);
  agent->col_ref[2] = 0.5f * (sinf(0.4f * ind + 4 * PI / 3) + 1);
  agent->col_ref[3] = 1.0f;
  for (int c = 0; c < 3; c++) { agent->col[c] = 0.0; }
  agent->col[4] = 1.0;

  // Initialize the genes
  ea_init_gene_nval(agent->gene_nval, x);
  ea_calc_gene_val(agent->gene_nval, agent->gene_val, x);

  // Initialize the bin counters
  for (int g = 0; g < BIN_CNT_X * BIN_CNT_Y + 1; g++) { agent->bin_xy[g] = 0; }
  for (int g = 0; g < BIN_CNT_VEL; g++) { agent->bin_vel[g] = 0; }

  return;
}

/****************************************************************
*  Set an agent to hatchin cntd iterations.
*  
*  @param cntd countdown to hatching,
*    or if CNTD_RAND, the countdown is set randomly
*/
void agent_seed(t_swarm *x, t_agent *agent, t_swr_cnt cntd)
{
  TRACE("agent_seed");

  if (IS_INACTIVE(agent)) {
    jit_object_error((t_object *)x, "agent_seed:  Invalid action. Agent is inactive.");
    return;
  }

  agent_init(x, agent);
  agent->state = ST_ACTIVE | ST_DORMANT;

  if (cntd != CNTD_RAND) { agent->cntd = cntd; }
  else { agent->cntd = rand_time(x->dorm_t, x->dorm_dt); }

  return;
}

/****************************************************************
*  Post information on an agent to the Max console.
*
*  @param x A pointer to the t_swarm structure.
*  @param agent A pointer to the t_agent structure.
*/
void agent_post(t_swarm *x, t_agent *agent)
{
  post("  %i:  Act: %i - Pos: %.2f %.2f %.2f - Vel: %.2f %.2f %.2f - NVel: %.2f - Size: %.2f - Col: %.2f %.2f %.2f %.2f",
    (int)(agent - x->agent_arr), agent->state, agent->pos[0], agent->pos[1], agent->pos[2],
    agent->vel[0], agent->vel[1], agent->vel[2], agent->nvel, agent->size,
    agent->col[0], agent->col[1], agent->col[2], agent->col[3]);

  char str[200];
  char tmp[20];

  // Print the normalized values, 4 per line
  int gene = 0;
  while (gene < GENE_CNT) {

    strcpy(str, "      ");
    for (int i = 0; (i < 4) && (gene < GENE_CNT); i++, gene++) {
      if (i != 0) { strcat(str, " - "); }
      strcat(str, x->gene_repr_arr[gene].sym->s_name);
      strcat(str, ": ");
      snprintf(tmp, 13, "%.2f / ", agent->gene_nval[gene]);
      strcat(str, tmp);
      snprintf(tmp, 10, "%.2f", agent->gene_val[gene]);
      strcat(str, tmp);
    }
    post("%s", str);
  }
  return;
}

/****************************************************************
*  Change the state of an agent.
*/
void agent_state_change(t_swarm *x, t_agent *agent)
{
  TRACE("agent_state_change");

  if (IS_INACTIVE(agent)) {
    jit_object_error((t_object *)x, "agent_state_change:  Invalid action. Agent is inactive.");
    return;
  }

  // Switch through the different state transitions
  switch (agent->state & ~ST_ACTIVE) {

  case ST_DORMANT:
    agent->state = ST_ACTIVE | ST_HATCHING;
    agent->age = 0;
    agent->cntd = rand_time(x->hatch_t, x->hatch_dt);
    break;

  case ST_HATCHING:
    agent->state = ST_ACTIVE | ST_FLYING;
    agent->age = 0;
    agent->cntd = rand_time(x->living_t, x->living_dt);
    agent->size = agent->size_ref;
    for (int c = 0; c < 3; c++) { agent->col[c] = agent->col_ref[c]; }
    break;

  case ST_FLYING:
    agent->state = ST_ACTIVE | ST_DYING;
    agent->age = 0;
    agent->cntd = rand_time(x->dying_t, x->dying_dt);
    for (int v = 0; v < 3; v++) { agent->vel_ref[v] = agent->vel[v]; }
    break;

  case ST_DYING:
    ea_agent_fitness(x, agent);
    ea_agent_store(x, agent);
    ea_agent_stat(x, agent);
    ea_swarm_stat(x);

    t_swr_ind ind1, ind2;
    agent_init(x, agent);
    x->select_fct(x, &ind1, &ind2, x->rank_cnt);
    x->cross_fct(x, agent, ind1, ind2);
    ea_mutate(x, agent);

    agent->state = ST_ACTIVE | ST_DORMANT;
    agent->age = 0;
    agent->cntd = rand_time(x->dorm_t, x->dorm_dt);
    break;

  default:
    jit_object_error((t_object *)x, "agent_state_change:  Invalid state:  %i", agent->state);
    break;
  }
  return;
}

/****************************************************************
*  Iteration function when the agent is hatching.
*
*  Ramp the size from 0, and ramp the color from white.
*/
void agent_iter_hatching(t_swarm *x, t_agent *agent)
{
  t_swr_float r = (t_swr_float)agent->age / (agent->age + agent->cntd);

  agent->size = agent->size_ref * r;
  
  for (int c = 0; c < 3; c++) {
    agent->col[c] = (agent->age < agent->cntd) ? 1.0f :
      2.0f * (agent->col_ref[c] - 1) * r + 2 - agent->col_ref[c];
  }
}

/****************************************************************
*  Iteration function when the agent is flying.
*
*  Apply acceleration to velocity vector.
*/
void agent_iter_flying(t_swarm *x, t_agent *agent)
{
  t_swr_float *X = agent->pos;
  t_swr_float *X2 = NULL;
  t_swr_float *dX = agent->vel;
  t_swr_float *Xa = x->pos_avg;
  t_swr_float *dXa = x->vel_avg;
  t_swr_float  mult, norm;
  t_agent     *agent2 = NULL;

  // Inertia
  mult = agent->gene_val[GENE_INERT_MUL];
  V2D_V_te_S(dX, mult);

  // Cohesion: Seek the average position
  mult = agent->gene_val[GENE_COHES_MUL] / COHES_MUL_BASE;
  V2D_V_pe_S_V_V(dX, mult, Xa, X);

  // Aligment: Match the average velocity
  mult = agent->gene_val[GENE_ALIGN_MUL] / ALIGN_MUL_BASE;
  V2D_V_pe_S_V_V(dX, mult, dXa, dX);

  // Second loop through all the agents
  for (t_swr_ind ind = 0; ind < x->agent_cnt; ind++) {
    agent2 = x->agent_arr + (x->agent_act[ind]);

    if ((agent == agent2) || (agent->state & ST_DORMANT)) { continue; }

    // Calculate the square distance
    X2 = agent2->pos;
    norm = V2D_N_V_V(X, X2);

    // Separation:  If the other agent is close enough
    if (norm < agent->gene_val[GENE_SEPAR_THR]) {
      mult = agent->gene_val[GENE_SEPAR_MUL] / SEPAR_MUL_BASE;
      V2D_V_pe_S_V_V(dX, mult, X, X2);
    }

    // Attraction:  If the other agent is close enough
    if (norm < agent->gene_val[GENE_ATTRA_THR]) {
      mult = agent->gene_val[GENE_ATTRA_MUL] * ATTRA_MUL_BASE / (norm * norm);
      V2D_V_pe_S_V_V(dX, mult, X2, X);
    }

  }   // End second loop through agents

  // Brownian motion
  mult = agent->gene_val[GENE_BROWN_MUL] * BROWN_MUL_BASE;
  if (rand() < RAND_MAX * 0.09f * agent->gene_val[GENE_BROWN_PROB] * agent->gene_val[GENE_BROWN_PROB]) {
    dX[0] += mult * ((t_swr_float)rand() / RAND_MAX - 0.5f);
    dX[1] += mult * ((t_swr_float)rand() / RAND_MAX - 0.5f);
  }

  // Boundary repulsion
  mult = agent->gene_val[GENE_BOUND_MUL] / BOUND_MUL_BASE;
  if (X[0] < agent->gene_val[GENE_BOUND_THR]) {
    dX[0] += mult * (agent->gene_val[GENE_BOUND_THR] - X[0]);
  }
  if (X[0] > x->length_x - agent->gene_val[GENE_BOUND_THR]) {
    dX[0] -= mult * (X[0] - x->length_x + agent->gene_val[GENE_BOUND_THR]);
  }
  if (X[1] < agent->gene_val[GENE_BOUND_THR]) {
    dX[1] += mult * (agent->gene_val[GENE_BOUND_THR] - X[1]);
  }
  if (X[1] > x->length_y - agent->gene_val[GENE_BOUND_THR]) {
    dX[1] -= mult * (X[1] - x->length_y + agent->gene_val[GENE_BOUND_THR]);
  }

  // Clip velocity
  norm = V2D_N(dX);
  if (norm > x->vel_max) {
    mult = x->vel_max / norm;
    V2D_V_te_S(dX, mult);
  }

  // Set other agent attributes
  agent->nvel = min(norm, x->vel_max);

  return;
}

/****************************************************************
*  Iteration function when the agent is dying.
*
*  Ramp the color to black, and ramp the velocity to 0.
*/
void agent_iter_dying(t_swarm *x, t_agent *agent)
{
  t_swr_float r = (t_swr_float)agent->age / (agent->age + agent->cntd);

  for (int c = 0; c < 3; c++) {
    agent->col[c] = -agent->col_ref[c] * r + agent->col_ref[c];
  }

  for (int v = 0; v < 3; v++) {
    agent->vel[v] = agent->vel_ref[v] * (1 - r);
  }
}

/****************************************************************
*  Utility function to get a random time, uniform distribution.
*
*  @param t the average time
*  @param dt +/- time value 
*/
t_swr_cnt rand_time(t_swr_float t, t_swr_float dt)
{
  t_swr_cnt r = (t_swr_cnt)(30 * t * (1 + (2.0f * rand() / RAND_MAX - 1) * dt));

  return max(r, 1);
}