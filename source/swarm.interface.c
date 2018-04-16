//******************************************************************************
//  Header files
//
#include "jit.swarm.h"

//******************************************************************************
//  Reset the swarm system.
//
//  Reset the swarm parameters and the agent parameters.
//
//  @param x A pointer to the t_swarm structure.
//
void swarm_reset(t_swarm* x) {

  TRACE("swarm_reset");

  ea_init_gene_rec(x);
  swarm_init_param(x);

  t_agent* agent;
  t_swr_float time = x->dorm_t + x->hatch_t + x->living_t + x->dying_t;

  for (t_swr_ind ind = 0; ind < x->agent_max; ind++) {
    agent = x->agent_arr + ind;
    if (IS_ACTIVE(agent)) { agent_seed(x, agent, rand_time(3.0f, 0.9f)); }
    else { agent_init(x, agent); }
  }
  return;
}

//******************************************************************************
//  Add an agent.
//
//  @param x A pointer to the t_swarm structure.
//
//  @mess none: add the first available agent
//  @mess int: add the nth agent
//
void swarm_add(t_swarm* x, t_symbol* sym, long argc, t_atom* argv) {

  TRACE("swarm_add");

  t_agent* agent = agent_add(x);
  agent_seed(x, agent, 0);

  t_atom atom[1];
  atom_setlong(atom, x->agent_cnt);
  max_jit_obex_dumpout(x->max_wrap, gensym("count"), 1, atom);

  return;
}

//******************************************************************************
//  Remove an agent.
//
//  @param x A pointer to the t_swarm structure.
//
void swarm_remove(t_swarm* x, t_symbol* sym, long argc, t_atom* argv) {

  TRACE("swarm_remove");

  agent_remove(x, NULL);

  t_atom atom[1];
  atom_setlong(atom, x->agent_cnt);
  max_jit_obex_dumpout(x->max_wrap, gensym("count"), 1, atom);

  return;
}

//******************************************************************************
//  Post information on the swarm object to the Max console.
//
//  @param x A pointer to the t_swarm structure.
//
void swarm_post(t_swarm* x) {

  post("Swarm:  Cnt: %i / %i - X_a: %.2f %.2f %.2f - V_a: %.2f %.2f %.2f - dt: %.2f - V_max: %.2f",
    x->agent_cnt, x->agent_max, x->pos_avg[0], x->pos_avg[1], x->pos_avg[2],
    x->vel_avg[0], x->vel_avg[1], x->vel_avg[2], x->dt, x->vel_max);
  post("  Mult:  Inert: %.2f - Cohes: %.2f - Align: %.2f - Attr: %.2f - Separ: %.2f - Bound: %.2f - Brown: %.2f",
    x->inert_mul, x->cohes_mul, x->align_mul, x->attra_mul, x->separ_mul, x->bound_mul, x->brown_mul);
  post("  Dim: %i %i %i - Thresh:  Attr: %.2f - Separ: %.2f - Bound: %.2f - Brown Prob: %.2f",
    x->length_x, x->length_y, x->length_z, x->attra_thr, x->separ_thr, x->bound_thr, x->brown_prob);

  for (t_swr_ind ind = 0; ind < x->agent_max; ind++) {
    agent_post(x, x->agent_arr + ind);
  }

  return;
}

//******************************************************************************
//  Post the GA parameters.
//
void swarm_ea_post(t_swarm* x) {

  ea_post(x);
}

//******************************************************************************
//  Initialize the object methods.
//
void swarm_init_method(void* swarm_class) {

  jit_class_addmethod(swarm_class, (method)swarm_reset, "reset", 0L);
  jit_class_addmethod(swarm_class, (method)swarm_iter, "iter", 0L);
  jit_class_addmethod(swarm_class, (method)swarm_add, "add", A_GIMME, 0L);
  jit_class_addmethod(swarm_class, (method)swarm_remove, "remove", A_GIMME, 0L);
  jit_class_addmethod(swarm_class, (method)swarm_post, "post", 0L);
  jit_class_addmethod(swarm_class, (method)swarm_ea_post, "post_ea", 0L);
  jit_class_addmethod(swarm_class, (method)ea_test, "test_ea", A_GIMME, 0L);

  jit_class_addmethod(swarm_class, (method)jit_swarm_matrix_calc, "matrix_calc", A_CANT, 0L);
}

//******************************************************************************
//  Initialize the object attributes.
//
void swarm_init_attr(void* swarm_class) {

  // Common flag for the attributes
  long attrflags = JIT_ATTR_GET_DEFER_LOW | JIT_ATTR_SET_USURP_LOW;
  t_jit_object* attr;

  // Attributes:  Basic - Swarm parameters

  CLASS_STICKY_CATEGORY((t_class*)swarm_class, 0, "1 - Swarm parameters");   // attribute category
  CLASS_STICKY_ATTR((t_class*)swarm_class, "basic", 0, "1");                 // attribute tab

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "agent_max", _jit_sym_long, attrflags,
    (method)0L, (method)swarm_set_agent_max, calcoffset(t_swarm, agent_max));
  swarm_init_attr_util(swarm_class, attr, "agent_max", "max agent", "1", 10, 100);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "dt", _jit_sym_float32, attrflags,
    (method)0L, (method)0L, calcoffset(t_swarm, dt));
  swarm_init_attr_util(swarm_class, attr, "dt", "dt", "2", 0, 10);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "vel_max", _jit_sym_float32, attrflags,
    (method)0L, (method)0L, calcoffset(t_swarm, vel_max));
  swarm_init_attr_util(swarm_class, attr, "vel_max", "veloc max", "3", 0, 50);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "size", _jit_sym_float32, attrflags,
    (method)0L, (method)swarm_set_size, calcoffset(t_swarm, size));
  swarm_init_attr_util(swarm_class, attr, "size", "size", "4", 0, 50);

  CLASS_STICKY_CATEGORY_CLEAR((t_class*)swarm_class);
  CLASS_STICKY_ATTR_CLEAR((t_class*)swarm_class, "basic");

  // Attributes: Basic - Evolutionary parameters

  CLASS_STICKY_CATEGORY((t_class*)swarm_class, 0, "2 - Evolutionary parameters");
  CLASS_STICKY_ATTR((t_class*)swarm_class, "basic", 0, "1");

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "fit_xy_mul", _jit_sym_float32, attrflags,
    (method)0L, (method)0L, calcoffset(t_swarm, fit_xy_mul));
  swarm_init_attr_util(swarm_class, attr, "fit_xy_mul", "xy fitness", "1", 0, 10);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "fit_vel_mul", _jit_sym_float32, attrflags,
    (method)0L, (method)0L, calcoffset(t_swarm, fit_vel_mul));
  swarm_init_attr_util(swarm_class, attr, "fit_vel_mul", "vel fitness", "2", 0, 10);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "mut_d", _jit_sym_float32, attrflags,
    (method)0L, (method)0L, calcoffset(t_swarm, mut_d));
  swarm_init_attr_util(swarm_class, attr, "mut_d", "mutation range", "3", 0, 1);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "mut_p", _jit_sym_float32, attrflags,
    (method)0L, (method)0L, calcoffset(t_swarm, mut_p));
  swarm_init_attr_util(swarm_class, attr, "mut_p", "mutation prob", "4", 0, 1);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "select_fi", _jit_sym_long, attrflags,
    (method)0L, (method)swarm_set_select_fct, calcoffset(t_swarm, select_fi));
  swarm_init_attr_util(swarm_class, attr, "select_fi", "selection fct", "5", 1, 3);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "cross_fi", _jit_sym_long, attrflags,
    (method)0L, (method)swarm_set_cross_fct, calcoffset(t_swarm, cross_fi));
  swarm_init_attr_util(swarm_class, attr, "cross_fi", "crossover fct", "6", 1, 3);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "rank_cnt", _jit_sym_long, attrflags,
    (method)0L, (method)0L, calcoffset(t_swarm, rank_cnt));
  swarm_init_attr_util(swarm_class, attr, "rank_cnt", "rank count", "7", 1, 20);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "init_d", _jit_sym_float32, attrflags,
    (method)0L, (method)0L, calcoffset(t_swarm, init_d));
  swarm_init_attr_util(swarm_class, attr, "init_d", "init range", "8", 0, 1);

  CLASS_STICKY_CATEGORY_CLEAR((t_class*)swarm_class);
  CLASS_STICKY_ATTR_CLEAR((t_class*)swarm_class, "basic");

  // Attributes: Basic - Timing

  CLASS_STICKY_CATEGORY((t_class*)swarm_class, 0, "3 - Timing");
  CLASS_STICKY_ATTR((t_class*)swarm_class, "basic", 0, "1");

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "dorm_t", _jit_sym_float32, attrflags,
    (method)0L, (method)0L, calcoffset(t_swarm, dorm_t));
  swarm_init_attr_util(swarm_class, attr, "dorm_t", "dormant t", "1", 0, 10);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "dorm_dt", _jit_sym_float32, attrflags,
    (method)0L, (method)0L, calcoffset(t_swarm, dorm_dt));
  swarm_init_attr_util(swarm_class, attr, "dorm_dt", "dormant dt", "2", 0, 1);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "hatch_t", _jit_sym_float32, attrflags,
    (method)0L, (method)0L, calcoffset(t_swarm, hatch_t));
  swarm_init_attr_util(swarm_class, attr, "hatch_t", "hatching t", "3", 0, 10);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "hatch_dt", _jit_sym_float32, attrflags,
    (method)0L, (method)0L, calcoffset(t_swarm, hatch_dt));
  swarm_init_attr_util(swarm_class, attr, "hatch_dt", "hatching dt", "4", 0, 1);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "living_t", _jit_sym_float32, attrflags,
    (method)0L, (method)0L, calcoffset(t_swarm, living_t));
  swarm_init_attr_util(swarm_class, attr, "living_t", "life t", "5", 0, 600);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "living_dt", _jit_sym_float32, attrflags,
    (method)0L, (method)0L, calcoffset(t_swarm, living_dt));
  swarm_init_attr_util(swarm_class, attr, "living_dt", "life dt", "6", 0, 1);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "dying_t", _jit_sym_float32, attrflags,
    (method)0L, (method)0L, calcoffset(t_swarm, dying_t));
  swarm_init_attr_util(swarm_class, attr, "dying_t", "death t", "7", 0, 10);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "dying_dt", _jit_sym_float32, attrflags,
    (method)0L, (method)0L, calcoffset(t_swarm, dying_dt));
  swarm_init_attr_util(swarm_class, attr, "dying_dt", "death dt", "8", 0, 1);

  CLASS_STICKY_CATEGORY_CLEAR((t_class*)swarm_class);
  CLASS_STICKY_ATTR_CLEAR((t_class*)swarm_class, "basic");

  // Attributes: Basic - Agent modifiers

  CLASS_STICKY_CATEGORY((t_class*)swarm_class, 0, "4 - Agent modifiers");
  CLASS_STICKY_ATTR((t_class*)swarm_class, "basic", 0, "1");

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "inert_mul", _jit_sym_float32, attrflags,
    (method)0L, (method)swarm_set_inert_mul, calcoffset(t_swarm, inert_mul));
  swarm_init_attr_util(swarm_class, attr, "inert_mul", "inertia", "1", 0, 2);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "cohes_mul", _jit_sym_float32, attrflags,
    (method)0L, (method)swarm_set_cohes_mul, calcoffset(t_swarm, cohes_mul));
  swarm_init_attr_util(swarm_class, attr, "cohes_mul", "cohesion", "2", 0, 2);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "align_mul", _jit_sym_float32, attrflags,
    (method)0L, (method)swarm_set_align_mul, calcoffset(t_swarm, align_mul));
  swarm_init_attr_util(swarm_class, attr, "align_mul", "alignment", "3", 0, 2);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "attra_mul", _jit_sym_float32, attrflags,
    (method)0L, (method)swarm_set_attra_mul, calcoffset(t_swarm, attra_mul));
  swarm_init_attr_util(swarm_class, attr, "attra_mul", "attraction", "4", 0, 2);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "attra_thr", _jit_sym_float32, attrflags,
    (method)0L, (method)swarm_set_attra_thr, calcoffset(t_swarm, attra_thr));
  swarm_init_attr_util(swarm_class, attr, "attra_thr", "attra thresh", "5", 0, 2);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "separ_mul", _jit_sym_float32, attrflags,
    (method)0L, (method)swarm_set_separ_mul, calcoffset(t_swarm, separ_mul));
  swarm_init_attr_util(swarm_class, attr, "separ_mul", "separation", "6", 0, 2);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "separ_thr", _jit_sym_float32, attrflags,
    (method)0L, (method)swarm_set_separ_thr, calcoffset(t_swarm, separ_thr));
  swarm_init_attr_util(swarm_class, attr, "separ_thr", "separ thresh", "7", 0, 2);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "bound_mul", _jit_sym_float32, attrflags,
    (method)0L, (method)swarm_set_bound_mul, calcoffset(t_swarm, bound_mul));
  swarm_init_attr_util(swarm_class, attr, "bound_mul", "boundaries", "8", 0, 2);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "bound_thr", _jit_sym_float32, attrflags,
    (method)0L, (method)swarm_set_bound_thr, calcoffset(t_swarm, bound_thr));
  swarm_init_attr_util(swarm_class, attr, "bound_thr", "bound thresh", "9", 0, 2);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "brown_mul", _jit_sym_float32, attrflags,
    (method)0L, (method)swarm_set_brown_mul, calcoffset(t_swarm, brown_mul));
  swarm_init_attr_util(swarm_class, attr, "brown_mul", "brownian", "10", 0, 2);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "brown_prob", _jit_sym_float32, attrflags,
    (method)0L, (method)swarm_set_brown_prob, calcoffset(t_swarm, brown_prob));
  swarm_init_attr_util(swarm_class, attr, "brown_prob", "brown prob", "11", 0, 2);

  CLASS_STICKY_CATEGORY_CLEAR((t_class*)swarm_class);
  CLASS_STICKY_ATTR_CLEAR((t_class*)swarm_class, "basic");

  // Attributes: Basic - Display

  CLASS_STICKY_CATEGORY((t_class*)swarm_class, 0, "5 - Display");
  CLASS_STICKY_ATTR((t_class*)swarm_class, "basic", 0, "1");

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "length_x", _jit_sym_long, attrflags,
    (method)0L, (method)0L, calcoffset(t_swarm, length_x));
  swarm_init_attr_util(swarm_class, attr, "length_x", "length x", "1", 100, 1920);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "length_y", _jit_sym_long, attrflags,
    (method)0L, (method)0L, calcoffset(t_swarm, length_y));
  swarm_init_attr_util(swarm_class, attr, "length_y", "length y", "2", 100, 1080);

  attr = (t_jit_object*)jit_object_new(_jit_sym_jit_attr_offset, "length_z", _jit_sym_long, attrflags,
    (method)0L, (method)0L, calcoffset(t_swarm, length_z));
  swarm_init_attr_util(swarm_class, attr, "length_z", "length z", "3", 100, 1080);

  CLASS_STICKY_CATEGORY_CLEAR((t_class*)swarm_class);
  CLASS_STICKY_ATTR_CLEAR((t_class*)swarm_class, "basic");

  return;
}

//******************************************************************************
//  Utility function to initialize the object attributes.
//
void swarm_init_attr_util(void* swarm_class, t_jit_object* attr,
  const char* name, const char* label, const char* order, float min_, float max_) {

  jit_class_addattr(swarm_class, attr);
  CLASS_ATTR_LABEL((t_class*)swarm_class, name, 0, label);           // label
  CLASS_ATTR_ORDER((t_class*)swarm_class, name, 0, order);           // order
  CLASS_ATTR_FILTER_CLIP((t_class*)swarm_class, name, min_, max_);   // min-max filter
  CLASS_ATTR_SAVE((t_class*)swarm_class, name, 0);                   // save with patcher
  CLASS_ATTR_SELFSAVE((t_class*)swarm_class, name, 0);               // display as saved

  return;
}

//******************************************************************************
//  Setter method for the agents_max attribute.
//
//  Custom setter method because it resizes the matrices and arrays.
//
//  @param x A pointer to the t_swarm structure.
//  @param attr A pointer to the attribute.
//  @param argc The length of the list of atom arguments.
//  @param argv A pointer to the list of atom arguments.
//
t_jit_err swarm_set_agent_max(t_swarm* x, void* attr, long argc, t_atom* argv) {

  // Get the new maximum number
  t_swr_ind max_ = AGENT_MAX_DEF;
  if (argc && argv) {
    max_ = (t_swr_ind)jit_atom_getlong(argv);
  }  // clipping at filter level

     // Get the MOP wrapper as an adornment
  void* mop = jit_class_adornment_get(_jit_swarm_class, _jit_sym_jit_mop);

  // Get the MOP wrapper outputs and change the maximum dimension
  void* output1 = jit_object_method(mop, _jit_sym_getoutput, 1);
  void* output2 = jit_object_method(mop, _jit_sym_getoutput, 2);
  void* output3 = jit_object_method(mop, _jit_sym_getoutput, 3);
  void* output4 = jit_object_method(mop, _jit_sym_getoutput, 4);

  jit_attr_setlong(output1, _jit_sym_maxdim, max_);
  jit_attr_setlong(output2, _jit_sym_maxdim, max_);
  jit_attr_setlong(output3, _jit_sym_maxdim, max_);
  jit_attr_setlong(output4, _jit_sym_maxdim, max_);

  swarm_new_agents(x, min(x->agent_cnt, max_), max_);

  return JIT_ERR_NONE;
}

//******************************************************************************
//  Setter method for the size attribute.
//
//  Updates the corresponding agent values.
//
t_jit_err swarm_set_size(t_swarm* x, void* attr, long argc, t_atom* argv) {

  t_swr_float size_raw, size;    // The modifier for all agents, between 0 and 2
  if (argc && argv && ((atom_gettype(argv) == A_LONG) || (atom_gettype(argv) == A_FLOAT))) {
    size_raw = (t_swr_float)jit_atom_getfloat(argv);
    size = CLAMP(size_raw, 0.0f, 100.0f);

  } else { return JIT_ERR_INVALID_INPUT;
}

  // Set the structure attribute
  x->size = size;

  // Loop through all the agents
  for (t_swr_ind ind = 0; ind < x->agent_max; ind++) {
    x->agent_arr[ind].size_ref = x->size;
    x->agent_arr[ind].size = x->size;
  }

  return JIT_ERR_NONE;
}

//******************************************************************************
//  Setter method for the selection function.
//
t_jit_err swarm_set_select_fct(t_swarm* x, void* attr, long argc, t_atom* argv) {

  t_swr_ind fct = (t_swr_ind)atom_getlong(argv);

  switch (fct) {
  case 1: x->select_fct = ea_select_rank; break;
  case 2: x->select_fct = ea_select_roulette; break;
  case 3: x->select_fct = ea_select_tournament; break;
  default:
    jit_object_error((t_object*)x, "swarm_set_select_fct:  Invalid value:  1, 2 or 3 expected.");
    break;
  }

  x->select_fi = fct;

  return JIT_ERR_NONE;
}

//******************************************************************************
//  Setter method for the crossover function.
//
t_jit_err swarm_set_cross_fct(t_swarm* x, void* attr, long argc, t_atom* argv) {

  t_swr_ind fct = (t_swr_ind)atom_getlong(argv);

  switch (fct) {
  case 1: x->cross_fct = ea_cross_single; break;
  case 2: x->cross_fct = ea_cross_double; break;
  case 3: x->cross_fct = ea_cross_uniform; break;
  default:
    jit_object_error((t_object*)x, "swarm_set_cross_fct:  Invalid value:  1, 2 or 3 expected.");
    break;
  }

  x->cross_fi = fct;

  return JIT_ERR_NONE;
}

//******************************************************************************
//  Setter method for the inert_mul attribute.
//
//  Updates the corresponding agent values.
//
t_jit_err swarm_set_inert_mul(t_swarm* x, void* attr, long argc, t_atom* argv) {

  return swarm_set_gene_val(x, GENE_INERT_MUL, &(x->inert_mul), argc, argv);
}

//******************************************************************************
//  Setter method for the cohes_mul attribute.
//
t_jit_err swarm_set_cohes_mul(t_swarm* x, void* attr, long argc, t_atom* argv) {

  return swarm_set_gene_val(x, GENE_COHES_MUL, &(x->cohes_mul), argc, argv);
}

//******************************************************************************
//  Setter method for the align_mul attribute.
//
t_jit_err swarm_set_align_mul(t_swarm* x, void* attr, long argc, t_atom* argv) {

  return swarm_set_gene_val(x, GENE_ALIGN_MUL, &(x->align_mul), argc, argv);
}

//******************************************************************************
//  Setter method for the attra_mul attribute.
//
t_jit_err swarm_set_attra_mul(t_swarm* x, void* attr, long argc, t_atom* argv) {

  return swarm_set_gene_val(x, GENE_ATTRA_MUL, &(x->attra_mul), argc, argv);
}

//******************************************************************************
//  Setter method for the attra_thr attribute.
//
t_jit_err swarm_set_attra_thr(t_swarm* x, void* attr, long argc, t_atom* argv) {

  return swarm_set_gene_val(x, GENE_ATTRA_THR, &(x->attra_thr), argc, argv);
}

//******************************************************************************
//  Setter method for the separ_mul attribute.
//
t_jit_err swarm_set_separ_mul(t_swarm* x, void* attr, long argc, t_atom* argv) {

  return swarm_set_gene_val(x, GENE_SEPAR_MUL, &(x->separ_mul), argc, argv);
}

//******************************************************************************
//  Setter method for the separ_thr attribute.
//
t_jit_err swarm_set_separ_thr(t_swarm* x, void* attr, long argc, t_atom* argv) {

  return swarm_set_gene_val(x, GENE_SEPAR_THR, &(x->separ_thr), argc, argv);
}

//******************************************************************************
//  Setter method for the bound_mul attribute.
//
t_jit_err swarm_set_bound_mul(t_swarm* x, void* attr, long argc, t_atom* argv) {

  return swarm_set_gene_val(x, GENE_BOUND_MUL, &(x->bound_mul), argc, argv);
}

//******************************************************************************
//  Setter method for the bound_thr attribute.
//
t_jit_err swarm_set_bound_thr(t_swarm* x, void* attr, long argc, t_atom* argv) {

  return swarm_set_gene_val(x, GENE_BOUND_THR, &(x->bound_thr), argc, argv);
}

//******************************************************************************
//  Setter method for the brown_mul attribute.
//
t_jit_err swarm_set_brown_mul(t_swarm* x, void* attr, long argc, t_atom* argv) {

  return swarm_set_gene_val(x, GENE_BROWN_MUL, &(x->brown_mul), argc, argv);
}

//******************************************************************************
//  Setter method for the brown_prob attribute.
//
t_jit_err swarm_set_brown_prob(t_swarm* x, void* attr, long argc, t_atom* argv) {

  return swarm_set_gene_val(x, GENE_BROWN_PROB, &(x->brown_prob), argc, argv);
}

//******************************************************************************
//  Setter method for the attributes, to update the corresponding agent values.
//
//  Each agent has a [0, 2] normalized value (nval) used for the genetic algorithm.
//  The swarm as a whole has a [0, 2] modifier value (mod).
//  A modified normalized value (nval_mod) is scaled:
//  - mod = 0 --> nval = 0
//  - mod = 1 --> nval unchanged
//  - mod = 2 --> nval = 2
//
//  Finally a scaled value is calculated, for the swarm algorithm,
//  according to the scaling parameters defined in the gene representation.
//
//  @param x A pointer to the t_swarm structure.
//  @param gene_ind the index of the gene, defined in an enum.
//  @param pmod A pointer to the modifier attribute in the structure.
//  @param argc The length of the list of atom arguments.
//  @param argv A pointer to the list of atom arguments.
//
t_jit_err swarm_set_gene_val(t_swarm* x, e_gene_ind gene_ind, t_swr_float* pmod, long argc, t_atom* argv) {

  t_swr_float mod_raw, mod;    // The modifier for all agents, between 0 and 2
  if (argc && argv && ((atom_gettype(argv) == A_LONG) ||(atom_gettype(argv) == A_FLOAT))) {
    mod_raw = (t_swr_float)jit_atom_getfloat(argv);
    mod = CLAMP(mod_raw, 0.0f, 2.0f);

  } else { return JIT_ERR_INVALID_INPUT;
}

  // Set the structure attribute
  *pmod = mod;

  // Calculate the gene values
  t_swr_float nval;       // the unmodified normalized value, between 0 and 2
  t_swr_float nval_mod;   // the modified normalized value, between 0 and 2
  t_gene_repr* gene_repr = x->gene_repr_arr + gene_ind;    // the gene representation

  // Loop through all the agents
  for (t_swr_ind ind = 0; ind < x->agent_max; ind++) {
    nval = x->agent_arr[ind].gene_nval[gene_ind];

    nval_mod = (mod < 1) ? nval * mod :   // mod < 1: Scale between 0 and nval
      (2 - nval) * mod + 2 * nval - 2;    // mod > 1: Scale between nval and 2

    x->agent_arr[ind].gene_val[gene_ind] = gene_repr->scale_fct(nval_mod, gene_repr);
  }

  return JIT_ERR_NONE;
}
