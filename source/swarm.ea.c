/****************************************************************
*  Header files
*/
#include "jit.swarm.h"

/****************************************************************
*  Initialize the gene representation array for the evolutionary algorithm.
*/
void ea_init_gene_repr(t_swarm *x)
{
  ea_gene_repr_set(x, GENE_INERT_MUL, 0.9f, 1, 1.2f, ea_fct_scale_exp, true, gensym("inert"));
  ea_gene_repr_set(x, GENE_COHES_MUL, 0, 1, 10, ea_fct_scale_exp, true, gensym("cohes"));
  ea_gene_repr_set(x, GENE_ALIGN_MUL, 0, 1, 10, ea_fct_scale_exp, true, gensym("align"));
  ea_gene_repr_set(x, GENE_ATTRA_MUL, 0, 1, 10, ea_fct_scale_exp, true, gensym("attra"));
  ea_gene_repr_set(x, GENE_SEPAR_MUL, 0, 1, 10, ea_fct_scale_exp, true, gensym("separ"));
  ea_gene_repr_set(x, GENE_BOUND_MUL, 0, 1, 10, ea_fct_scale_exp, true, gensym("bound"));
  ea_gene_repr_set(x, GENE_BROWN_MUL, 0, 1, 10, ea_fct_scale_exp, true, gensym("brown"));

  ea_gene_repr_set(x, GENE_ATTRA_THR, 10, ATTRA_THR_BASE, 400, ea_fct_scale_exp, true, gensym("att_t"));
  ea_gene_repr_set(x, GENE_SEPAR_THR, 10, SEPAR_THR_BASE, 400, ea_fct_scale_exp, true, gensym("sep_t"));
  ea_gene_repr_set(x, GENE_BOUND_THR, 10, BOUND_THR_BASE, 400, ea_fct_scale_exp, true, gensym("bou_t"));
  ea_gene_repr_set(x, GENE_BROWN_PROB, 0, BROWN_PROB_BASE,  1, ea_fct_scale_exp, true, gensym("bro_p"));

  return;
}

/****************************************************************
*  Initialize the gene representation array for the evolutionary algorithm.
*/
void ea_init_gene_rec(t_swarm *x)
{
  t_gene_bank *gene_rec;
  for (int g = 0; g < GENE_REC_CNT; g++) {
    gene_rec = x->gene_bank_arr + g;

    gene_rec->fitness = 0.000001f;    // Starting at 0 causes problem with roulette
    gene_rec->fit_xy  = 0.000001f;
    gene_rec->fit_vel = 0.000001f;

    x->init_fct = ea_fct_init_lin;
    ea_init_gene_nval(gene_rec->gene_nval, x);
    
    x->gene_bank_sort[g] = g;
  }

  x->gene_bank_ptr = x->gene_bank_arr;

  return;
}

/****************************************************************
*  Calculate the agent's fitness.
*/
void ea_agent_fitness(t_swarm *x, t_agent *agent)
{
  TRACE("ea_agent_fitness");

  t_swr_float fit_xy = 0;
  t_swr_float fit_vel = 0;
  t_swr_float avg, fitness;

  // Calculate the square of the euclidian distance to a uniform spatial distribution
  avg = (t_swr_float)agent->age / (BIN_CNT_X * BIN_CNT_Y);
  for (int XY = 0; XY < BIN_CNT_X * BIN_CNT_Y; XY++) {
    fit_xy += (agent->bin_xy[XY] - avg) * (agent->bin_xy[XY] - avg);
  }
  // Penalize out-of-bounds count
  fit_xy += agent->bin_xy[BIN_CNT_X * BIN_CNT_Y] * agent->bin_xy[BIN_CNT_X * BIN_CNT_Y];  
  fit_xy = sqrtf(1 / (1 + fit_xy / (BIN_CNT_X * BIN_CNT_Y + 1)));

  // Calculate the square of the euclidian distance to a uniform velocity distribution
  avg = (t_swr_float)agent->age / (BIN_CNT_VEL);
  for (int v = 0; v < BIN_CNT_VEL; v++) {
    fit_vel += (agent->bin_vel[v] - avg) * (agent->bin_vel[v] - avg);
  }
  fit_vel = sqrtf(1 / (1 + fit_vel / BIN_CNT_VEL));

  // Calculate a linear combination of both components, avoid 0 division
  fitness = (x->fit_xy_mul * fit_xy + x->fit_vel_mul * fit_vel) /
    ((x->fit_xy_mul + x->fit_vel_mul) == 0 ? 1 :
      (x->fit_xy_mul + x->fit_vel_mul));

  // Set the values
  agent->fitness = fitness;
  agent->fit_xy = fit_xy;
  agent->fit_vel = fit_vel;

  // Send a list out of dumpout outlet
  t_atom argv[3];
  atom_setfloat(argv + 0, fitness);
  atom_setfloat(argv + 1, fit_xy);
  atom_setfloat(argv + 2, fit_vel);
  max_jit_obex_dumpout(x->max_wrap, gensym("fitness"), 3, argv);

  return;
}

/****************************************************************
*  Store the agent after it died.
*/
void ea_agent_store(t_swarm *x, t_agent *agent)
{
  TRACE("ea_agent_store");

  t_gene_bank *rec = x->gene_bank_arr;
  t_swr_ind *sort = x->gene_bank_sort;
  t_swr_ind rec_ind = (t_swr_ind)(x->gene_bank_ptr - x->gene_bank_arr);

  t_swr_ind ins = 0;   // important if new fitness lower than all others
  t_swr_ind del = 0xFFFF;

  // We need to insert into and shift the sorting array
  // First pass to get the indexes at which to insert and delete
  for (t_swr_ind ind = 0; ind < GENE_REC_CNT; ind++) {

    if (rec[sort[ind]].fitness > agent->fitness) { ins = ind + 1; }
    if (sort[ind] == rec_ind) { del = ind; }
  }

  if (del == 0xFFFF) { post("ERROR:  record index not found"); return; }

  // Second pass to shift upward or downward
  if (ins <= del) {
    for (t_swr_ind ind = del; ind > ins; ind--) { sort[ind] = sort[ind - 1]; }
    sort[ins] = rec_ind;
  } else {
    for (t_swr_ind ind = del; ind < ins - 1; ind++) { sort[ind] = sort[ind + 1]; }
    sort[ins - 1] = rec_ind;
  }

  // Copy the agent gene values into the record
  rec = x->gene_bank_ptr;
  rec->fitness = agent->fitness;
  rec->fit_xy = agent->fit_xy;
  rec->fit_vel = agent->fit_vel;

  ea_copy_gene_nval(agent->gene_nval, rec->gene_nval, x);

  // Iterate and fold when necessary the gene record pointer
  x->gene_bank_ptr++;
  x->gene_bank_ptr = (x->gene_bank_ptr == x->gene_bank_arr + GENE_REC_CNT) ?
    x->gene_bank_arr :
    x->gene_bank_ptr;

  return;
}

/****************************************************************
*  Send statistics about the agent through the dumpout outlet.
*/
void ea_agent_stat(t_swarm *x, t_agent *agent)
{
  TRACE("ea_agent_stat");

  t_atom argv[max(BIN_CNT_X * BIN_CNT_Y, BIN_CNT_VEL)];
  t_atom *atom = argv;

  // Get the maximum count in the spatial bin
  t_swr_cnt max_cnt = 0;
  for (int XY = 0; XY < BIN_CNT_X * BIN_CNT_Y; XY++) {
    max_cnt = max(max_cnt, agent->bin_xy[XY]);
  }

  // Normalize the spatial bins, reverse the vertical order
  for (int Y = BIN_CNT_Y - 1; Y >= 0; Y--) {
    for (int X = 0; X < BIN_CNT_X; X++) {
      atom_setfloat(atom++, (double)agent->bin_xy[X + Y * BIN_CNT_X] / max_cnt);
    }
  }

  // Send a list out of dumpout outlet
  max_jit_obex_dumpout(x->max_wrap, gensym("bin_xy"), BIN_CNT_X * BIN_CNT_Y, argv);

  // Get the maximum count in the velocity bins
  max_cnt = 0;
  for (int v = 0; v < BIN_CNT_VEL; v++) {
    max_cnt = max(max_cnt, agent->bin_vel[v]);
  }

  // Normalize the velocity bins
  for (int v = 0; v < BIN_CNT_VEL; v++) {
    atom_setfloat(argv + v, (double)agent->bin_vel[v] / max_cnt);
  }

  // Send a list out of dumpout outlet
  max_jit_obex_dumpout(x->max_wrap, gensym("bin_vel"), BIN_CNT_VEL, argv);

  return;
}

/****************************************************************
*  Send statistics about the swarm through the dumpout outlet.
*/
void ea_swarm_stat(t_swarm *x)
{
  TRACE("ea_swarm_stat");

  t_agent *agent;
  t_swr_float nval_avg[GENE_CNT];
  t_swr_float nval_std[GENE_CNT];
  t_atom argv[GENE_CNT];

  // Initialize to 0
  for (t_swr_ind g = 0; g < GENE_CNT; g++) {
    nval_avg[g] = 0;
    nval_std[g] = 0;
  }
  
  // Calculate the average
  for (t_swr_ind ind = 0; ind < x->agent_cnt; ind++) {
    agent = x->agent_arr + (x->agent_act[ind]);
    for (t_swr_ind g = 0; g < GENE_CNT; g++) { nval_avg[g] += agent->gene_nval[g]; }
  }

  // Divide, then set in the array of atoms
  for (t_swr_ind g = 0; g < GENE_CNT; g++) {
    nval_avg[g] /= x->agent_cnt;
    atom_setfloat(argv + g, nval_avg[g]);
  }

  // Send a list out of dumpout outlet
  max_jit_obex_dumpout(x->max_wrap, gensym("average"), GENE_CNT, argv);

  // Calculate the standard deviation
  for (t_swr_ind ind = 0; ind < x->agent_cnt; ind++) {
    agent = x->agent_arr + (x->agent_act[ind]);
    for (t_swr_ind g = 0; g < GENE_CNT; g++) {
      nval_std[g] += (agent->gene_nval[g] - nval_avg[g]) * (agent->gene_nval[g] - nval_avg[g]);
    }
  }

  // Divide, take the square root, then set in the array of atoms
  for (t_swr_ind g = 0; g < GENE_CNT; g++) {
    nval_std[g] = sqrtf(nval_std[g] / x->agent_cnt);
    atom_setfloat(argv + g, nval_std[g]);
  }

  // Send a list out of dumpout outlet
  max_jit_obex_dumpout(x->max_wrap, gensym("std_dev"), GENE_CNT, argv);

  return;
}

/****************************************************************
*  Selection by rank.
*
*  Select 1 or 2 randomly out of a set number of the best ones.
*
*  @param ind1 A pointer to the first index.
*  @param ind2 A pointer to the second index. Pass NULL if only one needed.
*  @param out_of The number of best agents to draw from.
*/
void ea_select_rank(t_swarm *x, t_swr_ind *ind1, t_swr_ind *ind2, t_swr_ind out_of)
{
  *ind1 = rand() / (RAND_MAX / out_of + 1);

  // If a second index is requested
  if (ind2) {
    if (out_of > 1) {
      *ind2 = rand() / (RAND_MAX / (out_of - 1) + 1);
      *ind2 += (*ind2 >= *ind1) ? 1 : 0;
    } else { *ind2 = 0; }
  }
  return;
}

/****************************************************************
*  Selection by roulette.
*
*  Select 1 or 2 using the fitnesses as a probability distribution.
*
*  @param ind1 A pointer to the first index.
*  @param ind2 A pointer to the second index. Pass NULL if only one needed.
*  @param N Not used, necessary for function pointer typedef.
*/
void ea_select_roulette(t_swarm *x, t_swr_ind *ind1, t_swr_ind *ind2, t_swr_ind N)
{
  *ind1 = ea_select_roulette_util(x, NO_EXCLUDE);
  if (ind2) { *ind2 = ea_select_roulette_util(x, *ind1); }

  return;
}

/****************************************************************
*  Selection by tournament.
*
*  Select 2 by roulette and keep the best. Repeat if necessary.
*
*  @param ind1 A pointer to the first index.
*  @param ind2 A pointer to the second index. Pass NULL if only one needed.
*  @param N Not used, necessary for function pointer typedef.
*/
void ea_select_tournament(t_swarm *x, t_swr_ind *ind1, t_swr_ind *ind2, t_swr_ind N)
{
  t_gene_bank *rec = x->gene_bank_arr;
  t_swr_ind *sort = x->gene_bank_sort;

  t_swr_ind i, j;

  // Choose two agents by roulette
  i = ea_select_roulette_util(x, NO_EXCLUDE);
  j = ea_select_roulette_util(x, NO_EXCLUDE);

  // And keep the best
  *ind1 = min(i, j);

  // If a second agent is requested
  if (ind2) {
    i = ea_select_roulette_util(x, *ind1);
    j = ea_select_roulette_util(x, *ind1);
    *ind2 = min(i, j);
  }
  return;
}

/****************************************************************
*  Utility function to selection by roulette.
*
*  Chooses from the whole population,
*  with the fitnesses used as a probability distribution.
*
*  @param exclude An index to exclude from the selection.
*    Use NO_EXCLUDE if no exclusion.
*
*  @return An index.
*/
t_swr_ind ea_select_roulette_util(t_swarm *x, t_swr_ind exclude)
{
  t_gene_bank *rec = x->gene_bank_arr;
  t_swr_ind *sort = x->gene_bank_sort;

  // Calculate the sum of the fitness
  t_swr_float sum = 0;
  for (t_swr_ind ind = 0; ind < GENE_REC_CNT; ind++) {

    // Exclude index from sum if necessary
    sum += (ind != exclude) ? rec[sort[ind]].fitness : 0.0f;
  }

  // Draw a random value between 0 and sum
  t_swr_float r = sum * rand() / (RAND_MAX + 1.0f);

  // Select within whole population with the fitness as probability distribution
  t_swr_ind ind = 0;
  sum = 0;
  do {

    // Exclude index from choice if necessary
    sum += (ind != exclude) ? rec[sort[ind]].fitness : 0.0f;
    ind++;
  } while (sum < r);

  return (ind - 1);
}

/****************************************************************
*  Single crossover of two records.
*
*  @param agent A pointer to the agent.
*  @param ind1 The index of the first record (for sorted array).
*  @param ind2 The index of the first record (for sorted array).
*/
void ea_cross_single(t_swarm *x, t_agent *agent, t_swr_ind ind1, t_swr_ind ind2)
{
  t_gene_bank *rec = x->gene_bank_arr;
  t_swr_ind *sort = x->gene_bank_sort;
  t_swr_ind g1;

  // Choose a gene index at random
  g1 = rand() / (RAND_MAX / (GENE_CNT + 1) + 1);

  // Crossover switching twice at g1 and g2
  for (int g = 0; g < GENE_CNT; g++) {
    agent->gene_nval[g] = (g < g1) ?
      rec[sort[ind1]].gene_nval[g] : rec[sort[ind2]].gene_nval[g];
  }

  // Update the calculated gene values
  ea_calc_gene_val(agent->gene_nval, agent->gene_val, x);

  return;
}

/****************************************************************
*  Double crossover of two records.
*
*  @param agent A pointer to the agent.
*  @param ind1 The index of the first record (for sorted array).
*  @param ind2 The index of the first record (for sorted array).
*/
void ea_cross_double(t_swarm *x, t_agent *agent, t_swr_ind ind1, t_swr_ind ind2)
{
  t_gene_bank *rec = x->gene_bank_arr;
  t_swr_ind *sort = x->gene_bank_sort;
  t_swr_ind g1, g2, tmp;

  // Choose two non equal gene indexes at random
  g1 = rand() / (RAND_MAX / (GENE_CNT + 1) + 1);
  g2 = rand() / (RAND_MAX / GENE_CNT + 1);
  g2 += (g2 >= g1) ? 1 : 0;

  // Sort the two numbers
  tmp = g1;
  g1 = min(g1, g2);
  g2 = max(tmp, g2);

  // Crossover switching once at g1
  for (int g = 0; g < GENE_CNT; g++) {
    agent->gene_nval[g] = ((g < g1) || (g >= g2)) ?
      rec[sort[ind1]].gene_nval[g] : rec[sort[ind2]].gene_nval[g];
  }

  // Update the calculated gene values
  ea_calc_gene_val(agent->gene_nval, agent->gene_val, x);

  return;
}

/****************************************************************
*  Uniform crossover of two records.
*
*  @param agent A pointer to the agent.
*  @param ind1 The index of the first record (for sorted array).
*  @param ind2 The index of the first record (for sorted array).
*/
void ea_cross_uniform(t_swarm *x, t_agent *agent, t_swr_ind ind1, t_swr_ind ind2)
{
  t_gene_bank *rec = x->gene_bank_arr;
  t_swr_ind *sort = x->gene_bank_sort;
  t_swr_ind r;

  // Crossover switching randomly at every index
  for (int g = 0; g < GENE_CNT; g++) {
    r = rand() / (RAND_MAX / 2 + 1);
    agent->gene_nval[g] = r ?
      rec[sort[ind1]].gene_nval[g] : rec[sort[ind2]].gene_nval[g];
  }

  // Update the calculated gene values
  ea_calc_gene_val(agent->gene_nval, agent->gene_val, x);

  return;
}

/****************************************************************
*  Mutate the agent.
*/
void ea_mutate(t_swarm *x, t_agent *agent)
{
  t_swr_float nval;

  // Loop through the genes and mutate
  for (int g = 0; g < GENE_CNT; g++) {

    nval = agent->gene_nval[g];
    nval += ((t_swr_float)rand() / RAND_MAX < x->mut_p) ?       // if mutation
      (rand() / (t_swr_float)RAND_MAX * 2 - 1) * x->mut_d : 0;  // range of mutation

    nval = CLAMP(nval, 0, 2);
    agent->gene_nval[g] = nval;
  }

  // Update the calculated gene values
  ea_calc_gene_val(agent->gene_nval, agent->gene_val, x);

  return;
}

/****************************************************************
*  Initialize the normalized genome values at random.
*/
void ea_init_gene_nval(t_swr_float *nval, t_swarm *x)
{
  for (int g = 0; g < GENE_CNT; g++) { nval[g] = x->init_fct(x); }
}

/****************************************************************
*  Calculate the scaled genome values.
*
*  Uses the scaling function defined in the gene representation.
*/
void ea_calc_gene_val(t_swr_float *nval, t_swr_float *val, t_swarm *x)
{
  for (int g = 0; g < GENE_CNT; g++) {
    val[g] = x->gene_repr_arr[g].scale_fct(nval[g], x->gene_repr_arr + g);
  }
}

/****************************************************************
*  Copy the normalized genome values.
*/
void ea_copy_gene_nval(t_swr_float *nval1, t_swr_float *nval2, t_swarm *x)
{
  for (int g = 0; g < GENE_CNT; g++) { nval2[g] = nval1[g]; }
}

/****************************************************************
*  Set the parameters for a gene representation.
*/
void ea_gene_repr_set(t_swarm *x, int g, t_swr_float min, t_swr_float mid, t_swr_float max, t_ea_scale scale, t_bool var, t_symbol *sym)
{
  x->gene_repr_arr[g].min = min;
  x->gene_repr_arr[g].mid = mid;
  x->gene_repr_arr[g].max = max;
  x->gene_repr_arr[g].scale_fct = scale;
  x->gene_repr_arr[g].variable = var;
  x->gene_repr_arr[g].sym = sym;

  return;
}

/****************************************************************
*  Post the GA parameters.
*/
void ea_post(t_swarm *x)
{
  post("EA:  cnt: %i - init_d: %.2f - mut_d: %.2f - mut_p: %.2f",
    GENE_CNT, x->init_d, x->mut_d, x->mut_p);

  post("Gene representation:");
  t_gene_repr *gene_repr;
  for (int g = 0; g < GENE_CNT; g++) {
    gene_repr = x->gene_repr_arr + g;
    post("  %i - %s:  min: %.2f - mid: %.2f - max: %.2f - variable: %i",
      g, gene_repr->sym->s_name, gene_repr->min, gene_repr->mid, gene_repr->max, gene_repr->variable);
  }

  post("Gene record:  Normalized gene values - Ptr:  %i", x->gene_bank_ptr - x->gene_bank_arr);
  t_gene_bank *gene_rec;
  for (int g = 0; g < GENE_REC_CNT; g++) {
    gene_rec = x->gene_bank_arr + g;
    post("  %i:  %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f", g,
      gene_rec->gene_nval[0], gene_rec->gene_nval[1], gene_rec->gene_nval[2], gene_rec->gene_nval[3],
      gene_rec->gene_nval[4], gene_rec->gene_nval[5], gene_rec->gene_nval[6], gene_rec->gene_nval[7],
      gene_rec->gene_nval[8], gene_rec->gene_nval[9], gene_rec->gene_nval[10]);
  }

  t_swr_float sum = 0;
  post("Gene bank:  Sorted fitnesses");
  for (int g = 0; g < GENE_REC_CNT; g++) {
    gene_rec = x->gene_bank_arr + x->gene_bank_sort[g];
    sum += gene_rec->fitness;
    post("  %i - %i:  fit: %.4f - sum: %f - fit_xy: %.4f - fit_vel: %.4f",
      g, x->gene_bank_sort[g], gene_rec->fitness, sum, gene_rec->fit_xy, gene_rec->fit_vel);
  }
  return;
}

/****************************************************************
*  Scaling functions
*  
*  Calculate a scaled gene value from an nvalue between 0 and 2.
*  Used to convert the gene values into parameters.
*/
t_swr_float ea_fct_scale_plin(t_swr_float nval, t_gene_repr *gene_repr)
{
  if (nval <= 0) { return gene_repr->min; }
  if (nval >= 2) { return gene_repr->max; }

  if (nval <= 1) { return ((gene_repr->mid - gene_repr->min) * nval + gene_repr->min); }
  return ((gene_repr->max - gene_repr->mid) * nval - gene_repr->max + 2 * gene_repr->mid);
}

t_swr_float ea_fct_scale_exp(t_swr_float nval, t_gene_repr *gene_repr)
{
  if (nval <= 0) { return gene_repr->min; }
  if (nval >= 2) { return gene_repr->max; }

  t_swr_float u, v, a;
  u = gene_repr->mid - gene_repr->min;
  v = gene_repr->max - gene_repr->mid;

  if (!(u * v)) { return gene_repr->mid; }
  if (u == v) { return u * nval + gene_repr->min; }

  a = u * u / (v - u);
  return a * powf(v / u, nval) - a + gene_repr->min;
}

/****************************************************************
*  Initialization functions
*
*  Get a random nvalue between 0 and 2.
*/
t_swr_float ea_fct_init_mid(t_swarm *x)
{
  return 1;
}

t_swr_float ea_fct_init_lin(t_swarm *x)
{
  return (1 + x->init_d * (2.0f * rand() / RAND_MAX - 1));
}

/****************************************************************
*  Test the evolutionary algorithm
*/
void ea_test(t_swarm *x, t_symbol *sym, long argc, t_atom *argv)
{
  t_swr_ind ind1, ind2;
  t_agent *agent;
  t_swr_float *pf;

  if (atom_getsym(argv) == gensym("select")) {

    if (atom_getsym(argv + 1) == gensym("rank")) {
      ea_select_rank(x, &ind1, NULL, (t_swr_ind)atom_getlong(argv + 2));
      post("test_ea select rank 1:  %i", ind1);
      ea_select_rank(x, &ind1, &ind2, (t_swr_ind)atom_getlong(argv + 2));
      post("test_ea select rank 2:  %i  %i", ind1, ind2);

    } else if (atom_getsym(argv + 1) == gensym("roulette")) {
      ea_select_roulette(x, &ind1, NULL, 0);
      post("test_ea select roulette 1:  %i", ind1);
      ea_select_roulette(x, &ind1, &ind2, 0);
      post("test_ea select roulette 2:  %i  %i", ind1, ind2);

    } else if (atom_getsym(argv + 1) == gensym("tournament")) {
      ea_select_tournament(x, &ind1, NULL, 0);
      post("test_ea select tournament 1:  %i", ind1);
      ea_select_tournament(x, &ind1, &ind2, 0);
      post("test_ea select tournament 2:  %i  %i", ind1, ind2);
    }
  }

  if (atom_getsym(argv) == gensym("cross")) {
    ind1 = CLAMP((t_swr_ind)atom_getlong(argv + 2), 0, GENE_REC_CNT - 1);
    ind2 = CLAMP((t_swr_ind)atom_getlong(argv + 3), 0, GENE_REC_CNT - 1);
    agent = x->agent_arr + CLAMP(atom_getlong(argv + 4), 0, x->agent_max - 1);
    
    if (atom_getsym(argv + 1) == gensym("single")) {
      ea_cross_single(x, agent, ind1, ind2);
      post("test_ea cross single");

    } else if (atom_getsym(argv + 1) == gensym("double")) {
      ea_cross_double(x, agent, ind1, ind2);
      post("test_ea cross double");

    } else if (atom_getsym(argv + 1) == gensym("uniform")) {
      ea_cross_uniform(x, agent, ind1, ind2);
      post("test_ea cross uniform");
    }

    pf = x->gene_bank_arr[x->gene_bank_sort[ind1]].gene_nval;
    post("src1:  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f",
      pf[0], pf[1], pf[2], pf[3], pf[4], pf[5], pf[6], pf[7], pf[8], pf[9], pf[10]);
    pf = x->gene_bank_arr[x->gene_bank_sort[ind2]].gene_nval;
    post("src2:  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f",
      pf[0], pf[1], pf[2], pf[3], pf[4], pf[5], pf[6], pf[7], pf[8], pf[9], pf[10]);
    pf = agent->gene_nval;
    post("nval:  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f",
      pf[0], pf[1], pf[2], pf[3], pf[4], pf[5], pf[6], pf[7], pf[8], pf[9], pf[10]);
    pf = agent->gene_val;
    post("val:   %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f",
      pf[0], pf[1], pf[2], pf[3], pf[4], pf[5], pf[6], pf[7], pf[8], pf[9], pf[10]);
  }

  if (atom_getsym(argv) == gensym("mutate")) {
    agent = x->agent_arr + CLAMP(atom_getlong(argv + 1), 0, x->agent_max - 1);

    post("test_ea mutate");
    pf = agent->gene_nval;
    post("pre nval:  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f",
      pf[0], pf[1], pf[2], pf[3], pf[4], pf[5], pf[6], pf[7], pf[8], pf[9], pf[10]);
    pf = agent->gene_val;
    post("pre val:   %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f",
      pf[0], pf[1], pf[2], pf[3], pf[4], pf[5], pf[6], pf[7], pf[8], pf[9], pf[10]);

    ea_mutate(x, agent);

    pf = agent->gene_nval;
    post("post nval:  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f",
      pf[0], pf[1], pf[2], pf[3], pf[4], pf[5], pf[6], pf[7], pf[8], pf[9], pf[10]);
    pf = agent->gene_val;
    post("post val:   %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f",
      pf[0], pf[1], pf[2], pf[3], pf[4], pf[5], pf[6], pf[7], pf[8], pf[9], pf[10]);
  }
  return;
}
