#ifndef YC_JIT_SWARM_H_
#define YC_JIT_SWARM_H_

//******************************************************************************
//  Header files
//
#include "jit.common.h"   // General Jitter header file
#include "max.jit.mop.h"
#include <time.h>         // To initialize random function

//******************************************************************************
//  Preprocessor
//
#define PI 3.141593f

#define AGENT_MAX_DEF 20
#define GENE_REC_CNT 100

#define SWR_FLOAT_SYM _jit_sym_float32

// Tracing macro for debugging purposes
#define _TRACE false
#define TRACE(str, ...) do { if (_TRACE) object_post(NULL, "TRACE:  " str, __VA_ARGS__); } while (0)

// Tracing macro for repeatedly called functions
#define _TRACE_F false
#define TRACE_F(str, ...) do { if (_TRACE_F) object_post(NULL, "TRACE:  " str, __VA_ARGS__); } while (0)

#define COHES_MUL_BASE  200.0f
#define ALIGN_MUL_BASE  20.0f
#define ATTRA_MUL_BASE  10.0f
#define SEPAR_MUL_BASE  6.0f
#define BROWN_MUL_BASE  20.0f
#define BOUND_MUL_BASE  50.0f

#define ATTRA_THR_BASE  150.0f
#define SEPAR_THR_BASE  30.0f
#define BOUND_THR_BASE  150.0f
#define BROWN_PROB_BASE 0.3f

#define SIZE_MUL_BASE   0.015f

#define BIN_CNT_X    16
#define BIN_CNT_Y    9
#define BIN_CNT_VEL  16
#define BIN_CNT_CURV 16

#define CNTD_RAND  0xFFFF
#define NO_EXCLUDE 0xFFFF

//******************************************************************************
//  Macros for vector calculations.
//
#define V2D_V_e_S(V, S)       { V[0] = S;     V[1] = S;     }
#define V2D_V_e_SS(V, S0, S1) { V[0] = S0;    V[1] = S1;    }
#define V2D_V_e_V(V, V1)      { V[0] = V1[0]; V[1] = V1[1]; }

#define V2D_V_pe_S(V, S)       { V[0] += S;       V[1] += S;       }
#define V2D_V_pe_SS(V, S0, S1) { V[0] += S0;      V[1] += S1;      }
#define V2D_V_pe_V(V, V1)      { V[0] += V1[0];   V[1] += V1[1];   }
#define V2D_V_pe_S_V(V, S, V1) { V[0] += S*V1[0]; V[1] += S*V1[1]; }

#define V2D_V_te_S(V, S)       { V[0] *= S;     V[1] *= S;     }
#define V2D_V_te_SS(V, S0, S1) { V[0] *= S0;    V[1] *= S1;    }
#define V2D_V_te_V(V, V1)      { V[0] *= V1[0]; V[1] *= V1[1]; }

#define V2D_V_e_S_V_V(V, S, V1, V2)  { V[0]  = S * (V1[0] - V2[0]); V[1]  = S * (V1[1] - V2[1]); }
#define V2D_V_pe_S_V_V(V, S, V1, V2) { V[0] += S * (V1[0] - V2[0]); V[1] += S * (V1[1] - V2[1]); }

#define V2D_N2(V)     (V[0] * V[0] + V[1] * V[1])
#define V2D_N(V) sqrtf(V[0] * V[0] + V[1] * V[1])

#define V2D_N2_V_V(V1, V2)     ((V1[0]-V2[0]) * (V1[0]-V2[0]) + (V1[1]-V2[1]) * (V1[1]-V2[1]))
#define V2D_N_V_V(V1, V2) sqrtf((V1[0]-V2[0]) * (V1[0]-V2[0]) + (V1[1]-V2[1]) * (V1[1]-V2[1]))

#define V3D_V_e_S(V, S)            { V[0] = S;     V[1] = S;     V[2] = S;     }
#define V3D_V_e_SSS(V, S0, S1, S2) { V[0] = S0;    V[1] = S1;    V[2] = S2;    }
#define V3D_V_e_V(V, V1)           { V[0] = V1[0]; V[1] = V1[1]; V[2] = V1[2]; }

#define V3D_V_pe_S(V, S)            { V[0] += S;       V[1] += S;       V[2] += S;       }
#define V3D_V_pe_SSS(V, S0, S1, S2) { V[0] += S0;      V[1] += S1;      V[2] += S2;      }
#define V3D_V_pe_V(V, V1)           { V[0] += V1[0];   V[1] += V1[1];   V[2] += V1[2];   }
#define V3D_V_pe_S_V(V, S, V1)      { V[0] += S*V1[0]; V[1] += S*V1[1]; V[2] += S*V1[2]; }

#define V3D_V_te_S(V, S)            { V[0] *= S;     V[1] *= S;     V[2] *= S;     }
#define V3D_V_te_SSS(V, S0, S1, S2) { V[0] *= S0;    V[1] *= S1;    V[2] *= S2;    }
#define V3D_V_te_V(V, V1)           { V[0] *= V1[0]; V[1] *= V1[1]; V[2] *= V1[2]; }

#define V3D_V_e_S_V_V(V, S, V1, V2)  { V[0]  = S * (V1[0] - V2[0]); V[1]  = S * (V1[1] - V2[1]); V[2]  = S * (V1[2] - V2[2]); }
#define V3D_V_pe_S_V_V(V, S, V1, V2) { V[0] += S * (V1[0] - V2[0]); V[1] += S * (V1[1] - V2[1]); V[2] += S * (V1[2] - V2[2]); }

#define V3D_N2(V)     (V[0] * V[0] + V[1] * V[1] + V[2] * V[2])
#define V3D_N(V) sqrtf(V[0] * V[0] + V[1] * V[1] + V[2] * V[2])

#define V3D_N2_V_V(V1, V2)     ((V1[0]-V2[0]) * (V1[0]-V2[0]) + (V1[1]-V2[1]) * (V1[1]-V2[1]) + (V1[2]-V2[2]) * (V1[2]-V2[2]))
#define V3D_N_V_V(V1, V2) sqrtf((V1[0]-V2[0]) * (V1[0]-V2[0]) + (V1[1]-V2[1]) * (V1[1]-V2[1]) + (V1[2]-V2[2]) * (V1[2]-V2[2]))

//******************************************************************************
//  Jitter class pointer declaration as extern variable.
//
extern void* _jit_swarm_class;

//******************************************************************************
//  Typedef and structure declarations
//
// Variable types

typedef t_uint16 t_swr_ind;     // for agent indexes and arrays
typedef t_uint16 t_swr_cnt;     // for time counters, 0xFFFF at 30 fps is 36 min
typedef float    t_swr_float;   // all float calculations

// Swarm structures

typedef struct _agent     t_agent;
typedef struct _jit_swarm t_swarm;

// Evolutionary algorithm structures

typedef struct _gene_repr t_gene_repr;
typedef struct _gene_bank  t_gene_bank;

// Function pointers for evolutionary algorithm structures

typedef t_swr_float (*t_ea_scale) (t_swr_float nval, t_gene_repr* gene_repr);
typedef t_swr_float (*t_ea_init)  (t_swarm* x);

typedef void (*t_ea_select) (t_swarm* x, t_swr_ind* ind1, t_swr_ind* ind2, t_swr_ind N);
typedef void (*t_ea_cross)  (t_swarm* x, t_agent* agent, t_swr_ind ind1, t_swr_ind ind2);

//******************************************************************************
//  Enum:  Bit flags to keep track of agent state
//
typedef enum _swr_state {

  ST_NULL     = 0x00,
  ST_ACTIVE   = 0x01,
  ST_DORMANT  = 0x02,
  ST_HATCHING = 0x04,
  ST_FLYING   = 0x08,
  ST_DYING    = 0x10

} e_swr_state;

#define IS_ACTIVE(agent)   ((agent)->state & ST_ACTIVE)
#define IS_INACTIVE(agent) (~(agent)->state & ST_ACTIVE)

//******************************************************************************
//  Enum:  Indexes to access the genes
//
typedef enum _gene_ind {

  GENE_INERT_MUL,
  GENE_COHES_MUL,
  GENE_ALIGN_MUL,
  GENE_ATTRA_MUL,
  GENE_ATTRA_THR,
  GENE_SEPAR_MUL,
  GENE_SEPAR_THR,
  GENE_BOUND_MUL,
  GENE_BOUND_THR,
  GENE_BROWN_MUL,
  GENE_BROWN_PROB,
  GENE_CNT

} e_gene_ind;

//******************************************************************************
//  Function declarations:  jit.swarm.c
//
// Jitter object functions

t_jit_err  jit_swarm_init (void);
t_swarm*   jit_swarm_new  (void);
void       jit_swarm_free (t_swarm* x);

t_jit_err  jit_swarm_matrix_calc    (t_swarm* x, void* inputs, void* outputs);
void       jit_swarm_calculate_ndim (t_swarm* x, long dimcount, long* dim, long planecount,
             t_jit_matrix_info* in_minfo, char* bip, t_jit_matrix_info* out_minfo, char* bop);
void       jit_minfo_post           (t_jit_matrix_info* minfo);

// Swarm functions

void      swarm_init_param (t_swarm* x);
t_jit_err swarm_new_agents  (t_swarm* x, t_swr_ind cur, t_swr_ind max);

void swarm_iter        (t_swarm* x);
void swarm_average     (t_swarm* x);
void swarm_move        (t_swarm* x);
void swarm_track_stat  (t_swarm* x);
void swarm_refresh_act (t_swarm* x);

// Agent functions

t_agent* agent_add (t_swarm* x);
void agent_remove  (t_swarm* x, t_agent* agent);
void agent_init    (t_swarm* x, t_agent* agent);
void agent_seed    (t_swarm* x, t_agent* agent, t_swr_cnt cntd);
void agent_post    (t_swarm* x, t_agent* agent);

void agent_state_change  (t_swarm* x, t_agent* agent);
void agent_iter_hatching (t_swarm* x, t_agent* agent);
void agent_iter_flying   (t_swarm* x, t_agent* agent);
void agent_iter_dying    (t_swarm* x, t_agent* agent);

// Utility functions

t_swr_cnt rand_time(t_swr_float t, t_swr_float dt);

//******************************************************************************
//  Function declarations:  swarm.ea.c
//
void ea_test(t_swarm* x, t_symbol* sym, long argc, t_atom* argv);

void ea_init_gene_repr (t_swarm* x);
void ea_init_gene_rec  (t_swarm* x);

void ea_agent_fitness (t_swarm* x, t_agent* agent);
void ea_agent_store   (t_swarm* x, t_agent* agent);
void ea_agent_stat    (t_swarm* x, t_agent* agent);
void ea_swarm_stat    (t_swarm* x);

void ea_select_rank       (t_swarm* x, t_swr_ind* ind1, t_swr_ind* ind2, t_swr_ind out_of);
void ea_select_roulette   (t_swarm* x, t_swr_ind* ind1, t_swr_ind* ind2, t_swr_ind N);
void ea_select_tournament (t_swarm* x, t_swr_ind* ind1, t_swr_ind* ind2, t_swr_ind N);

t_swr_ind ea_select_roulette_util (t_swarm* x, t_swr_ind exclude);

void ea_cross_single  (t_swarm* x, t_agent* agent, t_swr_ind ind1, t_swr_ind ind2);
void ea_cross_double  (t_swarm* x, t_agent* agent, t_swr_ind ind1, t_swr_ind ind2);
void ea_cross_uniform (t_swarm* x, t_agent* agent, t_swr_ind ind1, t_swr_ind ind2);

void ea_mutate (t_swarm* x, t_agent* agent);

void ea_init_gene_nval (t_swr_float* nval, t_swarm* x);
void ea_calc_gene_val  (t_swr_float* nval, t_swr_float* val, t_swarm* x);
void ea_copy_gene_nval (t_swr_float* nval1, t_swr_float* nval2, t_swarm* x);

void ea_gene_repr_set (t_swarm* x, int g, t_swr_float min, t_swr_float mid, t_swr_float max, t_ea_scale scale, t_bool var, t_symbol* sym);
void ea_post (t_swarm* x);

// Scaling functions:  to be assigned to the function pointer in the structure

t_swr_float ea_fct_scale_plin (t_swr_float nval, t_gene_repr* gene_repr);
t_swr_float ea_fct_scale_exp  (t_swr_float nval, t_gene_repr* gene_repr);

// Initialization functions:  to be assigned to the function pointer in the structure

t_swr_float ea_fct_init_mid (t_swarm* x);
t_swr_float ea_fct_init_lin (t_swarm* x);

//******************************************************************************
//  Function declarations:  swarm.interface.c
//
// Message functions for the object

void swarm_reset   (t_swarm* x);
void swarm_add     (t_swarm* x, t_symbol* sym, long argc, t_atom* argv);
void swarm_remove  (t_swarm* x, t_symbol* sym, long argc, t_atom* argv);
void swarm_post    (t_swarm* x);
void swarm_ea_post (t_swarm* x);

// Object methods and attributes

void swarm_init_method (void* swarm_class);
void swarm_init_attr   (void* swarm_class);
void swarm_init_attr_util(void* swarm_class, t_jit_object* attr,
  const char* name, const char* label, const char* order, float min_, float max_);

// Custom setter functions

t_jit_err swarm_set_agent_max  (t_swarm* x, void* attr, long argc, t_atom* argv);
t_jit_err swarm_set_size       (t_swarm* x, void* attr, long argc, t_atom* argv);
t_jit_err swarm_set_select_fct (t_swarm* x, void* attr, long argc, t_atom* argv);
t_jit_err swarm_set_cross_fct  (t_swarm* x, void* attr, long argc, t_atom* argv);

t_jit_err swarm_set_inert_mul  (t_swarm* x, void* attr, long argc, t_atom* argv);
t_jit_err swarm_set_cohes_mul  (t_swarm* x, void* attr, long argc, t_atom* argv);
t_jit_err swarm_set_align_mul  (t_swarm* x, void* attr, long argc, t_atom* argv);
t_jit_err swarm_set_attra_mul  (t_swarm* x, void* attr, long argc, t_atom* argv);
t_jit_err swarm_set_attra_thr  (t_swarm* x, void* attr, long argc, t_atom* argv);
t_jit_err swarm_set_separ_mul  (t_swarm* x, void* attr, long argc, t_atom* argv);
t_jit_err swarm_set_separ_thr  (t_swarm* x, void* attr, long argc, t_atom* argv);
t_jit_err swarm_set_bound_mul  (t_swarm* x, void* attr, long argc, t_atom* argv);
t_jit_err swarm_set_bound_thr  (t_swarm* x, void* attr, long argc, t_atom* argv);
t_jit_err swarm_set_brown_mul  (t_swarm* x, void* attr, long argc, t_atom* argv);
t_jit_err swarm_set_brown_prob (t_swarm* x, void* attr, long argc, t_atom* argv);
t_jit_err swarm_set_gene_val   (t_swarm* x, e_gene_ind gene, t_swr_float* mod, long argc, t_atom* argv);

//******************************************************************************
//  Agent structure.
//
//  Used for now in a static array in the main Jitter structure.
//
typedef struct _agent {

  t_uint32 state;

  t_swr_cnt age;
  t_swr_cnt cntd;

  t_swr_float pos[3];     // position
  t_swr_float vel[3];     // velocity
  t_swr_float nvel;       // velocity norm
  t_swr_float size;       // size
  t_swr_float col[4];     // colour

  t_swr_float vel_ref[3]; // velocity reference
  t_swr_float size_ref;   // size reference
  t_swr_float col_ref[4]; // colour reference

  // Genome
  t_swr_float gene_nval[GENE_CNT];  // normalized genome
  t_swr_float gene_val[GENE_CNT];   // scaled genome

  // Bin counters for fitness
  t_swr_cnt bin_xy[BIN_CNT_X * BIN_CNT_Y + 1];
  t_swr_cnt bin_vel[BIN_CNT_VEL];
  t_swr_cnt bin_curv[BIN_CNT_CURV];

  // Fitness values, calculated when the agent dies
  t_swr_float fitness;    // fitness
  t_swr_float fit_xy;     // spatial fitness
  t_swr_float fit_vel;    // velocity fitness
  t_swr_float fit_curv;   // curvature fitness

} t_agent;

//******************************************************************************
//  Gene representation structure.
//
//  One static array for the whole population.
//  Describes each gene:
//  - how it is scaled
//  - whether it varies or not
//
typedef struct _gene_repr {

  t_swr_float min;    // the minimum value
  t_swr_float mid;    // the neutral value
  t_swr_float max;    // the maximum value

  t_ea_scale scale_fct;   // the scaling function from the nvalue to the value
  t_bool variable;        // whether the gene is variable or fixed
  t_symbol* sym;          // a symbol for the gene

} t_gene_repr;

//******************************************************************************
//  Gene bank structure.
//
//  One static array for the whole population.
//  Stores the genes of agents after they die.
//
typedef struct _gene_bank {

  t_swr_float gene_nval[GENE_CNT];  // normalized genome
  t_swr_float fitness;              // fitness
  t_swr_float fit_xy;               // spatial fitness
  t_swr_float fit_vel;              // velocity fitness
  t_swr_float fit_curv;             // curvature fitness

} t_gene_bank;

//******************************************************************************
//  Jitter object structure.
//
typedef struct _jit_swarm {

  t_object ob;

  t_agent*    agent_arr;    // array of agents
  t_swr_ind*  agent_act;    // array of active agents
  t_swr_ind   agent_cnt;    // current number of agents

  // Attributes:  Basic - Swarm parameters

  long        agent_max;    // maximum number of agents - leave long for attribute size
  t_swr_float dt;           // rate of temporal change
  t_swr_float vel_max;      // maximum velocity, used for clipping
  t_swr_float size;         // default size for the agents

  // Attributes:  Evolutionary parameters

  t_swr_float fit_xy_mul;   // multiplier for the spatial distribution
  t_swr_float fit_vel_mul;  // multiplier for the velocity distribution

  // Attributes:  Basic - Timing

  t_swr_float dorm_t;       // average dormant time in s
  t_swr_float dorm_dt;      // +/- dormant time (fraction)
  t_swr_float hatch_t;      // average incubation time in s
  t_swr_float hatch_dt;     // +/- incubation time (fraction)
  t_swr_float living_t;     // average life time in s
  t_swr_float living_dt;    // +/- life time (fraction)
  t_swr_float dying_t;      // average incubation time in s
  t_swr_float dying_dt;     // +/- incubation time (fraction)

  // Attributes:  Basic - Agent modifiers

  t_swr_float inert_mul;    // inertia
  t_swr_float cohes_mul;    // cohesion: match the average position
  t_swr_float align_mul;    // alignment: match the average velocity
  t_swr_float attra_mul;    // attraction between agents
  t_swr_float attra_thr;    // attraction threshold
  t_swr_float separ_mul;    // separation between agents
  t_swr_float separ_thr;    // separation threshold
  t_swr_float bound_mul;    // boundary repulsion
  t_swr_float bound_thr;    // boundary threshold
  t_swr_float brown_mul;    // brownian motion
  t_swr_float brown_prob;   // probability of brownian motion

  // Attributes:  Basic - Display

  long length_x;
  long length_y;
  long length_z;

  // Evolutionary algorithm

  t_swr_float init_d;       // max distance from neutral for init: 0-1
  t_ea_init   init_fct;     // function to use for init

  t_swr_ind   rank_cnt;     // number to choose from for rank selection
  t_ea_select select_fct;   // function to use for selection
  long        select_fi;    // index to select the function
  t_ea_cross  cross_fct;    // function to use for crossover
  long        cross_fi;     // index to select the function

  t_swr_float mut_d;        // max distance for mutation
  t_swr_float mut_p;        // probability for mutation

  t_gene_repr gene_repr_arr[GENE_CNT];      // array of gene representations

  t_gene_bank  gene_bank_arr[GENE_REC_CNT];   // array for the gene bank
  t_gene_bank* gene_bank_ptr;                 // pointer to the next element to replace in the bank
  t_uint16     gene_bank_sort[GENE_REC_CNT];  // array of indexes to sort the gene bank

  // Other structure members

  t_swr_float pos_avg[3];   // store the average position
  t_swr_float vel_avg[3];   // store the average velocity

  void* max_wrap;   // pointer to the max wrapper, stored at creation

} t_swarm;

#endif
