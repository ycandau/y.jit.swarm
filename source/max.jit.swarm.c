/**
  jit.swarm - a Max MSP object for swarm sonification
  Yves Candau - ycandau@gmail.com

  Max wrapper object definition
*/

// ========  PREPROCESSOR DEFINITIONS  ========

#define C74_NO_DEPRECATION    // disable deprecation warnings (C4996)

// ========  HEADER FILES  ========

#include "jit.common.h"
#include "max.jit.mop.h"
#include "jit.swarm.h"

// ========  MAX WRAPPER STRUCTURE  ========

typedef struct _max_jit_swarm
{
	t_object  ob;
	void     *obex;   // extra information and resources to wrap the Jitter class

} t_max_jit_swarm;

// ========  FUNCTION PROTOTYPES  ========

t_jit_err jit_swarm_init(void);

void *max_jit_swarm_new(t_symbol *s, long argc, t_atom *argv);
void  max_jit_swarm_free(t_max_jit_swarm *x);

void  max_jit_swarm_outputmatrix(t_max_jit_swarm *x);

// ========  CLASS POINTER  ========

t_messlist *max_jit_swarm_class;

// ========  INITIALIZATION  ========

void ext_main(void *r)
{
	void *p,*q;

  // Initialize the Jitter class, from jit.swarm.c
	jit_swarm_init();

  // Create the Max class
	setup(&max_jit_swarm_class,         // Max class pointer
    (method)max_jit_swarm_new,        // Max class constructor
    (method)max_jit_swarm_free,       // Max class destructor
    (short)sizeof(t_max_jit_swarm),   // Max class structure size
		0L, A_GIMME, 0);                  

  // Specify a byte offset to keep additional information
	p = max_jit_classex_setup(calcoffset(t_max_jit_swarm, obex));

  // Look up the Jitter class in the class registry
	q = jit_class_findbyname(gensym("jit_swarm"));

  // Add default methods and attributes to the MOP Max wrapper class
	max_jit_classex_mop_wrap(p, q,
    //MAX_JIT_MOP_FLAGS_OWN_BANG |        // override default bang method
    MAX_JIT_MOP_FLAGS_OWN_OUTPUTMATRIX |  // override default outputmatrix method
    MAX_JIT_MOP_FLAGS_OWN_JIT_MATRIX |
    MAX_JIT_MOP_FLAGS_OWN_PLANECOUNT |    // remove default MOP attributes
    MAX_JIT_MOP_FLAGS_OWN_TYPE |
    MAX_JIT_MOP_FLAGS_OWN_DIM |
    MAX_JIT_MOP_FLAGS_OWN_NAME |
    MAX_JIT_MOP_FLAGS_OWN_ADAPT |
    MAX_JIT_MOP_FLAGS_OWN_OUTPUTMODE);

  // Wrap the Jitter class with the standard methods for Jitter objects
	max_jit_classex_standard_wrap(p, q, 0);
	
  //max_addmethod_usurp_low((method)max_jit_swarm_outputmatrix, "bang");
  max_addmethod_usurp_low((method)max_jit_swarm_outputmatrix, "outputmatrix");

  // Add an inlet/outlet assistance method
	addmess((method)max_jit_mop_assist, "assist", A_CANT, 0);
}

// ========  MAX CLASS CONSTRUCTOR  ========

void *max_jit_swarm_new(t_symbol *s, long argc, t_atom *argv)
{
	t_max_jit_swarm *x;
	void *jitter_obj;

  // Create a new Max wrapper
	if (x = (t_max_jit_swarm *)max_jit_obex_new(
    max_jit_swarm_class,                        // class pointer
    gensym("jit_swarm"))) {                     // class name
		
    // Create a new Jitter object
    if (jitter_obj = jit_object_new(gensym("jit_swarm"))) {
			max_jit_mop_setup_simple(x, jitter_obj, argc, argv);    // setup standard MOP Max wrapper
			max_jit_attr_args(x, (short)argc, argv);                // process attribute arguments

    // Allocation error
		} else {
			jit_object_error((t_object *)x, "jit.swarm:  Object allocation failed.");
			freeobject((t_object *) x);
			x = NULL;
		}
	}

  // Store the Max wrapper pointer in the Jitter object
  t_swarm *jo = (t_swarm *)jitter_obj;
  jo->max_wrap = x;

	return (x);
}

// ========  MAX CLASS DESTRUCTOR  ========

void max_jit_swarm_free(t_max_jit_swarm *x)
{
  // Free the MOP wrapper
  max_jit_mop_free(x);

  // Get the Jitter object from the Max wrapper and free it
  jit_object_free(max_jit_obex_jitob_get(x));

  // Free the Max wrapper object
  max_jit_obex_free(x);
}

/*void max_jit_swarm_bang(t_max_jit_swarm *x)
{
  jit_object_method(max_jit_obex_jitob_get(x), 

}*/

// ========  OUPUTMATRIX  ========

void max_jit_swarm_outputmatrix(t_max_jit_swarm *x)
{
  // Get the outputmode
  long outputmode = max_jit_mop_getoutputmode(x);

  // Get the MOP wrapper
  void *mop = max_jit_obex_adornment_get(x, _jit_sym_jit_mop);

  t_jit_err err;

  // Always output unless output mode is none
  if (outputmode && mop) {
    if (outputmode == 1) {
      if (err = (t_jit_err)jit_object_method(
        max_jit_obex_jitob_get(x),
        _jit_sym_matrix_calc,
        jit_object_method(mop, _jit_sym_getinputlist),
        jit_object_method(mop, _jit_sym_getoutputlist))) {
        jit_error_code(x, err);

      } else {
        max_jit_mop_outputmatrix(x);
      }
    } else {
      max_jit_mop_outputmatrix(x);
    }
  }
}