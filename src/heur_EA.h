
#ifndef __HEAR_EA_H__
#define __HEAR_EA_H__

#include "objscip/objscip.h"

namespace EA
{

/** C++ evolutionary algorithm heuristic for MIPs */
class HeurEA : public  scip::ObjHeur
{

   /* primal heuristic data */
   int submip_node_limit; // maximum number of nodes that are explored in subproblem
   int solution_pool_pop; // size of solution pool
   SCIP_Real fixing_fraction; // fraction of variables that are fixed 
   SCIP_Real fixing_increment; // initial fixing increment 
   SCIP_Real fixing_increment_reduction; // percentage reduction of the fixing increment
   SCIP_Real min_fixing_increment; // minimum fixing increment 
   int number_of_mutations; // number of mutation candidates generated
   int number_of_combinations; // number of combination steps
   int seed; // seed for ranom number generation
   SCIP_RANDNUMGEN* randnumgen;

public:
      
   /** default constructor */
   HeurEA(
      SCIP* scip,              
      int submip_node_limit = 500,
      int solution_pool_pop = 40,
      SCIP_Real fixing_fraction = 0.5, 
      SCIP_Real fixing_increment = 0.2, 
      SCIP_Real fixing_increment_reduction = 0.25, 
      SCIP_Real min_fixing_increment = 0.01, 
      int number_of_mutations = 20,
      int number_of_combinations = 40,     
      int seed = 1
      )
      : ObjHeur(scip, "EA", "evolutionary algorithm heuristic for MIPs", 'E', 1210000, 1, 1, -1,
            SCIP_HEURTIMING_BEFORENODE, TRUE),
      submip_node_limit(submip_node_limit),
      solution_pool_pop(solution_pool_pop),
      fixing_fraction(fixing_fraction), 
      fixing_increment(fixing_increment), 
      fixing_increment_reduction(fixing_increment_reduction), 
      min_fixing_increment(min_fixing_increment), 
      number_of_mutations(number_of_mutations),
      number_of_combinations(number_of_combinations),
      seed(seed)
   {
   
   }


   /** destructor */
   virtual ~HeurEA()
   {

   }

   /** destructor of primal heuristic to free user data (called when SCIP is exiting)
    *
    *  @see SCIP_DECL_HEURFREE(x) in @ref type_heur.h
    */
   virtual SCIP_DECL_HEURFREE(scip_free);

   /** initialization method of primal heuristic (called after problem was transformed)
    *
    *  @see SCIP_DECL_HEURINIT(x) in @ref type_heur.h
    */
   virtual SCIP_DECL_HEURINIT(scip_init);

   /** deinitialization method of primal heuristic (called before transformed problem is freed)
    *
    *  @see SCIP_DECL_HEUREXIT(x) in @ref type_heur.h
    */
   virtual SCIP_DECL_HEUREXIT(scip_exit);

   /** solving process initialization method of primal heuristic (called when branch and bound process is about to begin)
    *
    *  @see SCIP_DECL_HEURINITSOL(x) in @ref type_heur.h
    */
   virtual SCIP_DECL_HEURINITSOL(scip_initsol);

   /** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
    *
    *  @see SCIP_DECL_HEUREXITSOL(x) in @ref type_heur.h
    */
   virtual SCIP_DECL_HEUREXITSOL(scip_exitsol);

   /** execution method of primal heuristic
    *
    *  @see SCIP_DECL_HEUREXEC(x) in @ref type_heur.h
    */
   virtual SCIP_DECL_HEUREXEC(scip_exec);

   /** clone method which will be used to copy a objective plugin */
   virtual SCIP_DECL_HEURCLONE(ObjCloneable* clone);

   /** returns whether the objective plugin is copyable */
   virtual SCIP_DECL_HEURISCLONEABLE(iscloneable)
   {
      return TRUE;
   }

};
} /* namespace EA */

#endif
