#include <assert.h>
#include "heur_EA.h"
#include <iostream>

using namespace EA;
#define SCIP_DEBUG
/*
 * Callback methods of primal heuristic
 */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
SCIP_DECL_HEURFREE(HeurEA::scip_free)
{  
   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
SCIP_DECL_HEURINIT(HeurEA::scip_init)
{ 
   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen,
         SCIPinitializeRandomSeed(scip, seed), TRUE) );

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
SCIP_DECL_HEUREXIT(HeurEA::scip_exit)
{  

   SCIPfreeRandom(scip, &randnumgen);
   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
SCIP_DECL_HEURINITSOL(HeurEA::scip_initsol)
{  
   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
SCIP_DECL_HEUREXITSOL(HeurEA::scip_exitsol)
{  
   return SCIP_OKAY;
}

static
SCIP_RETCODE performMutation(SCIP* scip, SCIP_HEUR* heur,
      int submip_node_limit,
      int solution_pool_pop,
      SCIP_Real* fixing_fraction,
      SCIP_Real fixing_increment, 
      SCIP_RANDNUMGEN* randnumgen
)
{

   SCIP* submip; // the subproblem created by mutation
   SCIP_VAR** fixedvars; // contains fixed variables 
   SCIP_Real* fixedvals; // contains values of fixed variables 
   SCIP_VAR** vars;// stores original variables
   int number_bin_vars;
   int number_int_vars;
   int number_vars;

   SCIP_Bool success;
   assert(scip != NULL);

   int num_sols = SCIPgetNSols(scip);

   // need at least 1 feasible solution to perform a mutation
   if( num_sols == 0 )
      return SCIP_OKAY;

   SCIP_SOL** sols;
   SCIP_SOL* sol;

   sols = SCIPgetSols(scip);

   assert(sols != NULL);

   int seed_sol = SCIPrandomGetInt(randnumgen, 0, 
      ((num_sols> solution_pool_pop)? solution_pool_pop: num_sols)-1);

   sol = sols[seed_sol];

   assert(sol != NULL);

   // get solution data
   SCIP_CALL(SCIPgetVarsData(scip, &vars, &number_vars, &number_bin_vars, &number_int_vars, NULL, NULL));
   number_int_vars += number_bin_vars;

   // need integer vars to perform a mutation
   if (number_int_vars == 0)
      return SCIP_OKAY;

   // calculate the number of fixed vars
   int number_fixed_vars = (int)((*fixing_fraction) * number_int_vars);


   // can not mutate if none/all of the variables are fixed
   if (number_int_vars == 0 || number_fixed_vars >= number_int_vars ){
      (*fixing_fraction) -= fixing_increment;
      return SCIP_OKAY;
   }

   // memory stuff
   SCIP_CALL( SCIPallocBufferArray(scip, &fixedvars, number_fixed_vars));
   SCIP_CALL( SCIPallocBufferArray(scip, &fixedvals, number_fixed_vars));
   BMScopyMemoryArray(fixedvars, vars, number_fixed_vars);


   assert(fixedvars != NULL);
   assert(fixedvals != NULL);

   // permute variable array to simulate picking fixed variables at random
   SCIPrandomPermuteArray(randnumgen, (void **)fixedvars, 0, number_fixed_vars);



   // fix first number_fixed_vars of array   
   for(int i = 0; i < number_fixed_vars; i++ )
   {

      SCIP_Real solval;
      SCIP_Real lb;
      SCIP_Real ub;

      solval = SCIPgetSolVal(scip, sol, fixedvars[i]);
      lb = SCIPvarGetLbGlobal(fixedvars[i]);
      ub = SCIPvarGetUbGlobal(fixedvars[i]);

      assert(SCIPisLE(scip, lb, ub));

      // this case was checked in the example
      if( SCIPisLT(scip, solval, lb) )
         solval = lb;
      else if( SCIPisGT(scip, solval, ub) )
         solval = ub;

      // do not fix to infinity
      if( SCIPisInfinity(scip, REALABS(solval)) )
      {
         SCIPfreeBufferArray(scip, &fixedvals);
         SCIPfreeBufferArray(scip, &fixedvars);
         return SCIP_OKAY;
      }

      fixedvals[i] = solval;
   }


   // creating the subproblem 
   SCIP_CALL(SCIPcreate(&submip));


   assert(submip != NULL);
   assert(fixedvars != NULL);
   assert(fixedvals != NULL);
 
   // creating subMIP code
   // why does it need a hashmap?
   SCIP_VAR** submip_vars; //sub MIP vars
   SCIP_HASHMAP* varmap; // mapping of original vars to subMIP vars

   SCIP_CALL(SCIPallocBufferArray(scip, &submip_vars, number_vars));
   SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(submip), number_vars) );

   SCIP_CALL(SCIPcopyLargeNeighborhoodSearch(scip, submip, varmap, "mutation", 
         fixedvars, fixedvals, number_fixed_vars,
         FALSE, TRUE,
         &success, NULL) );

   if (!success){
      SCIPfreeBufferArray(scip, &submip_vars);
      SCIP_CALL(SCIPfree(&submip));
      SCIPfreeBufferArray(scip, &fixedvals);
      SCIPfreeBufferArray(scip, &fixedvars);         
      SCIPhashmapFree(&varmap);

      return SCIP_OKAY;
   }

   for( int i = 0; i < number_vars; i++ )
      submip_vars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, vars[i]);

   SCIPhashmapFree(&varmap);

   // make subMIP less spammy
   SCIP_CALL( SCIPsetIntParam(submip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetBoolParam(submip, "timing/statistictiming", FALSE) );

   // subMIP limits
   SCIP_CALL(SCIPcopyLimits(scip, submip));
   SCIP_CALL(SCIPsetLongintParam(submip, "limits/nodes", submip_node_limit));
   SCIP_CALL(SCIPsetSubscipsOff(submip, TRUE) );
   SCIP_CALL( SCIPsetBoolParam(submip, "lp/checkdualfeas", FALSE) );
   /* solve the subproblem */
   SCIP_CALL_ABORT(SCIPsolve(submip));

   // stats from sub mip
   SCIP_CALL(SCIPmergeVariableStatistics(submip, scip, submip_vars, vars, number_vars));
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(submip, NULL)));

   // see if there is a solution to add
   SCIP_SOL* oldbestsol = SCIPgetBestSol(scip);

   SCIP_SOL*  newsol; // solution to be created for the original problem
   SCIP_SOL* subsol = SCIPgetBestSol(submip); // sub mip solution  

   
   // see if there is a solution to add
   assert(heur!=NULL);
   
   // see if there is a solution to add
   if (SCIPgetNSols(submip) == 0){
      success = FALSE;
   } else {
      // try to add new solution to scip
      SCIP_CALL(SCIPtranslateSubSol(scip, submip, subsol, heur, submip_vars, &newsol ));
      SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, &success) );

   }

   if (!success || oldbestsol == SCIPgetBestSol(scip)) {
      (*fixing_fraction) += fixing_increment;
   } else if (oldbestsol == newsol){
         (*fixing_fraction) -= fixing_increment;
   }
   // memory stuff
   SCIPfreeBufferArray(scip, &submip_vars);
   SCIP_CALL(SCIPfree(&submip));
   SCIPfreeBufferArray(scip, &fixedvals);
   SCIPfreeBufferArray(scip, &fixedvars);

   return SCIP_OKAY;
}

static
SCIP_RETCODE performCombination(SCIP* scip, SCIP_HEUR* heur,
      int submip_node_limit,
      int solution_pool_pop,
      SCIP_RANDNUMGEN* randnumgen,
      bool combine_all
)
{  
   SCIP* submip; // the subproblem created by mutation
   SCIP_VAR** fixedvars; // contains fixed variables 
   SCIP_Real* fixedvals; // contains values of fixed variables 
   SCIP_VAR** vars;// stores original variables
   int number_bin_vars;
   int number_int_vars;
   int number_vars;

   SCIP_Bool success;

   assert(scip != NULL);

   int num_sols = SCIPgetNSols(scip);

   // need at least 2 feasible solution to perform a combination
   if( num_sols <= 1 )
      return SCIP_OKAY;

   SCIP_SOL** sols;
   SCIP_SOL* sol;

   sols = SCIPgetSols(scip);
   assert(sols != NULL);

   // get var data
   SCIP_CALL(SCIPgetVarsData(scip, &vars, &number_vars, &number_bin_vars, &number_int_vars, NULL, NULL));
   number_int_vars += number_bin_vars;

   // can not combine
   if (number_int_vars == 0){
      return SCIP_OKAY;
   }

   int* sols_comb_ids;
   int number_sols_used = 2;

   if (combine_all){
      number_sols_used = num_sols;
   } 
   assert(number_sols_used <= num_sols);

   SCIP_CALL(SCIPallocBufferArray(scip, &sols_comb_ids, number_sols_used));
   
   if (combine_all){
      for (int i = 0; i < number_sols_used; i++){
         sols_comb_ids[i]=i;
      }
   } else {
      sols_comb_ids[0] = SCIPrandomGetInt(randnumgen, 1, 
         ((num_sols> solution_pool_pop) ? solution_pool_pop: num_sols)-1);
      sols_comb_ids[1] = SCIPrandomGetInt(randnumgen, 0, sols_comb_ids[0]);
   }

   // do not continue if all solutions have been found by the same heuristic at the same node
   SCIP_HEUR* solheur;
   SCIP_Longint solnodenum;
   SCIP_Bool allsame;

   solheur = SCIPsolGetHeur(sols[0]);
   solnodenum = SCIPsolGetNodenum(sols[0]);
   allsame = TRUE;

   for (int i = 1; i < number_sols_used; i++){
      if( SCIPsolGetHeur(sols[i]) != solheur || SCIPsolGetNodenum(sols[i]) != solnodenum )
         allsame = FALSE;
   }

   if (allsame){
      SCIPfreeBufferArray(scip, &sols_comb_ids);
      return SCIP_OKAY;
   }

   // memory stuff
   SCIP_CALL( SCIPallocBufferArray(scip, &fixedvars, number_int_vars));
   SCIP_CALL( SCIPallocBufferArray(scip, &fixedvals, number_int_vars));

   int number_fixed_vars = 0;
   // fix int vars with common vals   
   for(int i = 0; i < number_int_vars; i++ )
   {
      SCIP_Real solval;
      SCIP_Bool fixable;

      bool fixed = TRUE;

      solval = SCIPgetSolVal(scip, sols[sols_comb_ids[0]], vars[i]);

      // var has same val in all sols
      for(int j = 1; j < number_sols_used; j++ )
      {
         SCIP_Real sol_var_val;
         sol_var_val= SCIPgetSolVal(scip, sols[sols_comb_ids[j]], vars[i]);
         if( !SCIPisEQ(scip,solval, sol_var_val) )
         {
            fixed = FALSE;
            break;
         }
      }

      // bounds
      fixed = fixed && SCIPvarGetLbGlobal(vars[i]) <= solval && solval <= SCIPvarGetUbGlobal(vars[i]);


      if(fixed)
      {
         fixedvars[number_fixed_vars] = vars[i];
         fixedvals[number_fixed_vars] = solval;
         number_fixed_vars++;
      }
   }

   // creating the subproblem 
   SCIP_CALL(SCIPcreate(&submip));
   assert(submip != NULL);

   // creating subMIP code
   // why does it need a hashmap?

   SCIP_VAR** submip_vars; //sub MIP vars
   SCIP_HASHMAP* varmap; // mapping of original vars to subMIP vars

   SCIP_CALL(SCIPallocBufferArray(scip, &submip_vars, number_vars));
   SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(submip), number_vars) );
   SCIP_CALL(SCIPcopyLargeNeighborhoodSearch(scip, submip, varmap, "combination", 
         fixedvars, fixedvals, number_fixed_vars,
         FALSE, TRUE,
         &success, NULL) );

   for( int i = 0; i < number_vars; i++ )
      submip_vars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, vars[i]);

   SCIPhashmapFree(&varmap);

   // make subMIP less spammy
   SCIP_CALL(SCIPsetIntParam(submip, "display/verblevel", 0) );
   SCIP_CALL(SCIPsetBoolParam(submip, "timing/statistictiming", FALSE) );

   // subMIP limits
   SCIP_CALL(SCIPcopyLimits(scip, submip));
   SCIP_CALL(SCIPsetLongintParam(submip, "limits/nodes", submip_node_limit));
   SCIP_CALL(SCIPsetSubscipsOff(submip, TRUE) );
   SCIP_CALL( SCIPsetBoolParam(submip, "lp/checkdualfeas", FALSE) );

   //
   SCIP_CALL_ABORT(SCIPsolve(submip));

   // stats from sub mip
   SCIP_CALL(SCIPmergeVariableStatistics(submip, scip, submip_vars, vars, number_vars));
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(submip, NULL)));

   // see if there is a solution to add
   assert(heur!=NULL);
   // see if there is a solution to add
   if (SCIPgetNSols(submip) > 0){
      // try to add new solution to scip
      SCIP_CALL( SCIPtranslateSubSols(scip, submip, heur, submip_vars, &success, NULL) );
   }

   // memory stuff
   SCIPfreeBufferArray(scip, &submip_vars);
   SCIP_CALL(SCIPfree(&submip));
   SCIPfreeBufferArray(scip, &fixedvals);
   SCIPfreeBufferArray(scip, &fixedvars);
   SCIPfreeBufferArray(scip, &sols_comb_ids);

   return SCIP_OKAY;
}


/** execution method of primal heuristic  */
SCIP_DECL_HEUREXEC(HeurEA::scip_exec)
{

   SCIP* subsmip; // the subproblem created by mutation
   SCIP_VAR** fixedvars; // contains fixed variables 
   SCIP_Real* fixedvals; // contains values of fixed variables 

   SCIP_Bool success;

   assert(scip != NULL);
   assert(result != NULL);
   assert(heur!=NULL);


   *result = SCIP_DELAYED;
   if(SCIPgetNNodes(scip) < 20000 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   int old_sol_num = SCIPgetNSols(scip);
   // need at least 1 feasible solution to perform a mutation
   if( old_sol_num <= 0 )
      return SCIP_OKAY;



   int number_bin_vars;
   int number_int_vars;

   SCIP_CALL(SCIPgetVarsData(scip, NULL, NULL, &number_bin_vars, &number_int_vars, NULL, NULL));
   number_int_vars += number_bin_vars;

   // to apply heuristic
   if (number_int_vars == 0)
      return SCIP_OKAY;

   SCIP_CALL( SCIPcheckCopyLimits(scip, &success) );
   if( !success )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) )
     return SCIP_OKAY;

   // SCIPprintMemoryDiagnostic(scip);

   *result = SCIP_DIDNOTFIND;

   // mutate 
   for (int i = 0; i < number_of_mutations; i++){
      /* check whether there is enough time and memory left */
      performMutation(scip, heur, submip_node_limit, solution_pool_pop, 
         &fixing_fraction, fixing_increment, randnumgen);
   }
   if (fixing_increment > 0.01) {
      fixing_increment = (1-fixing_increment_reduction)*fixing_increment;
      if (fixing_increment < 0.01) {   
         fixing_increment = 0.01;
      }
   }

   // combine
   for (int i = 1; i < number_of_combinations; i++){
      performCombination(scip, heur, submip_node_limit, solution_pool_pop, 
         randnumgen, FALSE);
   }

   performCombination(scip, heur,submip_node_limit, solution_pool_pop, 
         randnumgen, TRUE);
   
   if(old_sol_num < SCIPgetNSols(scip)){
         *result = SCIP_FOUNDSOL;
   }
   return SCIP_OKAY;

}

/** clone method which will be used to copy a objective plugin */
SCIP_DECL_HEURCLONE(scip::ObjCloneable* HeurEA::clone) 
{
   return new HeurEA(scip);
}