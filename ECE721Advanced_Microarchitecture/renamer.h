#include <inttypes.h>
#include <cassert>
#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>
#include <memory>

#include "svp_row_t.h"

using namespace std;

class renamer {
public:
	/////////////////////////////////////////////////////////////////////
	// Put private class variables here.
	/////////////////////////////////////////////////////////////////////

	uint64_t H_al, T_al;
	bool H_al_phase, T_al_phase;

	uint64_t H_fl, T_fl;
	bool H_fl_phase, T_fl_phase;

	uint64_t n_log_regs;
	uint64_t n_phys_regs;
	uint64_t n_branches;
	uint64_t n_active;


	//can move later
	bool st_reg;
	bool st_br;
	bool st_dp;

	/////////////////////////////////////////////////////////////////////
	// Structure 1: Rename Map Table
	// Entry contains: physical register mapping
	/////////////////////////////////////////////////////////////////////
		
	vector<uint64_t> RMT ;

	/////////////////////////////////////////////////////////////////////
	// Structure 2: Architectural Map Table
	// Entry contains: physical register mapping
	/////////////////////////////////////////////////////////////////////

	uint64_t *AMT;
	
	/////////////////////////////////////////////////////////////////////
	// Structure 3: Free List
	//
	// Entry contains: physical register number
	//
	// Notes:
	// * Structure includes head, tail, and their phase bits.
	/////////////////////////////////////////////////////////////////////

	uint64_t *FL;

	/////////////////////////////////////////////////////////////////////
	// Structure 4: Active List
	//
	// Entry contains:
	//
	// ----- Fields related to destination register.
	// 1. destination flag (indicates whether or not the instr. has a
	//    destination register)
	// 2. logical register number of the instruction's destination
	// 3. physical register number of the instruction's destination
	// ----- Fields related to completion status.
	// 4. completed bit
	// ----- Fields for signaling offending instructions.
	// 5. exception bit
	// 6. load violation bit
	//    * Younger load issued before an older conflicting store.
	//      This can happen when speculative memory disambiguation
	//      is enabled.
	// 7. branch misprediction bit
	//    * At present, not ever set by the pipeline. It is simply
	//      available for deferred-recovery Approaches #1 or #2.
	//      Project 1 uses Approach #5, however.
	// 8. value misprediction bit
	//    * At present, not ever set by the pipeline. It is simply
	//      available for deferred-recovery Approaches #1 or #2,
	//      if value prediction is added (e.g., research projects).
	// ----- Fields indicating special instruction types.
	// 9. load flag (indicates whether or not the instr. is a load)
	// 10. store flag (indicates whether or not the instr. is a store)
	// 11. branch flag (indicates whether or not the instr. is a branch)
	// 12. amo flag (whether or not instr. is an atomic memory operation)
	// 13. csr flag (whether or not instr. is a system instruction)
	// ----- Other fields.
	// 14. program counter of the instruction
	//
	// Notes:
	// * Structure includes head, tail, and their phase bits.
	/////////////////////////////////////////////////////////////////////


	typedef struct active_list{

		bool dest_f = 0;
		int log_reg = -1;
		int phys_reg = -1;
		bool completed_bit = 0;
		bool ex = 0;
		bool load_viol = 0;
		bool bmp = 0;
		bool vmp = 0;
		bool lf = 0;
		bool sf = 0;
		bool bf = 0;
		bool amof = 0;
		bool csr = 0;
		uint64_t pc = -1;
		

	}active_list;
	
	active_list *AL;
	

	/////////////////////////////////////////////////////////////////////
	// Structure 5: Physical Register File
	// Entry contains: value
	//
	// Notes:
	// * The value must be of the following type: uint64_t
	//   (#include <inttypes.h>, already at top of this file)
	/////////////////////////////////////////////////////////////////////
	

	uint64_t *PRF;
	

	/////////////////////////////////////////////////////////////////////
	// Structure 6: Physical Register File Ready Bit Array
	// Entry contains: ready bit
	/////////////////////////////////////////////////////////////////////

	uint64_t *PRFr;	

	/////////////////////////////////////////////////////////////////////
	// Structure 7: Global Branch Mask (GBM)
	//
	// The Global Branch Mask (GBM) is a bit vector that keeps track of
	// all unresolved branches. A '1' bit corresponds to an unresolved
	// branch. The "branch ID" of the unresolved branch is its position
	// in the bit vector.
	//
	// The GBM serves two purposes:
	//
	// 1. It provides a means for allocating checkpoints to unresolved
	//    branches. There are as many checkpoints as there are bits in
	//    the GBM. If all bits in the GBM are '1', then there are no
	//    free bits, hence, no free checkpoints. On the other hand, if
	//    not all bits in the GBM are '1', then any of the '0' bits
	//    are free and the corresponding checkpoints are free.
	//    
	// 2. Each in-flight instruction needs to know which unresolved
	//    branches it depends on, i.e., which unresolved branches are
	//    logically before it in program order. This information
	//    makes it possible to squash instructions that are after a
	//    branch, in program order, and not instructions before the
	//    branch. This functionality will be implemented using
	//    branch masks, as was done in the MIPS R10000 processor.
	//    An instruction's initial branch mask is the value of the
	//    the GBM when the instruction is renamed.
	//
	// The simulator requires an efficient implementation of bit vectors,
	// for quick copying and manipulation of bit vectors. Therefore, you
	// must implement the GBM as type "uint64_t".
	// (#include <inttypes.h>, already at top of this file)
	// The "uint64_t" type contains 64 bits, therefore, the simulator
	// cannot support a processor configuration with more than 64
	// unresolved branches. The maximum number of unresolved branches
	// is configurable by the user of the simulator, and can range from
	// 1 to 64.
	/////////////////////////////////////////////////////////////////////
	uint64_t GBM;   // sis-comment: this serves as the bit vector itself ----- needs to be initialized to 0 to indicate that there areno dependencies at the start

	/////////////////////////////////////////////////////////////////////
	// Structure 8: Branch Checkpoints
	//
	// Each branch checkpoint contains the following:
	// 1. Shadow Map Table (checkpointed Rename Map Table)
	// 2. checkpointed Free List head pointer and its phase bit
	// 3. checkpointed GBM
	/////////////////////////////////////////////////////////////////////
	
	 //sis-comment:SMT: a dangling pointer is there - use new and allocate whenever you checkpoint anything 
	 struct BC{
	
		vector<uint64_t> SMT;
		uint64_t H_fl_cp;
		bool H_fl_phase_cp;
		uint64_t GBM_cp = 0;
		uint64_t vpq_tail_cp;
		bool vpq_phase_cp;
		// vector<svp_row_t> svp_table_cp;

		 BC()
        : H_fl_cp(0),
          H_fl_phase_cp(false),
          GBM_cp(0),
          vpq_tail_cp(0),
          vpq_phase_cp(false) 
    {}
		
	};

	 vector<BC> branch_checkpoints;

	/////////////////////////////////////////////////////////////////////
	// Private functions.
	// e.g., a generic function to copy state from one map to another.
	/////////////////////////////////////////////////////////////////////


	void movea2b_cp (uint64_t *b, uint64_t *a) {
	
		for (int i=0; i < n_log_regs; i++) {
		
			b[i] = a[i];

		}
	}

	////////////////////////////////////////
	// Public functions.
	////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////
	// This is the constructor function.
	// When a renamer object is instantiated, the caller indicates:
	// 1. The number of logical registers (e.g., 32).
	// 2. The number of physical registers (e.g., 128).
	// 3. The maximum number of unresolved branches.
	//    Requirement: 1 <= n_branches <= 64.
	// 4. The maximum number of active instructions (Active List size).
	//
	// Tips:
	//
	// Assert the number of physical registers > number logical registers.
	// Assert 1 <= n_branches <= 64.
	// Assert n_active > 0.
	// Then, allocate space for the primary data structures.
	// Then, initialize the data structures based on the knowledge
	// that the pipeline is intially empty (no in-flight instructions yet).
	/////////////////////////////////////////////////////////////////////
	renamer(uint64_t n_log_regs,
		uint64_t n_phys_regs,
		uint64_t n_branches,
		uint64_t n_active){
	
		this->n_log_regs = n_log_regs;
		this->n_phys_regs = n_phys_regs;
		this->n_branches = n_branches;
		this->n_active = n_active;

		assert(n_phys_regs > n_log_regs && "SIS ASSERT: Physical registers provided are less than Logical Registers");
		assert(n_branches >= 1 && n_branches <= 64 && "SIS ASSERT: Invalid Max unresilved branches");
		assert(n_active > 0 && "SIS ASSERT: Need non-zero active list size");
	
		//RMT = new uint64_t[n_log_regs
		RMT.resize(n_log_regs);
		AMT = new uint64_t[n_log_regs];
		FL = new uint64_t[n_phys_regs - n_log_regs];
		AL = new active_list[n_active];
		PRF = new uint64_t[n_phys_regs];
		PRFr = new uint64_t[n_phys_regs];
		//branch_checkpoints = (BC*) malloc(n_branches*sizeof(uint64_t));
		branch_checkpoints.resize(n_branches);
		for (int i = 0 ; i < n_branches ; i++) { 
		//	branch_checkpoints[i].SMT =  new uint64_t[n_log_regs]; //remember to delete all of them once done
			//printf("	Address of checkpoint[%d]'s SMT is : %x\n", i , &branch_checkpoints[i].SMT);
		
			branch_checkpoints[i].SMT.resize(n_log_regs);

		
			for (int j = 0; j < n_log_regs; j++) {		

				branch_checkpoints[i].SMT[j] = 404;
	
			}

		//	for (int j = 0; j < n_log_regs; j++) {		
//
//				//printf("%d ", branch_checkpoints[i].SMT[j] );
//		
//			}
		
		//	printf("\n\n");
	
		}


		
		
		for (int i = 0 ; i<n_phys_regs ; i++) {
		
			PRFr[i] = 0;  // TA comment: set all true
			PRF[i] = 404; // like error 404 or something // so if any source is getting 404 value that means it has pull from a PR that had not been written to and likely source is wrongly renamed or something
		}

		for (int i = 0 ; i<n_log_regs ; i++) {
		
			RMT[i] = i;
			AMT[i] = i;
			PRFr[i] = 1;

		}

		//initialize free list -> issue with iq wakeup flag was each thing was getting assigned 0 dest and it was getting woken up again and again

		for (int i = 0 ; i < n_phys_regs - n_log_regs ;i++ ) {
		
			FL[i] = n_log_regs + i;
		
		}

		//active list is initially empty
		H_al = 0;
		H_al_phase = 0;
		T_al = 0;
		T_al_phase = 0;

		//free list is initially full
		H_fl = 0;
		H_fl_phase = 0;
		T_fl = 0;
		T_fl_phase = 1;

		GBM=0;

		


		


	};

	/////////////////////////////////////////////////////////////////////
	// This is the destructor, used to clean up memory space and
	// other things when simulation is done.
	// I typically don't use a destructor; you have the option to keep
	// this function empty.
	/////////////////////////////////////////////////////////////////////
	~renamer();


	//////////////////////////////////////////
	// Functions related to Rename Stage.   //
	//////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////
	// The Rename Stage must stall if there aren't enough free physical
	// registers available for renaming all logical destination registers
	// in the current rename bundle.
	//
	// Inputs:
	// 1. bundle_dst: number of logical destination registers in
	//    current rename bundle
	//
	// Return value:
	// Return "true" (stall) if there aren't enough free physical
	// registers to allocate to all of the logical destination registers
	// in the current rename bundle.
	/////////////////////////////////////////////////////////////////////
	bool stall_reg(uint64_t bundle_dst) {

		st_reg=0;

		if (H_fl == T_fl) {
			st_reg = (H_fl_phase == T_fl_phase) ? 1:0;
			return st_reg;
		}

			
		if (H_fl<T_fl) {
			st_reg =  ((T_fl - H_fl) < bundle_dst)?1:0;
			return st_reg;
		} else {
			st_reg = (((n_phys_regs - n_log_regs)-(H_fl - T_fl)) < bundle_dst)?1:0;
			return st_reg;
		}

	};

	/////////////////////////////////////////////////////////////////////
	// The Rename Stage must stall if there aren't enough free
	// checkpoints for all branches in the current rename bundle.
	//
	// Inputs:
	// 1. bundle_branch: number of branches in current rename bundle
	//
	// Return value:
	// Return "true" (stall) if there aren't enough free checkpoints
	// for all branches in the current rename bundle.
	/////////////////////////////////////////////////////////////////////
	bool stall_branch(uint64_t bundle_branch) {
	
 
		//number of ones in GBM should be greater than input yo continue
	
		st_br = 0;
		uint64_t bit_mask = 1;
    		uint64_t free = n_branches;
    		for (int i = 0 ; i<n_branches ;i++) {
       			 //printf("bit value = %d\n", (v>>i)&bit_mask);
        			if ((GBM>>i)&bit_mask)
           				 free--;
   		 }
		if (free < bundle_branch) 
			st_br = 1;

		return st_br;
	
	};

	/////////////////////////////////////////////////////////////////////
	// This function is used to get the branch mask for an instruction.
	/////////////////////////////////////////////////////////////////////
	uint64_t get_branch_mask() {
	
		return GBM;

	};

	/////////////////////////////////////////////////////////////////////
	// This function is used to rename a single source register.
	//
	// Inputs:
	// 1. log_reg: the logical register to rename
	//
	// Return value: physical register name
	/////////////////////////////////////////////////////////////////////
	uint64_t rename_rsrc(uint64_t log_reg) {
	
	// its coming from RMT
		//printf("\n\n	DEBUG: log_reg beging renamed: %d\n", log_reg);	
		//printf("	DEBUG: log_reg renamed to: %d\n\n", RMT[log_reg]);	
		return RMT[log_reg];

	}

	
	/////////////////////////////////////////////////////////////////////
	// This function is used to rename a single destination register.
	//
	// Inputs:
	// 1. log_reg: the logical register to rename
	//
	// Return value: physical register name
	/////////////////////////////////////////////////////////////////////
	uint64_t rename_rdst(uint64_t log_reg) {
	
		uint64_t PR = FL[H_fl];

				
		//sis comments
		//do i need to check if free list is empty ?

		H_fl++;
		if (H_fl == n_phys_regs - n_log_regs) {  // wrap around scenario
		     	H_fl = 0;              // force wrap around
	     		H_fl_phase = !H_fl_phase; // toggle phase bit
		}

	
		//SIS BUG: once you rename RMT must be updated as well ---> this line solved my bug of wrong source : it was dumb of mee -- sources come from RMT - i should have checked if RMT was correct at the first place
		RMT[log_reg] = PR;

		//sis exp - this line of code should be redundent
		//PRFr[PR] = 0;

		return PR;

	};

	/////////////////////////////////////////////////////////////////////
	// This function creates a new branch checkpoint.
	//
	// Inputs: none.
	//
	// Output:
	// 1. The function returns the branch's ID. When the branch resolves,
	//    its ID is passed back to the renamer via "resolve()" below.
	//
	// Tips:
	//
	// Allocating resources for the branch (a GBM bit and a checkpoint):
	// * Find a free bit -- i.e., a '0' bit -- in the GBM. Assert that
	//   a free bit exists: it is the user's responsibility to avoid
	//   a structural hazard by calling stall_branch() in advance.
	// * Set the bit to '1' since it is now in use by the new branch.
	// * The position of this bit in the GBM is the branch's ID.
	// * Use the branch checkpoint that corresponds to this bit.
	// 
	// The branch checkpoint should contain the following:
	// 1. Shadow Map Table (checkpointed Rename Map Table)
	// 2. checkpointed Free List head pointer and its phase bit
	// 3. checkpointed GBM
	/////////////////////////////////////////////////////////////////////
	uint64_t checkpoint( uint64_t vpq_tail_cp_i , bool vpq_phase_cp_i) {
		
		//sis-comment: when assigning GDB bit for check from bit to a branch start from 0 to n_branches
	
	       uint64_t bit_mask = 1;
    		uint64_t free = n_branches;

			uint64_t gbmt = GBM;

    		for (int i = 0 ; i<n_branches ;i++) {
        			if ((gbmt>>i)%2)
           				 free--;
   		 }
		
		assert(free != 0 && "There are no free bits available in GBM to assign to this branch\n");
		
		int free_bit = -1;	
	
		for (int i = 0 ; i<n_branches ;i++) {  //find free spot
		
			if ((gbmt>>i)%2 == 0) {
			
				free_bit = i;
				break; //free spot found 
			}
		}
	
		int branch_id = free_bit;
		//printf("Branch_ID: %d \n", branch_id);
		//branch_checkpoints[branch_id].SMT.resize(n_log_regs);
		for (uint64_t i = 0 ; i < n_log_regs ; i++) {
	//		printf("%d ", RMT[i]);
			branch_checkpoints[branch_id].SMT[i] = RMT[i];
		}
	
	//	printf("\n");	
			
		
		branch_checkpoints[branch_id].H_fl_cp = H_fl;
		branch_checkpoints[branch_id].H_fl_phase_cp = H_fl_phase;
	//	printf("	GBM before making checkpoint:  %d\n", gbmt);
		branch_checkpoints[branch_id].GBM_cp = gbmt;

		//SVP data
		branch_checkpoints[branch_id].vpq_tail_cp = vpq_tail_cp_i;
		branch_checkpoints[branch_id].vpq_phase_cp = vpq_phase_cp_i;

	//	printf("	GBM in checkpoint:  %d\n", branch_checkpoints[branch_id].GBM_cp);

		//setting GBM bit to reserve this branch's place
		//GBM = GBM + pow(2,branch_id);

		GBM |= (uint64_t) (1ULL << branch_id);
		
	//	 printf("	GBM after making checkpoint:  %d\n\n", GBM);

	//	 printf("-----------------------------\n\n");
	
		// branch_checkpoints[branch_id].svp_table_cp = svp_table_cp_i; 
	
			
		return branch_id;
	}; 

	//////////////////////////////////////////
	// Functions related to Dispatch Stage. //
	//////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////
	// The Dispatch Stage must stall if there are not enough free
	// entries in the Active List for all instructions in the current
	// dispatch bundle.
	//
	// Inputs:
	// 1. bundle_inst: number of instructions in current dispatch bundle
	//
	// Return value:
	// Return "true" (stall) if the Active List does not have enough
	// space for all instructions in the dispatch bundle.
	/////////////////////////////////////////////////////////////////////
	bool stall_dispatch(uint64_t bundle_inst){
	
		st_dp = 0 ;
		if (H_al == T_al) {
			st_dp = 1;
			return (H_al_phase == T_al_phase) ? 0:1;
		}

		if (H_al>T_al) {
			st_dp  = ( (H_al-T_al ) < bundle_inst)?1:0;
			return ( (H_al-T_al ) < bundle_inst)?1:0;
		} else {
			st_dp = ( (n_active - (T_al-H_al)) < bundle_inst)?1:0;
			return ( (n_active - (T_al-H_al)) < bundle_inst)?1:0;
		}
	

	};

	/////////////////////////////////////////////////////////////////////
	// This function dispatches a single instruction into the Active
	// List.
	//
	// Inputs:
	// 1. dest_valid: If 'true', the instr. has a destination register,
	//    otherwise it does not. If it does not, then the log_reg and
	//    phys_reg inputs should be ignored.
	// 2. log_reg: Logical register number of the instruction's
	//    destination.
	// 3. phys_reg: Physical register number of the instruction's
	//    destination.
	// 4. load: If 'true', the instr. is a load, otherwise it isn't.
	// 5. store: If 'true', the instr. is a store, otherwise it isn't.
	// 6. branch: If 'true', the instr. is a branch, otherwise it isn't.
	// 7. amo: If 'true', this is an atomic memory operation.
	// 8. csr: If 'true', this is a system instruction.
	// 9. PC: Program counter of the instruction.
	//
	// Return value:
	// Return the instruction's index in the Active List.
	//
	// Tips:
	//
	// Before dispatching the instruction into the Active List, assert
	// that the Active List isn't full: it is the user's responsibility
	// to avoid a structural hazard by calling stall_dispatch()
	// in advance.
	//
	//
	//	bool dest_f = 0;
	//		int log_reg = -1;
	//		int phys_reg = -1;
	//		bool completed_bit = 0;
	//		bool ex = 0;
	//		bool load_viol = 0;
	//		bool bmp = 0;
	//		bool vmp = 0;
	//		bool lf = 0;
	//		bool sf = 0;
	//		bool bf = 0;
	//		bool amof = 0;
	//		bool csr = 0;
	//		uint64_t pc = -1;
	//
	//
	//
	//
	/////////////////////////////////////////////////////////////////////
	uint64_t dispatch_inst(bool dest_valid,
	                       uint64_t log_reg,
	                       uint64_t phys_reg,
	                       bool load,
	                       bool store,
	                       bool branch,
	                       bool amo,
	                       bool csr,
	                       uint64_t PC) {
	
		
		// sis comment
		// do i need to check if AL is full ?
		assert(!((H_al==T_al)&&(H_al_phase != T_al_phase))&&"SIS ASSERT: Active list if full - cant dispatch instructions :/ \n");
		
		if (dest_valid) {
			
			AL[T_al].dest_f = dest_valid;
			AL[T_al].log_reg = log_reg;
			AL[T_al].phys_reg = phys_reg;

		} else {
		
			AL[T_al].dest_f = dest_valid;
			AL[T_al].log_reg = -1;
			AL[T_al].phys_reg = -1;
		}
		
		AL[T_al].lf = load;
		AL[T_al].completed_bit = false;
		AL[T_al].ex = false;
		AL[T_al].load_viol = false;
		AL[T_al].bmp = false;
		AL[T_al].vmp = false;
		AL[T_al].sf = store;
		AL[T_al].bf = branch;
		AL[T_al].amof = amo;
		AL[T_al].csr = csr;
		AL[T_al].pc = PC;

		uint64_t pointer = T_al;

		T_al++;
		if (T_al == n_active ) {  // wrap around scenario
		     	T_al = 0;              // force wrap around
	     		T_al_phase = !T_al_phase; // toggle phase bit
		}

		//printf("T_al status : %d\n", T_al);

		return pointer;

	};

	/////////////////////////////////////////////////////////////////////
	// Test the ready bit of the indicated physical register.
	// Returns 'true' if ready.
	/////////////////////////////////////////////////////////////////////
	bool is_ready(uint64_t phys_reg) {
	
		if (PRFr[phys_reg])
			return true;
		else
			return false;
	
	};

	/////////////////////////////////////////////////////////////////////
	// Clear the ready bit of the indicated physical register.
	/////////////////////////////////////////////////////////////////////
	void clear_ready(uint64_t phys_reg) {
	
		PRFr[phys_reg] = 0;

	};


	//////////////////////////////////////////
	// Functions related to the Reg. Read   //
	// and Execute Stages.                  //
	//////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////
	// Return the contents (value) of the indicated physical register.
	/////////////////////////////////////////////////////////////////////
	uint64_t read(uint64_t phys_reg){
	
		return PRF[phys_reg];
	
	};

	/////////////////////////////////////////////////////////////////////
	// Set the ready bit of the indicated physical register.
	/////////////////////////////////////////////////////////////////////
	void set_ready(uint64_t phys_reg) {
	
		PRFr[phys_reg] = 1;
	
	};


	//////////////////////////////////////////
	// Functions related to Writeback Stage.//
	//////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////
	// Write a value into the indicated physical register.
	/////////////////////////////////////////////////////////////////////
	void write(uint64_t phys_reg, uint64_t value){
	
		//printf("\n\n	DEBUG: overwriting PRF[%d] value's previous value of: %d with the value %d \n", phys_reg, PRF[phys_reg], value );
		PRF[phys_reg] = value;
		
	};	

	/////////////////////////////////////////////////////////////////////
	// Set the completed bit of the indicated entry in the Active List.
	/////////////////////////////////////////////////////////////////////
	void set_complete(uint64_t AL_index){
	
		AL[AL_index].completed_bit = 1;
	
	};

	/////////////////////////////////////////////////////////////////////
	// This function is for handling branch resolution.
	//
	// Inputs:
	// 1. AL_index: Index of the branch in the Active List.
	// 2. branch_ID: This uniquely identifies the branch and the
	//    checkpoint in question.  It was originally provided
	//    by the checkpoint function.
	// 3. correct: 'true' indicates the branch was correctly
	//    predicted, 'false' indicates it was mispredicted
	//    and recovery is required.
	//
	// Outputs: none.
	//
	// Tips:
	//
	// While recovery is not needed in the case of a correct branch,
	// some actions are still required with respect to the GBM and
	// all checkpointed GBMs:
	// * Remember to clear the branch's bit in the GBM.
	// * Remember to clear the branch's bit in all checkpointed GBMs.
	//
	// In the case of a misprediction:
	// * Restore the GBM from the branch's checkpoint. Also make sure the
	//   mispredicted branch's bit is cleared in the restored GBM,
	//   since it is now resolved and its bit and checkpoint are freed.
	// * You don't have to worry about explicitly freeing the GBM bits
	//   and checkpoints of branches that are after the mispredicted
	//   branch in program order. The mere act of restoring the GBM
	//   from the checkpoint achieves this feat.
	// * Restore the RMT using the branch's checkpoint.
	// * Restore the Free List head pointer and its phase bit,
	//   using the branch's checkpoint.
	// * Restore the Active List tail pointer and its phase bit
	//   corresponding to the entry after the branch's entry.
	//   Hints:
	//   You can infer the restored tail pointer from the branch's
	//   AL_index. You can infer the restored phase bit, using
	//   the phase bit of the Active List head pointer, where
	//   the restored Active List tail pointer is with respect to
	//   the Active List head pointer, and the knowledge that the
	//   Active List can't be empty at this moment (because the
	//   mispredicted branch is still in the Active List).
	// * Do NOT set the branch misprediction bit in the Active List.
	//   (Doing so would cause a second, full squash when the branch
	//   reaches the head of the Active List. We don’t want or need
	//   that because we immediately recover within this function.)
	/////////////////////////////////////////////////////////////////////
	void resolve(uint64_t AL_index,
		     uint64_t branch_ID,
		     bool correct) {

		// if branch is executed correctly
		// 	- clear the branch_ID bit of GBM
		// 	- clear the branch_ID bit in all of the checkpoints --> the bit in checkpoint just signifies to be prepared to be squashed if branch resolved wrong, but now as 
		// 	it has been resolved correct no need to care about that result, it can be used by another branch
		 
		// if branch was mispredicted
		// copy GBM    from     checkpoint ---> define these first so that it can be used and come back to code [ done ]
		// clear the branch_ID th bit from GBM
		// copy RMT    from      CP[branch_ID].SMT
		// copy H_fl   from      CP[branch_ID].F_fl_cp 
		// copy H_fl_phase   from      CP[branch_ID].F_fl_phase_cp 
		// T_al    =  AL_index + 1
		//
		// T AL phase
		// if T_al == H_al  --> full AL means opposite phases
		// if T_al > H_al  ---> same phase 
		// if H_al > T_al -----> opposite phases
		

		

		if (correct) {

			//printf("Branch was predicted correct \n");
		
			GBM = GBM - pow(2,branch_ID); //pow(2,branch_id)
			GBM &= ~(1ULL << branch_ID);

			//printf("NOW GBM = %d \n", GBM);
			
			//reference
			//	typedef struct BC{
			//		uint64_t *SMT;
			//		uint64_t H_fl_cp;
			//		uint64_t H_fl_phase_cp;
			//		uint64_t GBM_cp = ;
			//	}BC;
			//	BC *branch_checkpoints;


			uint64_t temp = 0 ;
		//	printf ("removing all branch dependencies \n");

			for (int i = 0 ; i < n_branches ; i++) {  //inspect each checkpoint

				temp = branch_checkpoints[i].GBM_cp  ;
				
				if ( (temp >> branch_ID)%2 ) { //check if bit corresponding to this branch is set in any of the checkpoints 
				
				//	branch_checkpoints[i].GBM_cp = branch_checkpoints[i].GBM_cp - pow(2,branch_ID); //if set then clear it
					//printf("	GBM in checkpoint before:  %d\n", branch_checkpoints[i].GBM_cp);
					branch_checkpoints[i].GBM_cp &= ~(1ULL << branch_ID);
				//	printf("	GBM in checkpoint:  %d\n", branch_checkpoints[i].GBM_cp);

				}		
			}
		} else {
		//	printf("	BRANCHHHHHH: Branch was predicted incorrect \n");
		
			GBM = (uint64_t) branch_checkpoints[branch_ID].GBM_cp;
		//	printf("NOW GBM = %d \n", GBM);
	

	
		//	movea2b_cp(RMT,branch_checkpoints[branch_ID].SMT );

			for (int i = 0 ; i < n_log_regs ; i++) {
			//	printf("%d ", RMT[i]);
				RMT[i] = branch_checkpoints[branch_ID].SMT[i];
			}	

			H_fl = branch_checkpoints[branch_ID].H_fl_cp;
			H_fl_phase = branch_checkpoints[branch_ID].H_fl_phase_cp;
		
			
		//	printf("AL_index given: %d \n", AL_index);
			if (AL_index == n_active - 1)
				T_al = 0;
			else
				T_al = AL_index + 1;

		//	printf("Recovered tail index : %d \n", T_al);

			if (H_al == T_al || H_al > T_al) {
				H_al_phase = !T_al_phase;
			} else {
				H_al_phase = T_al_phase;	
			}

		}
	
	}; //for phase 2

	void resolve2(uint64_t AL_index,
		     uint64_t branch_ID,
		     bool correct) {
	
		//if (correct) {

			//printf("Branch was predicted correct \n");
		
	//		uint64_t clear_mask = ~(1ULL << branch_ID);	
		//	GBM &= clear_mask;
			
		
	//		for (int i = 0 ; i < n_branches ; i++) {  //inspect each checkpoint
				 
				//	branch_checkpoints[i].GBM_cp &= clear_mask;

	//		}
		//}
	};

	//////////////////////////////////////////
	// Functions related to Retire Stage.   //
	//////////////////////////////////////////

	///////////////////////////////////////////////////////////////////
	// This function allows the caller to examine the instruction at the head
	// of the Active List.
	//
	// Input arguments: none.
	//
	// Return value:
	// * Return "true" if the Active List is NOT empty, i.e., there
	//   is an instruction at the head of the Active List.
	// * Return "false" if the Active List is empty, i.e., there is
	//   no instruction at the head of the Active List.
	//
	// Output arguments:
	// Simply return the following contents of the head entry of
	// the Active List.  These are don't-cares if the Active List
	// is empty (you may either return the contents of the head
	// entry anyway, or not set these at all).
	// * completed bit
	// * exception bit
	// * load violation bit
	// * branch misprediction bit
	// * value misprediction bit
	// * load flag (indicates whether or not the instr. is a load)
	// * store flag (indicates whether or not the instr. is a store)
	// * branch flag (indicates whether or not the instr. is a branch)
	// * amo flag (whether or not instr. is an atomic memory operation)
	// * csr flag (whether or not instr. is a system instruction)
	// * program counter of the instruction
	/////////////////////////////////////////////////////////////////////
	bool precommit(bool &completed,
                       bool &exception, bool &load_viol, bool &br_misp, bool &val_misp,
	               bool &load, bool &store, bool &branch, bool &amo, bool &csr,
		       uint64_t &PC) {

	
		/* bool dest_f = 0;
		int log_reg = -1;
		int phys_reg = -1;
		bool completed_bit = 0;
		bool ex = 0;
		bool load_viol = 0;
		bool bmp = 0;
		bool vmp = 0;
		bool lf = 0;
		bool sf = 0;
		bool bf = 0;
		bool amof = 0;
		bool csr = 0;
		uint64_t pc = -1; */



		
		if ( (H_al == T_al) && (H_al_phase == T_al_phase) )
			return false;
		else {
		
			completed = AL[H_al].completed_bit; 
			exception = AL[H_al].ex; 
			load_viol = AL[H_al].load_viol; 
			br_misp   = AL[H_al].bmp; 
			val_misp  = AL[H_al].vmp; 
			load      = AL[H_al].lf; 
			store     = AL[H_al].sf; 
			branch    = AL[H_al].bf; 
			amo       = AL[H_al].amof; 
			csr       = AL[H_al].csr; 
			PC        = AL[H_al].pc;


			return true;

		}
	
	} ;

	/////////////////////////////////////////////////////////////////////
	// This function commits the instruction at the head of the Active List.
	//
	// Tip (optional but helps catch bugs):
	// Before committing the head instruction, assert that it is valid to
	// do so (use assert() from standard library). Specifically, assert
	// that all of the following are true:
	// - there is a head instruction (the active list isn't empty)
	// - the head instruction is completed
	// - the head instruction is not marked as an exception
	// - the head instruction is not marked as a load violation
	// It is the caller's (pipeline's) duty to ensure that it is valid
	// to commit the head instruction BEFORE calling this function
	// (by examining the flags returned by "precommit()" above).
	// This is why you should assert() that it is valid to commit the
	// head instruction and otherwise cause the simulator to exit.
	/////////////////////////////////////////////////////////////////////
	uint64_t commit() {
	
		assert(AL[H_al].completed_bit && "SIS ASSERT: Active List H instruction not completed");
		assert(!AL[H_al].ex && "SIS ASSERT: Active List H instruction has an exception :''(   ");
		assert(!AL[H_al].load_viol && "SIS ASSERT: Active List H instruction has aload violation :''c   ");
		assert( !((H_al == T_al)&&(H_al_phase == T_al_phase)) && "SIS ASSERT: Active List is empty :/   ");

		//this commiting and freeing only happens when there is a valid destination


		if(AL[H_al].dest_f) {

			uint64_t tobefreed = AMT[AL[H_al].log_reg]; 
			FL[T_fl] = tobefreed; //freeing ----> tail points towards empty slot in FL


			//move FL tail forward and AL head forward 
			T_fl++;
			if (T_fl == n_phys_regs - n_log_regs) {
		
				T_fl = 0;
				T_fl_phase = !T_fl_phase;
			}

			AMT[AL[H_al].log_reg] = AL[H_al].phys_reg ; //commiting

			

		}

		H_al++;
		if (H_al == n_active) {
		
			H_al = 0;
			H_al_phase = !H_al_phase;
		}


		return AL[H_al].pc;

	
	};

	//////////////////////////////////////////////////////////////////////
	// Squash the renamer class.
	//
	// Squash all instructions in the Active List and think about which
	// sructures in your renamer class need to be restored, and how.
	//
	// After this function is called, the renamer should be rolled-back
	// to the committed state of the machine and all renamer state
	// should be consistent with an empty pipeline.
	/////////////////////////////////////////////////////////////////////
	void squash() {
	
		//Squashing steps:
		//1) AL squash : 
		//2) FL restore
		//3) Flash copy AMT to RMT
		//4) no branches in the pipeline anymore
		
		T_al = H_al; //bring tail back to head
		T_al_phase = H_al_phase;

		H_fl = T_fl;
		H_fl_phase = !T_fl_phase;

		for (int i = 0 ; i < n_phys_regs ;i++) {
		
			PRFr[i] = 0;

		}


		for(int i = 0 ; i < n_log_regs ; i++) {
		
			RMT[i] = AMT[i];
			PRFr[AMT[i]] = 1; //always comitted ones in PRF will be ready
		
		}

		GBM=0; // no branches anymore so this resets

		for (int i = 0 ; i < n_branches ; i++) {  //inspect each checkpoint
				
				
			branch_checkpoints[i].GBM_cp =  0;

			for (int j = 0 ; j < n_log_regs ; j++) {

				branch_checkpoints[i].SMT[j] = 0;

			}

				branch_checkpoints[i].H_fl_cp = 0;
				branch_checkpoints[i].H_fl_phase_cp = 0;

		}

			
	
	};

	//////////////////////////////////////////
	// Functions not tied to specific stage.//
	//////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////
	// Functions for individually setting the exception bit,
	// load violation bit, branch misprediction bit, and
	// value misprediction bit, of the indicated entry in the Active List.
	/////////////////////////////////////////////////////////////////////
	void set_exception(uint64_t AL_index) {
	
		AL[AL_index].ex = 1;

	};
	void set_load_violation(uint64_t AL_index){
	
		AL[AL_index].load_viol = 1;

	};
	void set_branch_misprediction(uint64_t AL_index) {
	
		AL[AL_index].bmp = 1;

	};
	void set_value_misprediction(uint64_t AL_index) {
	
		AL[AL_index].vmp = 1;

	};

	/////////////////////////////////////////////////////////////////////
	// Query the exception bit of the indicated entry in the Active List.
	/////////////////////////////////////////////////////////////////////
	bool get_exception(uint64_t AL_index) {
	
		return AL[AL_index].ex; 

	};


	//////////////////////////////////////////////////////////////////////
	// sis - Value Prediction
	//////////////////////////////////////////////////////////////////////
	
	
	//svp (0,0,0,PERF_VAL_PRED); 
	



};
